#ifndef monte_carlo_h
#define monte_carlo_h

#include <omp.h>
#include <algorithm>
#include <armadillo>
#include <cassert>
#include <chrono>
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <thread>

#include "../helper/utility.h"
#include "../helper/prepare_directory.hpp"
#include "../helper/constants.h"
#include "../helper/progress.hpp"

#include "../../lib/json.hpp"

#include "../exciton_transfer/cnt.h"
#include "../exciton_transfer/exciton_transfer.h"
#include "./particle.h"
#include "./scatterer.h"
#include "./scattering_struct.h"

namespace mc
{

class monte_carlo {

private:
  typedef std::experimental::filesystem::path path_t;
  typedef std::experimental::filesystem::directory_entry directory_t;
  typedef std::pair<arma::vec, arma::vec> domain_t;
  typedef std::vector<std::vector<scatterer*>> bucket_t;
  typedef std::vector<std::vector<scattering_struct>> scatt_t;
  typedef std::pair<double, double> limit_t;
  typedef std::vector<std::vector<int>> map_t;

  // elapsed simulation time
  double _time;

  // maximum hopping radius considered in the simulation
  double _max_hopping_radius;

  // maximum dissolving radius considered in the simulation
  double _max_dissolving_radius;

  // input properties of the whole mc simulation in json format
  nlohmann::json _json_prop;  

  // list of all scatterer object in the simulation
  std::vector<scatterer> _all_scat_list;

  // list of quenching sites in the simulation
  std::vector<scatterer> _quenching_list;

  // minimum and maximum coordinates of the simulation domain
  domain_t _domain;

  // number of particles in the contacts
  unsigned _c1_pop, _c2_pop;

  // this is the address of the output_directory and input_directory
  directory_t _output_directory, _input_directory, _scatter_table_directory;

  // scatter table directory
  std::string _scat_directory;

  // instantiation of the scattering table for discrete mesh points
  scatt_t _scat_tables;

  // pointers to scatterers to divide the scatterers into multiple buckets based on their position in space
  bucket_t _scat_buckets;  

  // pointers to quenching sites to divide the quenching sites into multiple buckets based on their position in space
  bucket_t _q_buckets;  

  // pointers to scatterers in contact 1 and 2
  std::vector<const scatterer*> _c1_scat, _c2_scat;

  // number of segments that defines contacts
  unsigned _n_seg=0;

  // area profile of the structure along the y-axis
  std::vector<double> _area;

  // list of all particles in the simulation
  std::vector<particle> _particle_list;

  // file objects for saving population profile and current data
  std::fstream _pop_file, _curr_file;

  // excitons' group velocity. Normally close to 0
  double _particle_velocity=0;

  // number of quenching sites
  int _quenching_sites_num = 0;

  // a mapping of index and chirality, used for constructing multi-dimension scatter tables
  map_t chirality_map;

  //**************************************************************
  // this section holds variables specific to the green-kubo method
  // of calculating diffusion coefficient
  //**************************************************************

  // list of scatterers in the injection region
  std::vector<const scatterer *> _inject_scats;

  // domain limits to remove the particles and inject them in the injection region
  domain_t _removal_domain;

  // maximum time for the kubo simulation
  double _max_time;

  // file to record particle dispalcements
  std::fstream _displacement_file_x, _displacement_file_y, _displacement_file_z;

  // file to record particle positions
  std::fstream _position_file_x, _position_file_y, _position_file_z;

  // file to record average of square of displacements in each direction
  std::fstream _displacement_squard_file;

  std::fstream _diffusion_tensor_file;

  std::fstream _diffusion_length_file;

  //**************************************************************
  //**************************************************************

 public:
  // default constructor
  monte_carlo() = delete;

  // constructure with json input file
  monte_carlo(const nlohmann::json& j) {

    std::cout << "\n"
              << "ready properties from json file"
              << "\n";

    // store the json properties for use in other methods
    _json_prop = j;

    // set the output directory
    std::string directory_path = j["output directory"];
    _scat_directory = j["scatter table directory"];
    bool        keep_old_data = true;
    if (j.count("keep old results") == 1) {
      keep_old_data = j["keep old results"];
    }
    _output_directory = prepare_directory(directory_path, keep_old_data);
    _scatter_table_directory = check_directory(_scat_directory, true);

    // set the input directory for mesh information
    directory_path = j["mesh input directory"];
    _input_directory = check_directory(directory_path, false);
  };

	// get the mc simulation time
  const double& time() const { return _time; };

  // get constant reference to the output_directory
  const directory_t& output_directory() const { return _output_directory; };

  // get constant reference to the input_directory
  const directory_t& input_directory() const { return _input_directory; };

  // get constant reference to the output_directory
  const path_t& output_path() const { return _output_directory.path(); };

  // get constant reference to the input_directory
  const path_t& input_path() const { return _input_directory.path(); };

  // returns the number of particles
  unsigned number_of_particles() const { return _particle_list.size(); };

  double calc_diam(int _m, int _n){
		double _a_cc = 1.42e-10; // carbon-carbon distance [nm]
 		double _a_l = std::sqrt(float(3.0))*_a_cc; // graphene lattice constants [nm]
		double _circum = _a_l*std::sqrt(float(_n*_n+_m*_m+_n*_m));
		double pi=3.141592;
 		return (_circum/pi);
  }

  // initialize the simulation condition
  void init() {
    _max_hopping_radius = double(_json_prop["max hopping radius [m]"]);
    std::cout << "maximum hopping radius: " << _max_hopping_radius * 1.e9 << " [nm]\n";

    _particle_velocity = _json_prop["exciton velocity [m/s]"];
    std::cout << "exciton velocity [m/s]: " << _particle_velocity << std::endl;

    _n_seg = _json_prop["number of segments"];
    std::cout << "number of segments: " << _n_seg << std::endl;

    _scat_tables = create_scattering_table(_json_prop);
    _all_scat_list = create_scatterers(_input_directory.path());

    limit_t xlim = _json_prop["trim limits"]["xlim"];
    limit_t ylim = _json_prop["trim limits"]["ylim"];
    limit_t zlim = _json_prop["trim limits"]["zlim"];

    trim_scats(xlim, ylim, zlim, _all_scat_list);
    
    std::cout << "total number of scatterers: " << _all_scat_list.size() << std::endl;

    set_scat_tables(_scat_tables,chirality_map, _all_scat_list);

    _domain = find_simulation_domain();

    _area = get_area(_n_seg);
    get_scatterer_statistics(_n_seg, _area);

    create_scatterer_buckets(_domain, _max_hopping_radius, _all_scat_list, _scat_buckets, _quenching_list, _q_buckets);
    set_max_rate(_max_hopping_radius, _all_scat_list);

    _c1_scat = contact_scats(_all_scat_list, _n_seg, 1, _domain);
    _c2_scat = contact_scats(_all_scat_list, _n_seg, _n_seg, _domain);

    _c1_pop = 1100;
    _c2_pop = 0;

    _particle_list = create_particles(_domain, _n_seg, _all_scat_list, _c1_pop, _c2_pop);
	};

  // read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage
  // particle hopping between the sites
  std::vector<scatterer> create_scatterers(const path_t& input_path){
    std::cout << std::endl << "create scatterers in fiber structure ... " << std::flush;

    std::ifstream pos_file;
    std::ifstream orient_file;
    std::ifstream chiral1_file;
    std::ifstream chiral2_file;

    // x axis
    pos_file.open(input_path / "single_cnt.pos.x.dat");
    orient_file.open(input_path / "single_cnt.orient.x.dat");

    arma::mat xcoor;
    xcoor.load(pos_file);
    xcoor *= 1.e-9;

    arma::mat xorient;
    xorient.load(orient_file);

    pos_file.close();
    orient_file.close();

    // y axis
    pos_file.open(input_path / "single_cnt.pos.y.dat");
    orient_file.open(input_path / "single_cnt.orient.y.dat");

    arma::mat ycoor;
    ycoor.load(pos_file);
    ycoor *= 1.e-9;

    arma::mat yorient;
    yorient.load(orient_file);

    pos_file.close();
    orient_file.close();

    // z axis
    pos_file.open(input_path / "single_cnt.pos.z.dat");
    orient_file.open(input_path / "single_cnt.orient.z.dat");

    arma::mat zcoor;
    zcoor.load(pos_file);
    zcoor *= 1.e-9;

    arma::mat zorient;
    zorient.load(orient_file);

    pos_file.close();
    orient_file.close();

    // chiral 1
    chiral1_file.open(input_path / "single_cnt.chiral.1.dat");
    chiral2_file.open(input_path / "single_cnt.chiral.2.dat");

    arma::mat chiral1;
    chiral1.load(chiral1_file);

    arma::mat chiral2;
    chiral2.load(chiral2_file);

    chiral1_file.close();
    chiral2_file.close();



    std::vector<scatterer> scat_list(xcoor.n_elem);

    for (unsigned i = 0; i < xcoor.n_rows; ++i) {
      for (unsigned j = 0; j < xcoor.n_cols; ++j) {
        unsigned n = i * xcoor.n_cols + j;

        scat_list[n].set_pos({xcoor(i, j), ycoor(i, j), zcoor(i, j)});
        scat_list[n].set_orientation({xorient(i, j), yorient(i, j), zorient(i, j)});
        scat_list[n].set_chirality({chiral1(i,j), chiral2(i,j)}); 
        if (j > 0) {
          scat_list[n].left = n - 1;
        }
        if (j + 1 < xcoor.n_cols) {
          scat_list[n].right = n + 1;
        }
      }
    }

    std::cout << "done!!!"
              << std::endl
              << std::endl
              << "total number of scatterers: " << scat_list.size()
              << std::endl;

    return scat_list;
  }

  // read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage
  // particle hopping between the sites
  std::vector<scatterer> create_quenching_sites(const std::vector<scatterer>& scat_list, int num_quenching){
    std::cout << std::endl << "create quenching sites in fiber structure ... " << std::flush;

    std::vector<scatterer> q_list(num_quenching);

    for (int n=0; n<num_quenching; n++){
      int dice = std::rand() % scat_list.size();
      const scatterer* s = &scat_list[dice];
      arma::vec pos = s->pos();
      arma::vec chirality = s->chirality();
      arma::vec orientation = s->orientation();
      double diameter = calc_diam(chirality[0],chirality[1]);

      arma::vec dia_vec = {diameter/2, 0, 0};
      arma::vec new_pos = pos + dia_vec;

      q_list[n].set_quenching();
      q_list[n].set_pos(new_pos);
    }

    std::cout << "done!!!"
              << std::endl
              << std::endl
              << "total number of quenching sites: " << q_list.size()
              << std::endl;

    return q_list;
  }

  // create particles with a linear density profile in y direction
  std::vector<particle> create_particles( const domain_t& domain, const unsigned n_seg,
      const std::vector<scatterer>& scat_list, int left_pop, int right_pop) {

    std::cout << "\n"
              << "create particles list:...";

    std::vector<particle> p_list;

    double y_min = domain.first(1);
    double y_max = domain.second(1);

    double dy = (y_max - y_min) / double(n_seg);
    double dp = double(right_pop - left_pop) / (double(n_seg) - 1);

    for (unsigned i=0; i<n_seg; ++i){


      int n_particle = std::round(left_pop + double(i) * dp);

      std::cout << "("<< i << "," << n_particle << ") ,";
      
      double y1 = y_min + double(i) * dy;
      double y2 = y1 + dy;

      std::vector<const scatterer*> s_list;
      for (const scatterer& s: scat_list){
        if (y1<=s.pos(1) && s.pos(1)<y2){
          s_list.emplace_back(&s);
        }
      }

      for (int n=0; n<n_particle; n++){
        int dice = std::rand()%s_list.size();
        const scatterer* s = s_list[dice];
        arma::vec pos = s->pos();
        p_list.push_back(particle(pos,s,_particle_velocity));
      }
    }

    std::cout << "...done!!!" << std::endl;

    return p_list;
  }

	// save the json properties that is read and parsed from the input_json file.
	void save_json_properties() {
		std::ofstream json_file;
		json_file.open(_output_directory.path() / "input.json", std::ios::out);
		json_file << std::setw(4) << _json_prop << std::endl;
		json_file.close();
	};

  // find minimum of the minimum coordinates of the scattering objects, this function will effectively give us the
  // simulation domain
  domain_t find_simulation_domain() const {
    arma::vec min_coor = _all_scat_list.front().pos();
    arma::vec max_coor = _all_scat_list.front().pos();

    for (const auto& s : _all_scat_list) {
      for (int i = 0; i < 3; ++i) {
        min_coor(i) = min_coor(i) > s.pos(i) ? s.pos(i) : min_coor(i);
        max_coor(i) = max_coor(i) < s.pos(i) ? s.pos(i) : max_coor(i);
      }
    }

    return {min_coor, max_coor};
  };

  // step the simulation in time
	void step(double dt) {

    # pragma omp parallel
    {
      #pragma omp for
      for (unsigned i=0; i<_particle_list.size(); ++i){
        _particle_list[i].step(dt, _all_scat_list, _max_hopping_radius, _max_dissolving_radius);
      }
    }

    // increase simulation time
    _time += dt;
	};

  // high level method to calculate proper scattering table
	scatt_t create_scattering_table(nlohmann::json j);

	// method to calculate scattering rate via forster method
	scattering_struct create_forster_scatt_table(double gamma_0, double r_0);

	// method to calculate scattering rate via davoody et al. method
	scattering_struct create_davoody_scatt_table(const cnt& d_cnt, const cnt& a_cnt);

  // divide scatterers into buckets based on their location, and set the pointers to enclosing and neighboring buckets
  // for each scatterer object
  void create_scatterer_buckets(const domain_t domain, const double radius, std::vector<scatterer>& scat_list,
                                bucket_t& scat_buckets, std::vector<scatterer>& q_list, bucket_t& q_buckets) {
    using namespace std;
    
    std::cout << "\n" 
              << "finding scatterer buckets: ";

    double xmin = (domain.first)(0);
    double xmax = (domain.second)(0);
    int nx = std::ceil((xmax - xmin) / radius) + 1;

    double ymin = (domain.first)(1);
    double ymax = (domain.second)(1);
    int    ny = std::ceil((ymax - ymin) / radius) + 1;

    double zmin = (domain.first)(2);
    double zmax = (domain.second)(2);
    int    nz = std::ceil((zmax - zmin) / radius) + 1;

    scat_buckets.resize(nx*ny*nz);
    q_buckets.resize(nx*ny*nz);

    for (scatterer& s : scat_list) {
      int ix = (s.pos(0) - xmin) / radius;
      int iy = (s.pos(1) - ymin) / radius;
      int iz = (s.pos(2) - zmin) / radius;
      int idx = ix + iy * nx + iz * nx * ny;
      scat_buckets[idx].push_back(&s);
    }

    for (scatterer& s : q_list) {
      int ix = (s.pos(0) - xmin) / radius;
      int iy = (s.pos(1) - ymin) / radius;
      int iz = (s.pos(2) - zmin) / radius;
      int idx = ix + iy * nx + iz * nx * ny;
      q_buckets[idx].push_back(&s);
    }


    for (scatterer& s : scat_list) {
      int ix = (s.pos(0) - xmin) / radius;
      int iy = (s.pos(1) - ymin) / radius;
      int iz = (s.pos(2) - zmin) / radius;

      for (int i : {ix - 1, ix, ix + 1}) {
        for (int j : {iy - 1, iy, iy + 1}) {
          for (int k : {iz - 1, iz, iz + 1}) {
            if (i > -1 && i < nx && j > -1 && j < ny && k > -1 && k < nz) {
              unsigned idx = i + j * nx + k * nx * ny;
              s.close_scats.push_back(&(scat_buckets[idx]));
              s.close_quenches.push_back(&(q_buckets[idx]));
            }
          }
        }
      }
    }

    std::cout << "done!\n";

  }

  // method to associate each scattering site to its scatter table.
  void set_scat_tables(scatt_t& _scat_tabs, map_t& _chirality_map, std::vector<scatterer>& scat_list) {
    int tube_size = size(_scat_tabs);
    for (auto& s : scat_list) {
      s.chirality_map = _chirality_map;
      s.scat_tab.resize(tube_size);
      for (int i = 0; i < tube_size; i++) {
		  s.scat_tab[i] = std::vector<scattering_struct*>(tube_size);
		  for (int j = 0; j < tube_size; j++) {
			  s.scat_tab[i][j] = &_scat_tabs[i][j];
		  }
	  }
    }
    std::cout <<"th tube's chirality: [" << _chirality_map[0][0] << ", " << _chirality_map[0][1] << "]" << std::endl;
  }

  // // set the pointer to scattering table struct for all scatterer objects
  // void set_scat_table(const scattering_struct& scat_tab, std::vector<scatterer>& scat_list) {
  //   for (auto& s : scat_list) {
  //     s.scat_tab = &scat_tab;
  //   }
  // }

  // set the max scattering rate for all the scatterers
  void set_max_rate(const double max_hopping_radius, std::vector<scatterer>& scat_list){
    
    progress_bar prog(scat_list.size(), "setting max rate in scatterers");

    # pragma omp parallel
    {
      # pragma omp for
      for (unsigned i=0; i<scat_list.size(); ++i) {
        scat_list[i].set_max_rate(max_hopping_radius);

        #pragma omp critical
        prog.step();
      }
    }
  }

  // repopulate contacts
  void repopulate_contacts() {
    double ymin = _domain.first(1);
    double ymax = _domain.second(1);
    double dy = (ymax - ymin) / double(_n_seg);

    double y1 = ymin;
    double y2 = ymin + dy;
    repopulate(y1, y2, _c1_pop, _c1_scat, _particle_list);

    y1 = ymin + double(_n_seg - 1) * dy;
    y2 = ymax;
    repopulate(y1, y2, _c2_pop, _c2_scat, _particle_list);
  };

  // take all the particles between ymin and ymax region and recycle them and populate the region with new particles
  void repopulate(const double ymin, const double ymax, const unsigned n_particle,
                  const std::vector<const scatterer*>& s_list,
                  std::vector<particle>& p_list) {
    
    unsigned j=p_list.size();

    for (unsigned i = 0; i < j;) {
      if (p_list[i].pos(1) >= ymin && p_list[i].pos(1) <= ymax) {
        --j;
        std::swap(p_list[i], p_list[j]);
      } else {
        ++i;
      }
    }

    int dice=0;
    unsigned n=0;
    unsigned final_size = j+n_particle;

    unsigned j_lim = std::min(int(p_list.size()), int(final_size));
    
    for (;j < j_lim; ++j) {
      dice = std::rand() % s_list.size();
      p_list[j] = particle(s_list[dice]->pos(), s_list[dice], _particle_velocity);
      ++n;
    }

    for (;n<n_particle; ++n){
      dice = std::rand() % s_list.size();
      p_list.emplace_back(particle(s_list[dice]->pos(), s_list[dice], _particle_velocity));
    }

    p_list.resize(final_size);
  }

  // create a list of scatterer pointers in the contact number i
  std::vector<const scatterer*> contact_scats(const std::vector<scatterer>& s_list, const unsigned n_seg, const unsigned i,
                                              const domain_t& domain) {
        
    assert(i>0);
    assert(i<=n_seg);
    
    double ymin = domain.first(1);
    double ymax = domain.second(1);
    double dy = (ymax-ymin)/double(n_seg);

    double y1 = ymin + double(i - 1) * dy;
    double y2 = ymin + double(i) * dy;

    std::vector<const scatterer*> c_list;

    for (auto& s: s_list){
      if (s.pos(1)>=y1 && s.pos(1)<=y2){
        c_list.push_back(&s);
      }
    }

    return c_list;
  }

  // calculate all the metrics needed from the experiment
  void save_metrics(double dt) {
    save_population_profile(_n_seg, dt);
    save_currents(_n_seg, dt);
  }

  // calculate and save population profile
  void save_population_profile(unsigned n, double dt) {
    assert(_area.size()==n);

    std::vector<int> pop(n, 0);

    double ymax = (_domain.second)(1);
    double ymin = (_domain.first)(1);
    double dy = (ymax - ymin) / double(n);

    if (!_pop_file.is_open()){
      _pop_file.open(_output_directory.path() / "population_profile.dat", std::ios::out);
      
      _pop_file << "area";
      for (unsigned i=0; i<_area.size(); ++i) {
        _pop_file << "," << std::scientific << std::showpos << _area[i];
      }
      _pop_file << std::endl
                << std::endl;

      _pop_file << "dy";
      for (unsigned i=0; i<_area.size(); ++i) {
        _pop_file << "," << std::scientific << std::showpos << dy;
      }
      _pop_file << std::endl
                << std::endl;

      _pop_file << "section pos";
      for (unsigned i=0; i<_area.size(); ++i) {
        _pop_file << "," << std::scientific << std::showpos << ymin+(double(i)+0.5)*dy;
      }
      _pop_file << std::endl
                << std::endl;

      _pop_file << "time";
      for (unsigned i = 0; i < _area.size(); ++i) {
        _pop_file << ",section" << i;
      }
      _pop_file << std::endl;
    }


    int i = 0;

    for (auto p : _particle_list) {
      i = (p.pos(1) - ymin) / dy;
      i = i < 0 ? 0 : (i < int(n) ? i : int(n) - 1);

      pop[i]++;
    }


    _pop_file << std::showpos << std::scientific << _time;
    for (unsigned j=0; j<pop.size(); ++j) {
      _pop_file << "," << double(pop[j])/(_area[j]*dy);
    }
    _pop_file << std::endl;
  }

  // calculate and save population profile
  void save_currents(int n, double dt) {
    assert(int(_area.size())==n);


    double ymax = _domain.second(1);
    double ymin = _domain.first(1);
    double dy = (ymax - ymin) / double(n);


    std::vector<double> y(n-1,0);
    std::vector<double> area_at_interface(n-1,0);
    for (int i=1; i<n; ++i){
      y[i-1] = ymin+dy*double(i);
      area_at_interface[i-1] = (_area[i-1]+_area[i])/2;
    }

    if (!_curr_file.is_open()){
      _curr_file.open(_output_directory.path() / "region_current.dat", std::ios::out);
      
      _curr_file << "interface area";
      for (const auto& a: area_at_interface){
        _curr_file << std::showpos << std::scientific << "," << a;
      }
      _curr_file << std::endl
                 << std::endl;

      _curr_file << "interface pos";
      for (const auto &yy : y)
      {
        _curr_file << std::showpos << std::scientific << "," << yy;
      }
      _curr_file << std::endl
                 << std::endl;

      _curr_file << "time";
      for (int i=1; i<n; ++i){
        _curr_file << ",interface" << (i-1);
      }
      _curr_file << std::endl;
    }


    std::vector<int> curr(n-1, 0);

    for (unsigned i = 0; i < y.size(); ++i) {
      for (auto p : _particle_list) {
        if (p.old_pos(1) < y[i] && p.pos(1) >= y[i]) {
          curr[i]++;
        } else if (p.old_pos(1) >= y[i] && p.pos(1) < y[i]) {
          curr[i]--;
        }
      }
    }    

    _curr_file << std::showpos << std::scientific << _time;
    for (unsigned i = 0; i<curr.size(); ++i) {
      _curr_file << std::showpos << std::scientific << "," <<  double(curr[i])/(area_at_interface[i]*dt);
    }
    _curr_file << std::endl;
  }

  // get the max area of the structure for n_seg segments along y-axis
  std::vector<double> get_area(unsigned n_seg){
    assert(n_seg > 0);

    double ymax = (_domain.second)(1);
    double ymin = (_domain.first)(1);
    double dy = (ymax - ymin) / double(n_seg);

    std::vector<double> xmax(n_seg,_domain.first(0));
    std::vector<double> xmin(n_seg,_domain.second(0));
    std::vector<double> zmax(n_seg,_domain.first(2));
    std::vector<double> zmin(n_seg,_domain.second(2));

    for (auto& s: _all_scat_list){
      int i = (s.pos(1)-ymin)/dy;
      i = i<0 ? 0 : (i<int(n_seg) ? i : n_seg-1);  // force i to be in range of 0 to n_seg-1

      if (xmin[i] > s.pos(0)) {
        xmin[i] = s.pos(0);
      } else if (xmax[i] < s.pos(0)) {
        xmax[i] = s.pos(0);
      }

      if (zmin[i] > s.pos(2)) {
        zmin[i] = s.pos(2);
      } else if (zmax[i] < s.pos(2)) {
        zmax[i] = s.pos(2);
      }
    }

    std::vector<double> area (n_seg, 0);
    for (unsigned i=0; i<n_seg; ++i){
      area[i] = (zmax[i] - zmin[i]) * (xmax[i] - xmin[i]);
    }

    // std::cout << "\nareas:";
    // for (unsigned i=0; i<n_seg; ++i){
    //   std::cout << i << " , " << area[i] << std::endl;
    // }

    // std::cin.ignore();

    return area;
  }

  // calculate and save statistics about all scatterer objects
  void get_scatterer_statistics(const unsigned n_seg, const std::vector<double>& area){
    assert(n_seg > 0);
    assert(n_seg == area.size());

    double ymin = _domain.first(1);
    double ymax = _domain.second(1);
    double dy = (ymax - ymin) / double(_n_seg);

    std::vector<long> pop(_n_seg, 0);

    for (auto& s:_all_scat_list){
      int i = int(std::abs(s.pos(1) - ymin) / dy)%_n_seg;
      pop[i]++;
    }

    std::vector<double> pos(_n_seg, 0);
    for (unsigned i=0; i<pos.size(); ++i){
      pos[i] = ymin + (double(i) + 0.5) * dy;
    }

    std::fstream f;
    f.open(_output_directory.path() / "scatterer_statistics.dat", std::ios::out);

    f << "position,distribution,population,density\n";
    for (unsigned i=0; i<pop.size(); ++i){
      f << std::scientific << pos[i] << "," << double(pop[i])/double(_all_scat_list.size()) << "," << pop[i] << "," << double(pop[i])/(area[i]*dy) << "\n";
    }
    f.close();
  }

  // trim all the scatterer objects outside a particular region.
  void trim_scats(const limit_t xlim, const limit_t ylim, const limit_t zlim,
                  std::vector<scatterer>& s_list) {
    
    std::cout << std::endl
              << "triming scattering list..."
              << std::flush;

    // swap two scatterer objects in scatterer_list and update the index of right and left scatterer objects
    auto swap_scatterers = [&s_list] (int i, int j){
      
      int iLeft = s_list[i].left;
      int iRight = s_list[i].right;

      int jLeft = s_list[j].left;
      int jRight = s_list[j].right;

      int new_i = j;
      int new_j = i;

      if (iLeft > -1) s_list[iLeft].right = new_i;
      if (iRight > -1) s_list[iRight].left = new_i;

      if (jLeft > -1) s_list[jLeft].right = new_j;
      if (jRight > -1) s_list[jRight].left = new_j;

      scatterer t = s_list[i];
      s_list[i] = s_list[j];
      s_list[j] = t;
    };

    int j = s_list.size();
    
    for (int i=0; i<j; ){

      if (s_list[i].pos(0) < xlim.first  || s_list[i].pos(1) < ylim.first  || s_list[i].pos(2) < zlim.first ||
          s_list[i].pos(0) > xlim.second || s_list[i].pos(1) > ylim.second || s_list[i].pos(2) > zlim.second) {
        --j;

        swap_scatterers(i,j);

        // delete the links to scatterer objects at location j.
        if (s_list[j].left > -1) s_list[s_list[j].left].right = -1;
        if (s_list[j].right > -1) s_list[s_list[j].right].left = -1;

      } else {
        ++i;
      }
    }

    s_list.resize(j);
    s_list.shrink_to_fit();

    unsigned count(0);
    for (unsigned i=0; i<s_list.size(); ++i){
      if (s_list[i].left == -1 && s_list[i].right == -1)
        count++;
    }

    std::cout << "...done!"
              << std::endl;
  }

  // this function, adds a particle from the left contact, tracks its position while it has not entered the right
  // contact and saves its position
  void track_particle(double dt, int fileNo){

    // assert(_particle_list.empty() && "particle list is not empty!");

    double ymin = _domain.first(1);
    double ymax = _domain.second(1);
    double dy = (ymax - ymin) / double(_n_seg);

    double y1 = ymin;
    double y2 = ymin + dy;
    
    unsigned c1_pop = 1;
    std::vector<particle> p_list;
    repopulate(y1, y2, c1_pop, _c1_scat, p_list);

    std::string base = _output_directory.path() / "particle_path.";
    std::stringstream filename;
    filename << base << fileNo << ".dat";
    std::ofstream file(filename.str().c_str(), std::ios::out);

    file << std::scientific << std::showpos << std::scientific;


    y1 = ymin + double(_n_seg - 1) * dy;
    y2 = ymax;
    while (p_list.front().pos(1)<y1){
      p_list.front().step(dt, _all_scat_list, _max_hopping_radius, _max_dissolving_radius);
      file << "   " << p_list.front().pos(0) << " " << p_list.front().pos(1) << " " << p_list.front().pos(2) << "\n";
    }
    file << std::endl;

    file.close();
  }

  // This method calculate mean square displacement of num_pop particles in the domain with respect to each time step.
  // It uses same methodology of track_particle that puts all partcles in the left contact and tracks them until it reaches right contact.
  // It outputs a file called mean_square_displacement.dat which records mean square displacement for every time step in the output folder
  // parameters are: dt - time step in second (used in step function)
  //                 num_pop - number of particles added initially to the domain
  void calc_diffusion(double dt, int num_pop) {

	  // assert(_particle_list.empty() && "particle list is not empty!");

    double ymin = _domain.first(1);
	  double ymax = _domain.second(1);
	  double dy = (ymax - ymin) / double(_n_seg);

	  double y1 = ymin;
	  double y2 = ymin + dy;

	  std::vector<particle> p_list;
	  repopulate(y1, y2, num_pop, _c1_scat, p_list);
	  
	  std::vector<double> orig_pos0;
	  std::vector<double> orig_pos1;
	  std::vector<double> orig_pos2;

	  for (unsigned i = 0; i < p_list.size();i++) {
		  orig_pos0.emplace_back(p_list[i].pos(0));
		  orig_pos1.emplace_back(p_list[i].pos(1));
		  orig_pos2.emplace_back(p_list[i].pos(2));
	  }


	  std::stringstream filename;
	  std::string base = _output_directory.path() / "mean_square_displacement.dat";
	  filename  << base;
	  std::ofstream file(filename.str().c_str(), std::ios::out);

	  file << std::scientific << std::showpos << std::scientific;

	  y1 = ymin + double(_n_seg - 1) * dy;
	  y2 = ymax;

	  unsigned num_left = p_list.size();
	  std::cout << num_left;
	  unsigned time = 0;
	  while (num_left > 0) {
		  double total_square_displace = 0;
		

		  for (unsigned i = 0; i < num_left;) {
			  if (p_list[i].pos(1) >= y1) {
				 /* std::swap(p_list[i], p_list[num_left-1]);
				  std::swap(orig_pos0[i], orig_pos0[num_left - 1]);
				  std::swap(orig_pos1[i], orig_pos1[num_left - 1]);
				  std::swap(orig_pos2[i], orig_pos2[num_left - 1]);*/
				  p_list.erase(p_list.begin()+i);
				  orig_pos0.erase(orig_pos0.begin()+i);
				  orig_pos1.erase(orig_pos1.begin()+i);
				  orig_pos2.erase(orig_pos2.begin()+i);
				  num_left=p_list.size();
			  }
			  else {
				  p_list[i].step(dt, _all_scat_list, _max_hopping_radius, _max_dissolving_radius);
				  double square_displace = (p_list[i].pos(0) - orig_pos0[i]) * (p_list[i].pos(0) - orig_pos0[i]) +
					  (p_list[i].pos(1) - orig_pos1[i]) * (p_list[i].pos(1) - orig_pos1[i]) +
					  (p_list[i].pos(2) - orig_pos2[i]) * (p_list[i].pos(2) - orig_pos2[i]);
				  total_square_displace = total_square_displace + square_displace;
				  i++;
			  }
		  }
		  time++;
		  std::cout << "\r" <<"Number of particle left: " <<num_left<<" Simulation time: "<< time;
		  file << "   " << double(total_square_displace/num_left) << "   \n";
	  }
	  
	  file << std::endl;

	  file.close();
	  std::cout << std::endl;
  }

  /*unsigned j = p_list.size();
  for (unsigned i = 0; i < j;) {
	  if (p_list[i].pos(1) >= ymin && p_list[i].pos(1) <= ymax) {
		  --j;
		  std::swap(p_list[i], p_list[j]);
	  }
	  else {
		  ++i;
	  }
  }
  */

  // initialize the simulation condition to calculate diffusion coefficient using green-kubo approach
  void kubo_init();

  // slice the domain into n sections in each direction, and return a list of scatterers in the center region as the injection region
  std::vector<const scatterer *> injection_region(const std::vector<scatterer>& all_scat, const domain_t domain, const int n);

  // slice the domain into n sections in each direction, and return the domain that leaves only 1 section from each side
  domain_t get_removal_domain(const domain_t domain, const int n);

  // create particles for kubo simulation
  void kubo_create_particles();

  // get maximum time for kubo simulation
  const double& kubo_max_time() const {return _max_time;};

  // step the kubo simulation in time
  void kubo_step(double dt);

  // save the displacement of individual particles in kubo simulation
  void kubo_save_individual_particle_dispalcements();

  void kubo_save_individual_particle_positions();

  // save the average displacement of particles in kubo simulation
  void kubo_save_avg_dispalcement_squared();

  void kubo_save_diffusion_tensor();

  void kubo_save_diffusion_length();

  bool check_scat_tab(std::experimental::filesystem::path path_ref);

  void print_exciton_scatter_times();

  scattering_struct recovery_scatt_table(std::experimental::filesystem::path path, const cnt& d_cnt, const cnt& a_cnt);

}; // end class monte_carlo

} // end namespace mc

#endif // monte_carlo_h