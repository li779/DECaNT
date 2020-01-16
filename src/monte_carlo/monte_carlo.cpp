#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <omp.h>
#include <cassert>

#include "../../lib/json.hpp"
#include "../exciton_transfer/cnt.h"
#include "../helper/prepare_directory.hpp"
#include "../helper/progress.hpp"
#include "monte_carlo.h"


namespace mc
{

  // high level method to calculate proper scattering table
  std::vector<std::vector<scattering_struct>> monte_carlo::create_scattering_table(nlohmann::json j) {
    assert(j.count("rate type")>0);

    std::string rate_type = j["rate type"].get<std::string>();

    std::cout << "\ninitializing scattering table with " << rate_type << "..." << std::endl;


    if (rate_type == "davoody") {

      // get the parent directory for cnts
      std::string parent_directory = j["cnts"]["directory"];
      j["cnts"].erase("directory");
      j["cnts"].erase("comment");

      // create excitons and calculate exciton dispersions
      std::vector<cnt> cnts;
      cnts.reserve(j["cnts"].size()); // this is reservation of space is crucial to ensure we do not move
                                      // cnts, since the move constructor is not implemented yet
      
      for (const auto& j_cnt : j["cnts"]) {
        cnts.emplace_back(cnt(j_cnt, parent_directory));
        cnts.back().calculate_exciton_dispersion();
      };

	  std::vector<std::vector<scattering_struct>> all_tables(size(cnts));
	  for (int i = 0; i < size(cnts); i++) {
		  all_tables[i] = std::vector<scattering_struct>(size(cnts));
		  for (int j = 0; j < size(cnts); j++) {
			  all_tables[i][j] = create_davoody_scatt_table(cnts[i], cnts[j]);
		  }
	  }

	  return all_tables;
    }

    /*if (rate_type == "forster") {
      return create_forster_scatt_table(1.e15, 1.4e9);
    }

    if (rate_type == "wong") {
      return create_forster_scatt_table(1.e13, 1.4e9);
    }*/
    
    throw std::invalid_argument("rate type must be one of the following: \"davoody\", \"forster\", \"wong\"");

  };

  // method to calculate scattering rate via davoody et al. method
  scattering_struct monte_carlo::create_davoody_scatt_table(const cnt& d_cnt, const cnt& a_cnt) {
    auto zshift_prop = _json_prop["zshift [m]"];
    arma::vec z_shift = arma::linspace<arma::vec>(zshift_prop[0], zshift_prop[1], zshift_prop[2]);

    auto axis_shift_prop_1 = _json_prop["axis shift 1 [m]"];
    arma::vec axis_shift_1 = arma::linspace<arma::vec>(axis_shift_prop_1[0], axis_shift_prop_1[1], axis_shift_prop_1[2]);

    auto axis_shift_prop_2 = _json_prop["axis shift 2 [m]"];
    arma::vec axis_shift_2 = arma::linspace<arma::vec>(axis_shift_prop_2[0], axis_shift_prop_2[1], axis_shift_prop_2[2]);

    auto theta_prop = _json_prop["theta [degrees]"];
    arma::vec theta = arma::linspace<arma::vec>(theta_prop[0], theta_prop[1], theta_prop[2])*(constants::pi/180);

    arma::field<arma::cube> rate(theta.n_elem);
    rate.for_each([&](arma::cube& c){c.zeros(z_shift.n_elem, axis_shift_1.n_elem, axis_shift_2.n_elem);});

    exciton_transfer ex_transfer(d_cnt, a_cnt);

    ex_transfer.save_atom_locations(_output_directory.path(), {0, 0}, 1.5e-9, 0, ".0_angle");
    ex_transfer.save_atom_locations(_output_directory.path(), {0, 0}, 1.5e-9, constants::pi / 2, ".90_angle");
    ex_transfer.save_atom_locations(_output_directory.path(), {0, 0}, 1.5e-9, constants::pi, ".180_angle");


    #ifdef DEBUG_CHECK_RATES_SYMMETRY
    {
      double zsh = 1.5e-9;
      double ash1 = 0;
      double ash2 = 0;
      
      double th = 0;
      double r = ex_transfer.first_order(zsh, {ash1, ash2}, th, false);
      std::cout << "rate(" << th << ") = " << r << std::endl;

      th = constants::pi;
      r = ex_transfer.first_order(zsh, {ash1, ash2}, th, false);
      std::cout << "rate(" << th << ") = " << r << std::endl;

      std::exit(0);
    }
    #endif

    // progress_bar prog(theta.n_elem*z_shift.n_elem*axis_shift_1.n_elem*axis_shift_2.n_elem,"create davoody scattering table");
    progress_bar prog(theta.n_elem * z_shift.n_elem * axis_shift_1.n_elem * axis_shift_2.n_elem, "create davoody scattering table");

    #pragma omp parallel
    {
      double th, zsh, ash1, ash2;

      #pragma omp for
      for (unsigned i_th = 0; i_th<theta.n_elem; ++i_th) {

        th = theta(i_th);
        for (unsigned i_zsh = 0; i_zsh < z_shift.n_elem; ++i_zsh) {
          zsh = z_shift(i_zsh);
          for (unsigned i_ash1 = 0; i_ash1 < axis_shift_1.n_elem; ++i_ash1) {
            ash1 = axis_shift_1(i_ash1);
            for (unsigned i_ash2 = 0; i_ash2 < axis_shift_2.n_elem; ++i_ash2) {
              ash2 = axis_shift_2(i_ash2);
              // prog.step();
              rate(i_th)(i_zsh, i_ash1, i_ash2) = ex_transfer.first_order(zsh, {ash1, ash2}, th, false);
              
              #pragma omp critical
              {
                prog.step();
              }
            }
          }
        }

      }
    }

    scattering_struct scat_table(rate,theta,z_shift,axis_shift_1,axis_shift_2);

    double max_rate = 0;
    double min_rate = 10e15;
    rate.for_each([&min_rate](arma::cube& c) { min_rate = min_rate < c.min() ? min_rate : c.min(); });
    rate.for_each([&max_rate](arma::cube& c) { max_rate = max_rate > c.max() ? max_rate : c.max(); });
    
    std::cout << std::endl
              << "max rate in davoody scattering table: " << max_rate << " [1/s]" << std::endl
              << "min rate in davoody scattering table: " << min_rate << " [1/s]"
              << std::endl
              << std::endl;

    // std::string filename(_output_directory.path() / "davoody_scat_rates.dat");
    scat_table.save(_output_directory.path());

    return scat_table;
  };

  // method to calculate scattering rate via forster method
  scattering_struct monte_carlo::create_forster_scatt_table(double gamma_0, double r_0) {
    auto zshift_prop = _json_prop["zshift [m]"];
    arma::vec z_shift = arma::linspace<arma::vec>(zshift_prop[0], zshift_prop[1], zshift_prop[2]);

    auto axis_shift_prop_1 = _json_prop["axis shift 1 [m]"];
    arma::vec axis_shift_1 = arma::linspace<arma::vec>(axis_shift_prop_1[0], axis_shift_prop_1[1], axis_shift_prop_1[2]);

    auto axis_shift_prop_2 = _json_prop["axis shift 2 [m]"];
    arma::vec axis_shift_2 = arma::linspace<arma::vec>(axis_shift_prop_2[0], axis_shift_prop_2[1], axis_shift_prop_2[2]);

    auto theta_prop = _json_prop["theta [degrees]"];
    arma::vec theta = arma::linspace<arma::vec>(theta_prop[0], theta_prop[1], theta_prop[2])*(constants::pi/180);

    arma::field<arma::cube> rate(theta.n_elem);
    rate.for_each([&](arma::cube& c){c.zeros(z_shift.n_elem, axis_shift_1.n_elem, axis_shift_2.n_elem);});

    progress_bar prog(theta.n_elem*z_shift.n_elem*axis_shift_1.n_elem*axis_shift_2.n_elem,"create forster scattering table");

    unsigned i_th=0;
    for (const auto& th: theta) {
      unsigned i_zsh=0;
      for (const auto& zsh: z_shift) {
        unsigned i_ash1=0;
        for (const auto& ash1: axis_shift_1) {
          unsigned i_ash2=0;
          for (const auto& ash2: axis_shift_2) {
            prog.step();
            arma::vec r1 = {ash1, 0, 0};
            arma::vec r2 = {ash2*std::cos(th), ash2*std::sin(th), zsh};
            arma::vec dR = r1-r2;
            double angle_factor = std::cos(th)-3*arma::dot(arma::normalise(r1),arma::normalise(dR))*arma::dot(arma::normalise(r2),arma::normalise(dR));
            rate(i_th)(i_zsh,i_ash1,i_ash2) = gamma_0*std::pow(angle_factor,2)*std::pow(1.e-9/arma::norm(dR),6);
            i_ash2++;
          }
          i_ash1++;
        }
        i_zsh++;
      }
      i_th++;
    }

    scattering_struct scat_table(rate,theta,z_shift,axis_shift_1,axis_shift_2);

    return scat_table;
  };

  // slice the domain into n sections in each direction, and return a list of scatterers in the center region as the injection region
  std::vector<const scatterer *> monte_carlo::injection_region(const std::vector<scatterer> &all_scat, const domain_t domain, const int n) {
    assert((n > 0) && (n % 2 == 1));

    double xmin = domain.first(0), ymin = domain.first(1), zmin = domain.first(2);
    double xmax = domain.second(0), ymax = domain.second(1), zmax = domain.second(2);
    double dx = (xmax - xmin) / double(n), dy = (ymax - ymin) / double(n), dz = (zmax - zmin) / double(n);

    std::vector<double> x, y, z;

    for (int i = 0; i <= n; ++i) {
      x.push_back(double(i) * dx + xmin);
      y.push_back(double(i) * dy + ymin);
      z.push_back(double(i) * dz + zmin);
    }

    std::vector<const scatterer *> inject_list;

    for (const auto& s : all_scat) {
      if (x[n / 2] <= s.pos(0) && s.pos(0) <= x[n / 2 + 1] &&
          y[n / 2] <= s.pos(1) && s.pos(1) <= y[n / 2 + 1] &&
          z[n / 2] <= s.pos(2) && s.pos(2) <= z[n / 2 + 1])
        inject_list.push_back(&s);
    }


    return inject_list;
  }

  // slice the domain into n sections in each direction, and return the domain that leaves only 1 section from each side
  monte_carlo::domain_t monte_carlo::get_removal_domain(const monte_carlo::domain_t domain, const int n) {
    assert((n > 1));

    double xmin = domain.first(0), ymin = domain.first(1), zmin = domain.first(2);
    double xmax = domain.second(0), ymax = domain.second(1), zmax = domain.second(2);
    double dx = (xmax - xmin) / double(n), dy = (ymax - ymin) / double(n), dz = (zmax - zmin) / double(n);

    std::vector<double> x, y, z;

    for (int i = 0; i <= n; ++i) {
      x.push_back(double(i) * dx + xmin);
      y.push_back(double(i) * dy + ymin);
      z.push_back(double(i) * dz + zmin);
    }

    domain_t removal_domain;
    removal_domain.first = {x[1], y[1], z[1]};
    removal_domain.second = {x[n-1], y[n-1], z[n-1]};

    return removal_domain;
  }

  // initialize the simulation condition to calculate diffusion coefficient using green-kubo approach
  void monte_carlo::kubo_init() {
    // set maximum hopping radius
    _max_hopping_radius = double(_json_prop["max hopping radius [m]"]);
    std::cout << "maximum hopping radius: " << _max_hopping_radius * 1.e9 << " [nm]\n";

    _particle_velocity = _json_prop["exciton velocity [m/s]"];
    std::cout << "exciton velocity [m/s]: " << _particle_velocity << std::endl;

    _scat_tables = create_scattering_table(_json_prop);
    _all_scat_list = create_scatterers(_input_directory.path());

    domain_t d = find_simulation_domain();
    std::ios::fmtflags f(std::cout.flags()); // save cout flags to be reset after printing
    std::cout << std::fixed << std::showpos;
    std::cout << "\n"
              << "simulation domain BEFORE trimming:\n"
              << "    x (" << d.first(0) * 1e9 << " , " << d.second(0) * 1e9 << ") [nm]\n"
              << "    y (" << d.first(1) * 1e9 << " , " << d.second(1) * 1e9 << ") [nm]\n"
              << "    z (" << d.first(2) * 1e9 << " , " << d.second(2) * 1e9 << ") [nm]\n"
              << std::endl;
    std::cout.flags(f); // reset the cout flags

    limit_t xlim = _json_prop["trim limits"]["xlim"];
    limit_t ylim = _json_prop["trim limits"]["ylim"];
    limit_t zlim = _json_prop["trim limits"]["zlim"];

    trim_scats(xlim, ylim, zlim, _all_scat_list);

    _domain = find_simulation_domain();
    f = std::cout.flags(); // save cout flags to be reset after printing
    std::cout << std::fixed << std::showpos;
    std::cout << "\n"
              << "simulation domain AFTER trimming:\n"
              << "    x (" << _domain.first(0) * 1e9 << " , " << _domain.second(0) * 1e9 << ") [nm]\n"
              << "    y (" << _domain.first(1) * 1e9 << " , " << _domain.second(1) * 1e9 << ") [nm]\n"
              << "    z (" << _domain.first(2) * 1e9 << " , " << _domain.second(2) * 1e9 << ") [nm]\n"
              << std::endl;
    std::cout.flags(f); // reset the cout flags

    std::cout << "total number of scatterers: " << _all_scat_list.size() << std::endl;

    _quenching_list = create_quenching_sites(_all_scat_list, 1000);
    set_scat_table(_scat_tables[0][0], _all_scat_list);

    create_scatterer_buckets(_domain, _max_hopping_radius, _all_scat_list, _scat_buckets);
    set_max_rate(_max_hopping_radius, _all_scat_list);

    int n = _json_prop["number of sections for injection region"];
    _inject_scats = injection_region(_all_scat_list, _domain, n);
    _removal_domain = get_removal_domain(_domain, n);

    _max_time = _json_prop["maximum time for kubo simulation [seconds]"];
  };

  // create particles for kubo simulation
  void monte_carlo::kubo_create_particles() {
    int n_particle = _json_prop["number of particles for kubo simulation"];
    for (int i=0; i<n_particle; ++i) {
      int dice = std::rand() % _inject_scats.size();
      const scatterer *s = _inject_scats[dice];
      arma::vec pos = s->pos();
      _particle_list.push_back(particle(pos, s, _particle_velocity));
    }
  }

  // step the simulation in time
  void monte_carlo::kubo_step(double dt) {
    #pragma omp parallel
    {
      #pragma omp for
      for (unsigned i = 0; i < _particle_list.size(); ++i) {
        particle& p = _particle_list[i];
        
        p.step(dt, _all_scat_list, _max_hopping_radius);
        
        p.update_delta_pos();

        if (arma::any(p.pos()<_removal_domain.first) || arma::any(_removal_domain.second < p.pos())){
          int dice = std::rand() % _inject_scats.size();
          const scatterer* s = _inject_scats[dice];
          arma::vec pos = s->pos();
          p.set_pos(pos);
          p.set_scatterer(s);
        }
      }
    }

    // increase simulation time
    _time += dt;
  };

  // save the displacement of individual particles in kubo simulation
  void monte_carlo::kubo_save_individual_particle_dispalcements() {
    if (! _displacement_file_x.is_open()) {
      _displacement_file_x.open(_output_directory.path() / "particle_dispalcement.x.dat", std::ios::out);
      _displacement_file_y.open(_output_directory.path() / "particle_dispalcement.y.dat", std::ios::out);
      _displacement_file_z.open(_output_directory.path() / "particle_dispalcement.z.dat", std::ios::out);

      _displacement_file_x << std::showpos << std::scientific;
      _displacement_file_y << std::showpos << std::scientific;
      _displacement_file_z << std::showpos << std::scientific;

      _displacement_file_x << "time";
      _displacement_file_y << "time";
      _displacement_file_z << "time";
      for (int i=0; i<int(_particle_list.size()); ++i){
        _displacement_file_x << "," << i;
        _displacement_file_y << "," << i;
        _displacement_file_z << "," << i;
      }
      _displacement_file_x << std::endl;
      _displacement_file_y << std::endl;
      _displacement_file_z << std::endl;
    }

    _displacement_file_x << time();
    _displacement_file_y << time();
    _displacement_file_z << time();

    for (const auto& p: _particle_list) {
      _displacement_file_x << "," << p.delta_pos(0);
      _displacement_file_y << "," << p.delta_pos(1);
      _displacement_file_z << "," << p.delta_pos(2);
    }
    _displacement_file_x << std::endl;
    _displacement_file_y << std::endl;
    _displacement_file_z << std::endl;
  };

  void monte_carlo::kubo_save_avg_dispalcement_squared() {
    if (!_displacement_squard_file.is_open()) {
      _displacement_squard_file.open(_output_directory.path() / "particle_dispalcement.avg.squared.dat", std::ios::out);

      _displacement_squard_file << std::showpos << std::scientific;
      
      _displacement_squard_file << "# this file contains the average of dx^2, dy^2, and dz^2 of the particle ensemble over time" << std::endl
                                << "# number of particles: " << _particle_list.size() << std::endl
                                << std::endl;

      _displacement_squard_file << "time,x,y,z" << std::endl;
    }

    
    double avg_x2=0, avg_y2=0, avg_z2=0;

    for (const auto& p : _particle_list) {
      avg_x2 += std::pow(p.delta_pos(0), 2);
      avg_y2 += std::pow(p.delta_pos(1), 2);
      avg_z2 += std::pow(p.delta_pos(2), 2);
    }

    avg_x2 /= double(_particle_list.size());
    avg_y2 /= double(_particle_list.size());
    avg_z2 /= double(_particle_list.size());

    _displacement_squard_file << time() << "," << avg_x2 << "," << avg_y2 << "," << avg_z2 << std::endl;
  }

  void monte_carlo::kubo_save_diffusion_tensor(){
    if (!_diffusion_tensor_file.is_open()) {
      _diffusion_tensor_file.open(_output_directory.path() / "particle_diffusion_tensor.dat", std::ios::out);

     _diffusion_tensor_file << std::showpos << std::scientific;
      
      _diffusion_tensor_file << "# this file contains the diffusion tensor Dij of the particle ensemble over time" << std::endl
                                << "# number of particles: " << _particle_list.size() << std::endl
                                << std::endl;

      _diffusion_tensor_file << "time,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz" << std::endl;
    }

    _diffusion_tensor_file << time();

    double Dij, D1, D2, D3;
    for(int i = 0; i < 3; i++){
      for(int j = i; j < 3; j++){
        for(const auto& p : _particle_list){
          D1 += p.pos(i)*p.pos(j);
          D2 += p.pos(i);
          D3 += p.pos(j);
        }
        Dij = (D1 + D2 * D3)/(double)(_particle_list.size());
        Dij /= double(2 * double(time()));
        _diffusion_tensor_file << "," << Dij;
      }
    }
    _diffusion_tensor_file << std::endl;
  }

} // end of namespace mc