#ifndef _exciton_transfer_h_
#define _exciton_transfer_h_

#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>
#include <type_traits>

#include "cnt.h"
#include "../helper/prepare_directory.hpp"

class exciton_transfer
{
private:
  typedef std::experimental::filesystem::directory_entry t_directory; // custom directory type
  t_directory _directory; // this is the address of the directory that the simulation data is stored

  std::string _name;
  double _temperature; // temperature of the system to calculate thermal distribution for excitons and other particles
  double _broadening_factor; // broadening factor  
  std::array<const cnt*,2> _cnts = {nullptr,nullptr}; // array of pointers to the target excitons

  enum simulation_mode {ex_trans_vs_angle, ex_trans_vs_zshift, ex_trans_vs_axis_shift_1, ex_trans_vs_axis_shift_2};
  simulation_mode _sim_mode;

  nlohmann::json _j_prop; // json data structure to hold properties of the simulation

  // function to return the lorentzian based on the broadening factor
  const double lorentzian(const double& energy)
  {
    return constants::inv_pi*_broadening_factor/(energy*energy + _broadening_factor*_broadening_factor);
  };

public:
  // constructor
  exciton_transfer(const cnt& cnt1, const cnt& cnt2) {
    _directory = prepare_directory("~/research/exciton_transfer",true);

    _cnts = {&cnt1, &cnt2};
    _temperature = 300;
    _broadening_factor = 4.e-3*constants::eV;

    std::cout << "\n...exciton transfer parameters:\n";
    std::cout << "temperature: " << _temperature << " [Kelvin]\n";
    std::cout << "energy broadening factor: " << _broadening_factor/constants::eV << " [eV]\n";
  };

  // constructor
  exciton_transfer(nlohmann::json j, const std::vector<cnt>& cnts, std::string parent_directory) {
    namespace fs = std::experimental::filesystem;
    
    // set the the initial and final cnt
    for (const auto& cnt: cnts){
      if (cnt.name() == j["cnt 1"]){
        _cnts[0] = &cnt;
        std::cout << "donor cnt set to " << cnt.name() << std::endl;
      }
    }
    for (const auto& cnt: cnts){
      if (cnt.name() == j["cnt 2"]){
        _cnts[1] = &cnt;
        std::cout << "acceptor cnt set to " << cnt.name() << std::endl;
      }
    }

    // set the name
    _name = _cnts[0]->name() + "_" + _cnts[1]->name();

    // prepare the directory
    fs::path directory_path(parent_directory);
    directory_path /= _name;
    bool keep_old_results = false;
    if (j.count("keep old results")==1){
      keep_old_results = j["keep old results"];
    }
    _directory = prepare_directory(directory_path, keep_old_results);

    // set tempareture
    if ((j["temperature"][1] != "Kelvin" and j["temperature"][1] != "kelvin")){
      std::cout << "temperature units:..." << j["temperature"][1] << "\n";
      throw std::invalid_argument("temperature must be in \"Kelvin\" units");
    }
    _temperature = j["temperature"][0];
    
    // set the broadening factor
    if (j["broadening factor"][1] != "meV"){
      throw std::invalid_argument("broadening factor must be in \"meV\" units");
    }
    _broadening_factor = double(j["broadening factor"][0])*1.e-3*constants::eV;

    std::cout << "\n...exciton transfer parameters:\n";
    std::cout << "temperature: " << _temperature << " [Kelvin]\n";
    std::cout << "energy broadening factor: " << _broadening_factor*1.e3/constants::eV << " [meV]\n";

    _j_prop = j;
  };

  // struct to bundle information about the excitonic states that are relevant
  struct ex_state {
    ex_state(const cnt::exciton_struct& m_exciton, const int& m_ik_cm_idx, const int& m_i_principal) {
      exciton = &m_exciton;
      cnt_obj = m_exciton.cnt_obj;
      elec_struct = m_exciton.elec_struct;
      ik_cm = m_ik_cm_idx+m_exciton.ik_cm_range[0];
      i_principal = m_i_principal;
      ik_cm_idx = m_ik_cm_idx;
      energy = m_exciton.energy(ik_cm_idx,i_principal);
      mu_cm=0;
    };

    const cnt::exciton_struct* exciton=nullptr; // reference to the exciton struct that owns the state
    const cnt* cnt_obj=nullptr; // reference to the cnt object owning the exciton state
    const cnt::el_energy_struct* elec_struct=nullptr; // reference to the el_energy_struct that is used to calculate the exciton dipersion
    int ik_cm_idx=0; // index of the ik_cm state in the exciton.energy and exciton.psi matrices
    int i_principal=0; // the principal quantum number of the state in the exciton.energy and exciton.psi matrices

    double energy=0; // energy of the exciton state
    int ik_cm=0; // value of the ik_cm for the state
    int mu_cm=0; // value of the mu for the state

    // access to the whole exciton state wavefunction
    arma::cx_vec psi() const {
      return exciton->psi.slice(ik_cm_idx).col(i_principal);
    };

    // access to individual elements of exciton state wavefunction
    std::complex<double> psi(int ik_c_idx) const {
      return exciton->psi(ik_c_idx,i_principal,ik_cm_idx);
    };

    const arma::umat& ik_idx() const {
      return exciton->ik_idx.slice(ik_cm_idx);
    }

    unsigned int ik_idx(const int& j, const int& i_n_principal) const {
      return exciton->ik_idx(j,i_n_principal, ik_cm_idx);
    }

    const arma::vec& dk_l() const {
      return *(exciton->dk_l);
    }

    const arma::vec K_cm() const {
      return ik_cm*(*(exciton->dk_l));
    } 

  };

  // struct to bundle information about initial and final states that match energetically
  struct matching_states {
    matching_states(const ex_state& d_state, const ex_state& a_state): i(d_state), f(a_state) {};
    const ex_state i; // initial exciton state
    const ex_state f; // final exciton state
  };

  // get the energetically relevant states in the form a vector of ex_state structs
  std::vector<ex_state> get_relevant_states(const cnt::exciton_struct& exciton, const double min_energy);

  // calculate Q()
  std::complex<double> calculate_Q(const matching_states& pair) const;

  // calculate J()
  std::complex<double> calculate_J(const matching_states& pair, const std::array<double,2>& shifts_along_axis, const double& z_shift, const double& angle) const;

  // match states based on energies
  std::vector<matching_states> match_states(const std::vector<ex_state>& d_relevant_states, const std::vector<ex_state>& a_relevant_states) {

    // lambda function to check if a state is energetically matched using a lorentzian
    auto is_matched = [this](const ex_state &d_state, const ex_state &a_state) {
      double delta_e = d_state.energy - a_state.energy;
      if (lorentzian(delta_e) > 1.e-2 * lorentzian(0))
        return true;
      return false;
    };

    std::vector<matching_states> matched;
    
    for (const auto& d_state: d_relevant_states) {
      for (const auto& a_state: a_relevant_states) {
        if (is_matched(d_state, a_state)) {
          matched.emplace_back(matching_states(d_state,a_state));
        }
      }
    }

    // std::cout << "\n...calculated pairs of donor and acceptor states\n";
    // std::cout << "number of pairs: " << matched.size() << " out of " << a_relevant_states.size()*d_relevant_states.size() << "\n\n";
    // for (const auto& pair:matched)
    // { 
    //   double delta_e = pair.i.energy-pair.f.energy;
    //   std::cout << "delta energy:" << delta_e/constants::eV << " [eV]   , lorentzian:" << lorentzian(delta_e)/lorentzian(0) << "\n";
    // }
    
    return matched;
  };

  // match states based on energies
  std::vector<matching_states> match_all_states(const std::vector<ex_state>& d_relevant_states, const std::vector<ex_state>& a_relevant_states) {
    std::vector<matching_states> matched;
    
    for (const auto& d_state: d_relevant_states) {
      for (const auto& a_state: a_relevant_states) {
        matched.emplace_back(matching_states(d_state,a_state));
      }
    }

    std::cout << "\n...matched all possible pairs of donor and acceptor states\n";
    std::cout << "number of pairs: " << matched.size() << "\n\n";
    
    return matched;
  };

  // calculate and plot Q matrix element between two exciton bands
  void save_Q_matrix_element(const int i_n_principal, const int f_n_principal);

  // calculate and plot J matrix element between two exciton bands
  void save_J_matrix_element(const int i_n_principal, const int f_n_principal);

  // calculate first order transfer rate
  double first_order(const double& z_shift, const std::array<double,2> axis_shifts, const double& theta, const bool& show_results=false);

  // calculate first order transfer rate for varying angle
  void calculate_first_order_vs_angle(const arma::vec& angle_vec ,const double& z_shift, const std::array<double,2> axis_shifts);

  // calculate first order transfer rate for center to center distance
  void calculate_first_order_vs_zshift(const arma::vec& z_shift_vec, const std::array<double,2> axis_shifts, const double& theta);

  // calculate first order transfer rate for varying axis shift for initial cnt
  void calculate_first_order_vs_axis_shift_1(const arma::vec& axis_shift_vec_1, const double axis_shift_2, const double z_shift, const double& theta);

  // calculate first order transfer rate for varying axis shift for final cnt
  void calculate_first_order_vs_axis_shift_2(const arma::vec& axis_shift_vec_2, const double axis_shift_1, const double z_shift, const double& theta);

  void run() {
    // if flag skip is set do not run this simulation
    if (_j_prop.count("skip")==1){
      if (_j_prop["skip"]) return;
    }

    // determine the execution policy
    if ((_j_prop["angle"].size()==4) and (_j_prop["zshift"].size()==2) and (_j_prop["axis shift 1"].size()==2) and (_j_prop["axis shift 2"].size()==2))
    {
      std::cout << "\nexciton transfer versus angle" << std::endl;

      // check distance units
      if ((_j_prop["zshift"][1] != "nm") or (_j_prop["axis shift 1"][1] != "nm") or (_j_prop["axis shift 2"][1] != "nm")){
        throw std::invalid_argument("distance units must be \"nm\"");
      }
      double z_shift = double(_j_prop["zshift"][0])*1.e-9;
      std::array<double,2> axis_shifts={double(_j_prop["axis shift 1"][0])*1.e-9, double(_j_prop["axis shift 2"][0])*1.e-9};

      // check angle units
      if (_j_prop["angle"][3]!= "degrees"){
        throw std::invalid_argument("angle units must be \"degrees\"");
      }
      arma::vec angle_vec = arma::linspace<arma::vec>(_j_prop["angle"][0],_j_prop["angle"][1],_j_prop["angle"][2])*constants::pi/180;

      calculate_first_order_vs_angle(angle_vec, z_shift, axis_shifts);
    }
    else if ((_j_prop["angle"].size()==2) and (_j_prop["zshift"].size()==4) and (_j_prop["axis shift 1"].size()==2) and (_j_prop["axis shift 2"].size()==2))
    {
      std::cout << "\nexciton transfer versus z_shift" << std::endl;

      // check distance units
      if ((_j_prop["zshift"][3] != "nm") or (_j_prop["axis shift 1"][1] != "nm") or (_j_prop["axis shift 2"][1] != "nm")){
        throw std::invalid_argument("distance units must be \"nm\"");
      }
      arma::vec z_shift_vec = arma::linspace<arma::vec>(_j_prop["zshift"][0],_j_prop["zshift"][1],_j_prop["zshift"][2])*1.e-9;
      std::array<double,2> axis_shifts={double(_j_prop["axis shift 1"][0])*1.e-9, double(_j_prop["axis shift 2"][0])*1.e-9};

      // check angle units
      if (_j_prop["angle"][3]!= "degrees"){
        throw std::invalid_argument("angle units must be \"degrees\"");
      }
      double angle = double(_j_prop["angle"][0])*constants::pi/180;
      calculate_first_order_vs_zshift(z_shift_vec, axis_shifts, angle);
    }
    else if ((_j_prop["angle"].size()==2) and (_j_prop["zshift"].size()==2) and (_j_prop["axis shift 1"].size()==4) and (_j_prop["axis shift 2"].size()==2))
    {
      std::cout << "\nexciton transfer versus axis shift 1" << std::endl;

      // check distance units
      if ((_j_prop["zshift"][1] != "nm") or (_j_prop["axis shift 1"][3] != "nm") or (_j_prop["axis shift 2"][1] != "nm")){
        throw std::invalid_argument("distance units must be \"nm\"");
      }
      arma::vec axis_shift_vec_1 = arma::linspace<arma::vec>(_j_prop["axis shift 1"][0],_j_prop["axis shift 1"][1],_j_prop["axis shift 1"][2])*1.e-9;
      double axis_shift_2 = double(_j_prop["axis shift 2"][0])*1.e-9;
      double z_shift = double(_j_prop["zshift"][0])*1.e-9;

      // check angle units
      if (_j_prop["angle"][3]!= "degrees"){
        throw std::invalid_argument("angle units must be \"degrees\"");
      }
      double angle = double(_j_prop["angle"][0])*constants::pi/180;
      calculate_first_order_vs_axis_shift_1(axis_shift_vec_1, axis_shift_2, z_shift, angle);
    }
    else if ((_j_prop["angle"].size()==2) and (_j_prop["zshift"].size()==2) and (_j_prop["axis shift 1"].size()==2) and (_j_prop["axis shift 2"].size()==3))
    {
      std::cout << "\nexciton transfer versus axis shift 2" << std::endl;

      // check distance units
      if ((_j_prop["zshift"][1] != "nm") or (_j_prop["axis shift 1"][1] != "nm") or (_j_prop["axis shift 2"][3] != "nm")){
        throw std::invalid_argument("distance units must be \"nm\"");
      }
      arma::vec axis_shift_vec_2 = arma::linspace<arma::vec>(_j_prop["axis shift 2"][0],_j_prop["axis shift 2"][1],_j_prop["axis shift 2"][2])*1.e-9;
      double axis_shift_1 = double(_j_prop["axis shift 1"][0])*1.e-9;
      double z_shift = double(_j_prop["zshift"][0])*1.e-9;

      // check angle units
      if (_j_prop["angle"][3]!= "degrees"){
        throw std::invalid_argument("angle units must be \"degrees\"");
      }
      double angle = double(_j_prop["angle"][0])*constants::pi/180;
      calculate_first_order_vs_axis_shift_2(axis_shift_vec_2, axis_shift_1, z_shift, angle);
    }
    else
    {
      throw std::logic_error("Invalid format for specifications of the exciton transfer simulation.");
    }
  };

  // save the coordinates of cnt atoms after rotation in 3d space
  typedef std::experimental::filesystem::path path_t;
  void save_atom_locations(path_t path, const std::array<double, 2> &shifts_along_axis, const double &z_shift, const double &angle, std::string prefix="");
};

#endif //_exciton_transfer_h_