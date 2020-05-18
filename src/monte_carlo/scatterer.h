#ifndef scatterer_h
#define scatterer_h

#include <iostream>
#include <array>
#include <fstream>
#include <experimental/filesystem>
#include <regex>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <tuple>
#include <armadillo>

#include "../helper/utility.h"

#include "./scattering_struct.h"

namespace mc {

class scatterer {

private:
  typedef std::vector<std::vector<scattering_struct*>> scatt_t;
  typedef std::vector<std::vector<int>> map_t;

	double _max_rate; // maximum scattering rate in the scattering table
	double _inverse_max_rate; // inverse of the maximum scattering rate which is the lifetime
	arma::vec _pos; // position of the scatterer
  arma::vec _orientation; // orientation of the scattering object site
  arma::vec _chirality; // chirality of scattering site
  bool _is_quench = false;
  
public:
  // pointer to the scatterer on the right side (one side) of the current scatterer
  int right=-1; 

  // pointer to the scatterer on the left side (other side) of the current scatterer
  int left=-1;

  // index of neighboring grid cells
  std::vector<std::vector<scatterer*>*> close_scats; 

  // pointer to nearest quenching sites
  std::vector<std::vector<scatterer*>*> close_quenches; 

  // pointer to the scattering struct
  scatt_t scat_tab;
  //const scattering_struct* scat_tab = nullptr;

  map_t chirality_map;

public:

	// default constructor
  //scatterer(): _max_rate(0), _inverse_max_rate(0), right(-1), left(-1), scat_tab(nullptr) {};
  scatterer(): _max_rate(0), _inverse_max_rate(0), right(-1), left(-1) {};

  // set position of the scatterer
  void set_pos(const arma::vec& position) { _pos = position; };

  // set a component of the scatterer position
  void set_pos(const unsigned& i, const double& value) { _pos(i) = value; };

  // get position of the scatterer
  const arma::vec& pos() const { return _pos; };

  // get i'th component of position of the scatterer
  const double& pos(const unsigned& i) const { return _pos(i); };

  // set the orientation of the scatterer object
  void set_orientation(const arma::vec& m_orientation) { _orientation = m_orientation; };

  // set the i'th element of the orientation of the scatterer object
  void set_orientation(const mc::t_uint& i, const double& value) { _orientation(i) = value; };

  // get the orientation of the scatterer object
  const arma::vec& orientation() const { return _orientation; };

  // get the i'th element of the orientation of the scatterer object
  const double& orientation(const unsigned& i) const {
    return _orientation(i);
  };

  // set quenching sites
  void set_quenching() {_is_quench = true; };

  // set chirality of the scatterer
  void set_chirality(const arma::vec& chirality) { _chirality = chirality; };

  // get chirality of the scatterer
  const arma::vec& chirality() const { return _chirality; };

  // get random free flight time
  double ff_time() const {
    int r;
    while ((r = rand()) == 0) {
    }

    return -_inverse_max_rate * std::log(double(r) / double(RAND_MAX));
  };

  // find a new scattering object by looking at the scattering table
  const scatterer* update_state(const double& max_hopping_radius) const;

  // return the maximum scattering rate
  const double& max_rate() const { return _max_rate; };

  // set the max scattering rate by finding the neighbors
  // TODO: there might be a problem when there is nothing in neighbor list
  void set_max_rate(const double& max_hopping_radius){
    auto neighbors = find_neighbors(max_hopping_radius);
    if (neighbors.empty()){
      _max_rate = 5;
    } else {
      _max_rate = neighbors.back().first;
    }
    _inverse_max_rate = 1./_max_rate;
  };

  // find neighbors of the current scatterer and their scattering rates
  std::vector < std::pair<double, scatterer*>> find_neighbors(const double& max_hopping_radius) const;

  bool check_quenching(const double& max_disolving_radius) const;

  // count number of scatterer neighbors
  int no_of_neighbors(const double& max_hopping_radius) const;

  //return index of chirality using map
  int chiral_index(const arma::vec& chiral) const {
    int tube_size = size(chirality_map);
    for(int i = 0; i<tube_size; i++){
      if(chirality_map[i][0]==chiral[0] && chirality_map[i][1]==chiral[1])
        return i;
    }
    std::cout<< "panic: cannot find corresponding chirality!!"<< std::endl;
    return -1;
  }

};  // end class scatterer

} // end namespace mc

#endif  // scatterer_h
