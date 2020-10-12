#ifndef particle_h
#define particle_h

#include <iostream>
#include <array>
#include <memory>
#include <armadillo>

#include "./free_flight.h"
#include "./scatterer.h"

namespace mc
{

class free_flight;

class particle
{

private:

  // pointer to scatter object for scattering the particle
  const scatterer* _scat_ptr=nullptr;

  // position of the particle
  arma::vec _pos;

  // position of the particle in the previous time step, this is used for boundary collision detection
  arma::vec _old_pos;

  // initial position of particle
  arma::vec _init_pos;

  // free flight time until next scattering event
  double _ff_time;

  // boolean determining the direction of the movement of particle in the cnt
  bool _heading_right;
  
  scatterer* _next_scat=nullptr;

  // velocity of the particle along the axis of the CNT
  double _velocity=0;

  // incremental displacement of the particle that is used for calculating diffusion coefficient through green-kubo method
  arma::vec _delta_pos{0,0,0};

  // final diffusion length
  arma::vec _diff_len{0,0,0};

  bool isDisolved = false;

public:
  particle() : _scat_ptr(nullptr), _pos({0, 0, 0}), _old_pos({0, 0, 0}), _ff_time(0), _heading_right(true), _velocity(0), _delta_pos({0,0,0}) {};

  particle(const arma::vec& pos, const scatterer* s, const double& velocity)
      : _scat_ptr(s), _pos(pos), _old_pos(pos), _velocity(velocity), _delta_pos({0,0,0}) {
    _ff_time = scat_ptr()->ff_time();
    _heading_right = std::rand()%2;
  };

  // perform free flight within the simulation domain
  void fly(double dt, const std::vector<scatterer>& s_list);

  // set the pointer to the scatterer object
  void set_scatterer(const scatterer* s) { _scat_ptr = s; };

  // return the pointer to the scatterer object
  const scatterer* scat_ptr() const { return _scat_ptr; };

  void set_disolved() {isDisolved = true; };

  const bool disolved() {return isDisolved; }

  // get position of the particle
  const arma::vec& pos() const { return _pos; };

  // get initial position of particle
  const arma::vec& init_pos() const { return _init_pos; };

  void set_init_pos(const arma::vec& pos) {_init_pos = pos; };

  // get position of the particle
  const double& pos(const double& i) const { return _pos(i); };

  // get old position of the particle
  const arma::vec& old_pos() const { return _old_pos; };

  // get old position of the particle
  const double& old_pos(const int& i) const { return _old_pos(i); };

  // set position of the particle and set the old position into _old_pos
  void set_pos(const arma::vec& pos) { _pos = pos; };

  // // set an element of particle position and set the old position into _old_pos
  void set_pos(const int& i, const double& value) { _pos(i) = value; };

  // // set the old position of the particle.
  void set_old_pos(const arma::vec& old_pos) { _old_pos = old_pos; };

  // return the free flight time until the next scattering
  const double& ff_time() const { return _ff_time; };

  // return the free flight time until the next scattering
  void set_ff_time(const double& value) { _ff_time = value; };

  // update the _ff_time by calling the underlying scatterer
  void update_ff_time() { _ff_time = _scat_ptr->ff_time(); };

  // step particle state for dt in time
  void step(double dt, const std::vector<scatterer>& s_list, const double& max_hop_radius);

  // update incremental displacement of the particle.
  void update_delta_pos() { _delta_pos += pos() - old_pos(); };

  // update diffusion length after trapped in a quenching site
  void update_diff_len() {_diff_len = pos() - init_pos();};

  const arma::vec& diff_len() {return _diff_len;};

  const double& diff_len(const int& i) const {return _diff_len(i);};

  // get the incremental dispalcement of the particle
  const arma::vec& delta_pos() const { return _delta_pos; };

  // get the i'th element of the incremental dispalcement of the particle
  const double& delta_pos(const int& i) const { return _delta_pos(i); };

}; //particle class

} //mc namespace

#endif // particle_h
