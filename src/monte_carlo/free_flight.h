#ifndef free_flight_h
#define free_flight_h

#include <iostream>
#include <array>

#include "../helper/utility.h"
#include "./particle.h"

namespace mc
{
class particle;

class free_flight {
private:
public:
  free_flight(){};  // constructor
  void fly(particle* p, const double& dt){};  // perform free_flight
  void check_boundary(particle* p, const double& dt,
                      const std::pair<arma::vec, arma::vec>& domain) const;  // check for collision to boundaries
}; // end class free_flight

} // end namespace mc

#endif  // free_flight_h
