#ifndef scattering_struct_h
#define scattering_struct_h

#include <iostream>
#include <experimental/filesystem>
#include <armadillo>

namespace mc {

// struct to bundle information about scattering rate on discrete mesh points that are calculated through an arbitrary scattering mechanism
struct scattering_struct {
  // default constructor
  scattering_struct() {};

  // constructor
  scattering_struct(const arma::field<arma::cube>& m_rate,
                    const arma::vec& m_theta,
                    const arma::vec& m_z_shift,
                    const arma::vec& m_axis_shift_1,
                    const arma::vec& m_axis_shift_2) {
    rate = m_rate;
    theta = m_theta;
    z_shift = m_z_shift;
    axis_shift_1 = m_axis_shift_1;
    axis_shift_2 = m_axis_shift_2;
  };

  arma::field<arma::cube> rate;
  arma::vec theta;
  arma::vec z_shift;
  arma::vec axis_shift_1;
  arma::vec axis_shift_2;

  // get the scattering rate based on the precalculated rates for
  // discrete mesh points. Right now this equation only gets the
  // closes point on the mesh but I should implement some sort of
  // interpolation
  double get_rate(const double& m_theta, const double& m_z_shift,
                  const double& m_axis_shift_1,
                  const double& m_axis_shift_2) const {
    arma::vec tmp;
    tmp = arma::abs(theta - m_theta);
    unsigned i_th = tmp.index_min();
    tmp = arma::abs(z_shift - m_z_shift);
    unsigned i_zsh = tmp.index_min();
    tmp = arma::abs(axis_shift_1 - m_axis_shift_1);
    unsigned i_ash1 = tmp.index_min();
    tmp = arma::abs(axis_shift_2 - m_axis_shift_2);
    unsigned i_ash2 = tmp.index_min();

    return rate(i_th)(i_zsh, i_ash1, i_ash2);
  }

  // save the scattering table to file
  typedef std::experimental::filesystem::path path_t;
  void save(path_t path) {
    path /= "scat_table";

    std::ofstream file(std::string(path) + ".theta.dat", std::ios::ate);
    file << theta;
    file.close();

    file.open(std::string(path) + ".z_shift.dat", std::ios::ate);
    file << z_shift;
    file.close();

    file.open(std::string(path) + ".axis_shift_1.dat", std::ios::ate);
    file << axis_shift_1;
    file.close();

    file.open(std::string(path) + ".axis_shift_2.dat", std::ios::ate);
    file << axis_shift_2;
    file.close();

    file.open(std::string(path) + ".rates.dat", std::ios::ate);
    file << "sizes:" << std::endl;
    file << "theta, z_shift, axis_shift_1, axis_shift_2" << std::endl;
    file << theta.n_elem << ","<< z_shift.n_elem << "," << axis_shift_1.n_elem << "," << axis_shift_2.n_elem << std::endl;
    file << std::endl;
    

    for (unsigned i_th=0; i_th<theta.n_elem; i_th++) {
      for (unsigned i_zsh=0; i_zsh < z_shift.n_elem; i_zsh++) {
        for (unsigned i_ash1=0; i_ash1 < axis_shift_1.n_elem; i_ash1++) {
          for (unsigned i_ash2=0; i_ash2 < axis_shift_2.n_elem; i_ash2++) {
            file << rate(i_th)(i_zsh, i_ash1, i_ash2) << "\n";
          }
        }
      }
    }

    file.close();

  }
};

}
#endif // scattering_struct_h
