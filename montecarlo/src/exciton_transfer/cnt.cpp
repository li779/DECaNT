/**
cnt.cpp
Stores all relevant information for a carbon nanotube
*/

#include <iostream>
#include <numeric>
#include <armadillo>
#include <complex>
#include <stdexcept>

#include "constants.h"
#include "cnt.h"
#include "../helper/progress.hpp"

void cnt::get_parameters()
{
  // graphen unit vectors and reciprocal lattice vectors
  _a1 = arma::vec({_a_l*std::sqrt(3.0)/2.0, +_a_l/2.0});
  _a2 = arma::vec({_a_l*std::sqrt(3.0)/2.0, -_a_l/2.0});
  _b1 = arma::vec({1.0/sqrt(3.0)*2.0*constants::pi/_a_l, +2.0*constants::pi/_a_l});
  _b2 = arma::vec({1.0/sqrt(3.0)*2.0*constants::pi/_a_l, -2.0*constants::pi/_a_l});

  // carbon-carbon translation vector
  _aCC_vec = 1.0/3.0*(_a1+_a2);

  // cnt chirality vector and its length
	_ch_vec = double(_n) * _a1 + double(_m) * _a2;
	_ch_len = arma::norm(_ch_vec,2);

  // cnt radius
  _radius = _ch_len/2.0/constants::pi;

  // calculate cnt t_vector
  int dR = std::gcd(2*_n+_m,_n+2*_m);
	_t1 = +(2*_m+_n)/dR;
	_t2 = -(2*_n+_m)/dR;
  _t_vec = double(_t1)*_a1 + double(_t2)*_a2;

  // number of hexagons in cnt unit cells which is equal to number of carbon atoms divided by half
  _Nu = 2*(std::pow(_n,2)+std::pow(_m,2)+_n*_m)/dR;


  // rotate basis vectors so that ch_vec is along the x_axis and t_vec is along y_axis
	double cos_theta = _ch_vec(0)/arma::norm(_ch_vec);
	double sin_theta = _ch_vec(1)/arma::norm(_ch_vec);
	arma::mat rot = {{+cos_theta, +sin_theta},
                   {-sin_theta, +cos_theta}}; // rotation matrix

	_ch_vec = rot*_ch_vec;
	_t_vec = rot*_t_vec;
	_a1 = rot*_a1;
	_a2 = rot*_a2;
	_b1 = rot*_b1;
	_b2 = rot*_b2;
	_aCC_vec = rot*_aCC_vec;

  	//make 3d t_vec where the cnt axis is parallel to y-axis
	_t_vec_3d = arma::vec(3,arma::fill::zeros);
	_t_vec_3d(1) = _t_vec(1);


  std::cout << "\n...graphene unit cell vectors:\n";
  _a1.print("a1:");
  _a2.print("a2:");

  std::cout << "\n...graphene reciprocal lattice vectors:\n";
  _b1.print("b1:");
  _b2.print("b2:");

  std::cout << "\n...vector connecting basis carbon atoms:\n";
  _aCC_vec.print("aCC vector:");

  _ch_vec.print("chirality vector:");
  std::cout << "ch_vec length:\n   " << _ch_len << std::endl;

  _t_vec.print("t_vec:");
  _t_vec_3d.print("3d t_vec:");


	// calculate reciprocal lattice of CNT
	_K1 = (-double(_t2)*_b1 + double(_t1)*_b2)/(double(_Nu));
	_K2 = (double(_m)*_b1-double(_n)*_b2)/(double(_Nu));
  _K2_normed = arma::normalise(_K2);
  _nk_K1 = _number_of_cnt_unit_cells;
	_dk_l = _K2/(double(_nk_K1));

  std::cout << "\n...cnt reciprocal lattice vectors:\n";
  _K1.print("K1:");
  _K2.print("K2:");

	// calculate K2-extended representation parameters
  {
    double p_min = (1./double(_t1)+1./double(_n))/(double(_m)/double(_n)-double(_t2)/double(_t1));
    double p_max = (1./double(_t1)+double(_Nu)/double(_n))/(double(_m)/double(_n)-double(_t2)/double(_t1));

    bool found = false;

    for (int p=std::ceil(p_min); p<std::ceil(p_max); p++)
    {
      if (((1+_t2*p) % _t1) == 0)
      {
        int q = (1+_t2*p)/_t1;
        _M = _m*p - _n*q;
        _Q = std::gcd(_Nu,_M);
        std::cout << "\n...K2-extended representation parameters:\n M: " << _M << " ,Q: " << _Q << "\n";
        found = true;
        break;
      }
    }
    if (not found)
    {
      std::cout << "Failed to calculate p and q for K2-extended representation .... investigate! .... aborting the simulation!!!\n";
    }
  }

}

// calculates position of atoms and reciprocal lattice vectors
void cnt::get_atom_coordinates()
{

	// calculate positions of atoms in the cnt unit cell
	_pos_a = arma::mat(_Nu,2,arma::fill::zeros);
	_pos_b = arma::mat(_Nu,2,arma::fill::zeros);

	int k = 0;

	for (int i=0; i<=_t1+_n; i++)
	{
		for (int j=_t2; j<=_m; j++)
		{
			bool flag1 = double(_t2*i)/(double)_t1 <= double(j);
			bool flag2 = double(_m*i)/(double)_n >= double(j);
			bool flag3 = double(_t2*(i-_n))/double(_t1) > double(j-_m);
			bool flag4 = double(_m*(i-_t1))/double(_n) < double(j-_t2);

			if(flag1 && flag2 && flag3 && flag4)
			{
        _pos_a.row(k) = double(i)*_a1.t() + double(j)*_a2.t();
        _pos_b.row(k) = _pos_a.row(k) + _aCC_vec.t();

				if(_pos_a(k,0) > _ch_vec(0))
          _pos_a(k,0) -= _ch_vec(0);
				if(_pos_a(k,0) < 0.0)
          _pos_a(k,0) += _ch_vec(0);
				if(_pos_a(k,1) > _ch_vec(1))
          _pos_a(k,1) -= _ch_vec(1);
				if(_pos_a(k,1) < 0.0)
          _pos_a(k,1) += _ch_vec(1);

				if(_pos_b(k,0) > _ch_vec(0))
          _pos_b(k,0) -= _ch_vec(0);
				if(_pos_b(k,0) < 0.0)
          _pos_b(k,0) += _ch_vec(0);
				if(_pos_b(k,1) > _ch_vec(1))
          _pos_b(k,1) -= _ch_vec(1);
				if(_pos_b(k,1) < 0.0)
          _pos_b(k,1) += _ch_vec(1);

				k++;
			}
		}
	}

  std::cout << "\n...atom coordinates:\n";
  _pos_a.print("pos_a:");
  _pos_b.print("pos_b:");

	if (k != _Nu)
	{
		std::cout << "error in finding position of atoms in cnt unit cell!!!" << std::endl;
		std::cout << "Nu = " << _Nu << "  ,  k = " << k << std::endl;
		exit(1);
	}

	// put position of all atoms in a single variable in 2d space(unrolled graphene sheet)
	_pos_2d = arma::mat(2*_Nu,2,arma::fill::zeros);
  _pos_2d.rows(0,_Nu-1) = _pos_a;
  _pos_2d.rows(_Nu,2*_Nu-1) = _pos_b;

	// calculate position of all atoms in the 3d space (rolled graphene sheet)
	_pos_3d = arma::mat(2*_Nu,3,arma::fill::zeros);
	for (unsigned int i=0; i<_pos_3d.n_rows; i++)
	{
		_pos_3d(i,0) = _radius*cos(_pos_2d(i,0)/_radius);
		_pos_3d(i,1) = _pos_2d(i,1);
		_pos_3d(i,2) = _radius*sin(_pos_2d(i,0)/_radius);
	}

	// save coordinates of atoms in 2d space
  std::string filename = _directory.path() / "pos_2d.dat";
  _pos_2d.save(filename, arma::arma_ascii);

  // save coordinates of atoms in 3d space
  filename = _directory.path() / "pos_3d.dat";
  _pos_3d.save(filename, arma::arma_ascii);

  // put position of all graphene unit cells in 2d (unrolled graphene sheet) and 3d space (rolled graphene sheet)
	_pos_u_2d = _pos_a;
	_pos_u_3d = arma::mat(_Nu,3,arma::fill::zeros);
	for (unsigned int i=0; i<_pos_u_3d.n_rows; i++)
	{
		_pos_u_3d(i,0) = _radius*cos(_pos_u_2d(i,0)/_radius);
		_pos_u_3d(i,1) = _pos_u_2d(i,1);
		_pos_u_3d(i,2) = _radius*sin(_pos_u_2d(i,0)/_radius);
	}

}

// calculate electron energy dispersions in the K1-extended representation using full unit cell (2*Nu atoms)
void cnt::electron_full_unit_cell()
{

	// make the list of 1st nearest neighbor atoms
	arma::umat nn_list(2*_Nu,3,arma::fill::zeros); // contains index of the nearest neighbor atom
	arma::imat nn_tvec_index(2*_Nu,3,arma::fill::zeros); // contains the index of the cnt unit cell that the nearest neigbor atom is in.
	for (unsigned int i=0; i<_pos_3d.n_rows; i++)
	{
		int k=0;
		for (unsigned int j=0; j<_pos_3d.n_rows; j++)
		{
			for (int l=-1; l<=1; l++)
			{
        double dR = arma::norm(_pos_3d.row(i)-_pos_3d.row(j)+double(l)*_t_vec_3d.t());
				if ( (i!=j) && (dR<(1.4*_a_cc)) )
				{
					nn_list(i,k) = j;
					nn_tvec_index(i,k) = l;
					k++;
				}
			}
		}
    if (k != 3)
    {
      std::cout << "error: nearest neighbors partially found!!!\n";
      std::exit(1);
    }
	}

	arma::cx_mat H(2*_Nu, 2*_Nu, arma::fill::zeros);
	arma::cx_mat S(2*_Nu, 2*_Nu, arma::fill::zeros);
  arma::vec E;
  arma::cx_mat C;

	int NK = _nk_K1;

	arma::mat el_energy_full(2*_Nu, NK, arma::fill::zeros);
	arma::cx_cube el_psi_full(2*_Nu, 2*_Nu, NK, arma::fill::zeros);

	double t_len = arma::norm(_t_vec_3d);

	for (int n=0; n<NK; n++)
	{
		double wave_vec = double(n-_nk_K1/2)*arma::norm(_dk_l);

		H.zeros();
		S.zeros();

		for (int i=0; i<2*_Nu; i++)
		{
			H(i,i) = std::complex<double>(_e2p,0.e0);
			for (int k=0; k<3; k++)
			{
				int j = nn_list(i,k);
				int l = nn_tvec_index(i,k);

				H(i,j) += arma::cx_double(_t0,0.e0)*exp(arma::cx_double(0.0,wave_vec*double(l)*t_len));
				S(i,j) += arma::cx_double(_s0,0.e0)*exp(arma::cx_double(0.0,wave_vec*double(l)*t_len));
			}
		}

		arma::eig_sym(E, C, H);

		// fix the phase of the eigen vectors
		for (unsigned int i=0; i<C.n_cols; i++)
		{
			arma::cx_double phi = std::conj(C(0,i))/std::abs(C(0,i));
			for (unsigned int j=0; j<C.n_rows; j++)
			{
				C(j,i) *= phi;
			}
		}

    el_energy_full.col(n) = E;
    el_psi_full.slice(n) = C;

	}

  // save electron energy bands using full Brillouine zone
  std::cout << "saved electron energy dispersion in K1-extended representation\n";
  std::string filename = _directory.path() / "el_energy_full.dat";
  el_energy_full.save(filename, arma::arma_ascii);

  // // save electron wavefunctions using full Brillouine zone
  // filename = _directory.path()/"el_psi_full.dat";
  // el_psi_full.save(filename, arma::arma_ascii);

}

// calculate electron dispersion energies for an input range of ik and mu
cnt::el_energy_struct cnt::electron_energy(const std::array<int,2>& ik_range, const std::array<int,2>& mu_range, const std::string& name)
{
  int number_of_bands = 2;
  int number_of_atoms_in_graphene_unit_cell = number_of_bands;

  int nk = ik_range[1] - ik_range[0];
  int n_mu = mu_range[1] - mu_range[0];


  arma::cube energy(number_of_bands, nk, n_mu, arma::fill::zeros);
  arma::field<arma::cx_cube> wavefunc(n_mu); // this pretty weird order is chosen so than we can select each cutting line easier and more efficiently
  wavefunc.for_each([&](arma::cx_cube& c){c.zeros(number_of_atoms_in_graphene_unit_cell, number_of_bands, nk);}); // this pretty weird order is chosen so than we can select each cutting line easier and more efficiently

  const std::complex<double> i1(0,1);
  
  const int ic = 1;
  const int iv = 0;
  
  const int iA = 0;
  const int iB = 1;

  for (int mu=mu_range[0]; mu<mu_range[1]; mu++)
  {
    for (int ik=ik_range[0]; ik<ik_range[1]; ik++)
    {
      arma::vec k_vec = double(mu)*_K1 + double(ik)*_dk_l;
      std::complex<double> fk = std::exp(std::complex<double>(0,arma::dot(k_vec,(_a1+_a2)/3.))) + \
                                std::exp(std::complex<double>(0,arma::dot(k_vec,(_a1-2.*_a2)/3.))) + \
                                std::exp(std::complex<double>(0,arma::dot(k_vec,(_a2-2.*_a1)/3.)));
      
      energy(ic,ik-ik_range[0],mu-mu_range[0]) = +_t0*std::abs(fk);
      energy(iv,ik-ik_range[0],mu-mu_range[0]) = -_t0*std::abs(fk);

      (wavefunc(mu-mu_range[0]))(iA,ic,ik-ik_range[0]) = +1./std::sqrt(2.); // this pretty weird order is chosen so than we can select each cutting line easier and more efficiently
      (wavefunc(mu-mu_range[0]))(iA,iv,ik-ik_range[0]) = +1./std::sqrt(2.); // this pretty weird order is chosen so than we can select each cutting line easier and more efficiently
      (wavefunc(mu-mu_range[0]))(iB,ic,ik-ik_range[0]) = -1./std::sqrt(2.)*std::conj(fk)/std::abs(fk); // this pretty weird order is chosen so than we can select each cutting line easier and more efficiently
      (wavefunc(mu-mu_range[0]))(iB,iv,ik-ik_range[0]) = +1./std::sqrt(2.)*std::conj(fk)/std::abs(fk); // this pretty weird order is chosen so than we can select each cutting line easier and more efficiently
    }
  }

  // save electron energy bands using full Brillouine zone
  std::string filename = _directory.path() / (name +".el_energy.dat");
  energy.save(filename,arma::arma_ascii);

  std::cout << "\n...calculated " + name + " electron dispersion\n";

  el_energy_struct energy_s;
  energy_s.name = name;
  energy_s.energy = energy;
  energy_s.wavefunc = wavefunc;
  energy_s.ik_range = ik_range;
  energy_s.mu_range = mu_range;
  energy_s.nk = nk;
  energy_s.n_mu = n_mu;
  energy_s.no_of_atoms = number_of_atoms_in_graphene_unit_cell;
  energy_s.no_of_bands = number_of_bands;

  return energy_s;
}

void cnt::find_valleys(const cnt::el_energy_struct& elec_struct)
{


  // get the indicies of valleys
  std::vector<std::array<unsigned int,2>> ik_valley_idx;

  int iC = 1;
  for (int ik_idx=0; ik_idx<elec_struct.nk; ik_idx++)
  {
    for (int i_mu_idx=0; i_mu_idx<elec_struct.n_mu; i_mu_idx++)
    {
      int ik_idx_p1 = ik_idx + 1;
      int ik_idx_m1 = ik_idx - 1;
      while(ik_idx_p1 >= elec_struct.nk)
      {
        ik_idx_p1 -= elec_struct.nk;
      }
      while(ik_idx_m1 < 0)
      {
        ik_idx_m1 += elec_struct.nk;
      }

      if ((elec_struct.energy(iC,ik_idx,i_mu_idx) < elec_struct.energy(iC,ik_idx_m1,i_mu_idx)) \
            and (elec_struct.energy(iC,ik_idx,i_mu_idx) < elec_struct.energy(iC,ik_idx_p1,i_mu_idx)))
      {
        ik_valley_idx.push_back({(unsigned int)ik_idx, (unsigned int)i_mu_idx});
      }

    }
  }

  // sort valleys in order of their
  std::sort(ik_valley_idx.begin(),ik_valley_idx.end(), [&](const auto& s1, const auto& s2) {
                                    return elec_struct.energy(1,s1[0],s1[1]) < elec_struct.energy(1,s2[0],s2[1]);
                                  });

  // put them in a vector with each element containing the two equivalent valleys
  for (unsigned int i=0; i<ik_valley_idx.size()/2; i++)
  {
    std::array<std::array<unsigned int, 2>, 2> valley = {ik_valley_idx.at(2*i), ik_valley_idx.at(2*i+1)};
    _valleys_K2.push_back(valley);
  }

  std::cout << "\n...found and sorted indices of valleys:\n";
  for (const auto& valleys: _valleys_K2)
  {
    auto v1 = valleys.at(0);
    auto v2 = valleys.at(1);
    std::cout << "[" << v1[0] << "," << v1[1] << "] , [" << v2[0] << "," << v2[1] << "]\n";
  }

  std::cout << "number of valleys: " << _valleys_K2.size() << std::endl;

}

// find ik values that are energetically relevant around the bottom of the valley
void cnt::find_relev_ik_range(double delta_energy, const cnt::el_energy_struct& elec_struct)
{
  std::vector<std::vector<std::array<int,2>>> relev_ik_range(2);

  // first get ik for relevant states in the first valley
  int i_valley = 0;
  int ik_bottom = _valleys_K2[_i_sub][i_valley][0]+elec_struct.ik_range[0];
  int mu_bottom = _valleys_K2[_i_sub][i_valley][1]+elec_struct.mu_range[0];
  int iC = 1;
  double max_energy = elec_struct.energy(iC, ik_bottom-elec_struct.ik_range[0], mu_bottom-elec_struct.mu_range[0]) + delta_energy;

  relev_ik_range.at(i_valley).push_back({ik_bottom,mu_bottom});
  bool in_range = true;
  int count = 0;
  while (in_range)
  {
    in_range = false;
    count ++;
    int ik = ik_bottom + count;
    while(ik >= elec_struct.ik_range[1])
    {
      ik -= elec_struct.nk;
    }
    if (elec_struct.energy(iC, ik-elec_struct.ik_range[0], mu_bottom-elec_struct.mu_range[0]) < max_energy)
    {
      relev_ik_range.at(i_valley).push_back({ik,mu_bottom});
      in_range = true;
    }

    ik = ik_bottom - count;
    while(ik < elec_struct.ik_range[0])
    {
      ik += elec_struct.nk;
    }
    if (elec_struct.energy(iC, ik-elec_struct.ik_range[0], mu_bottom-elec_struct.mu_range[0]) < max_energy)
    {
      // std::array<int,2> relev_state = {ik,mu_bottom};
      relev_ik_range.at(i_valley).insert(relev_ik_range.at(i_valley).begin(), {ik,mu_bottom});
      in_range = true;
    }
  }

  // do the same thing for the second valley
  i_valley = 1;
  ik_bottom = _valleys_K2[_i_sub][i_valley][0]+elec_struct.ik_range[0];
  mu_bottom = _valleys_K2[_i_sub][i_valley][1]+elec_struct.mu_range[0];
  max_energy = elec_struct.energy(iC, ik_bottom-elec_struct.ik_range[0], mu_bottom-elec_struct.mu_range[0]) + delta_energy;

  relev_ik_range.at(i_valley).push_back({ik_bottom,mu_bottom});
  in_range = true;
  count = 0;
  while (in_range)
  {
    in_range = false;
    count ++;
    int ik = ik_bottom + count;
    while(ik >= elec_struct.ik_range[1])
    {
      ik -= elec_struct.nk;
    }
    if (elec_struct.energy(iC, ik-elec_struct.ik_range[0], mu_bottom-elec_struct.mu_range[0]) < max_energy)
    {
      relev_ik_range.at(i_valley).push_back({ik,mu_bottom});
      in_range = true;
    }

    ik = ik_bottom - count;
    while(ik < elec_struct.ik_range[0])
    {
      ik += elec_struct.nk;
    }
    if (elec_struct.energy(iC, ik-elec_struct.ik_range[0], mu_bottom-elec_struct.mu_range[0]) < max_energy)
    {
      relev_ik_range.at(i_valley).insert(relev_ik_range.at(i_valley).begin(), {ik,mu_bottom});
      in_range = true;
    }
  }


  std::cout << "\n...ik for relevant states calculated:\n";
  std::cout << "relev_ik_range has length of " << relev_ik_range[0].size() << std::endl;

  // i_valley = 0;
  // std::cout << "valley: " << i_valley << "\n";
  // for (const auto& state: relev_ik_range.at(i_valley))
  // {
  //   std::cout << "   [" << state.at(0) << "," << state.at(1) << "]\n";
  // }

  // i_valley = 1;
  // std::cout << "\nvalley: " << i_valley << "\n";
  // for (const auto& state: relev_ik_range.at(i_valley))
  // {
  //   std::cout << "   [" << state.at(0) << "," << state.at(1) << "]\n";
  // }

  _relev_ik_range = relev_ik_range;

}

// fourier transformation of the coulomb interaction a.k.a v(q)
cnt::vq_struct cnt::calculate_vq(const std::array<int,2> iq_range, const std::array<int,2> mu_range, unsigned int no_of_cnt_unit_cells)
{
  // primary checks for function input
  int nq = iq_range.at(1) - iq_range.at(0);
  if (nq <= 0) {
    throw "Incorrect range for iq!";
  }
  int n_mu = mu_range.at(1) - mu_range.at(0);
  if (n_mu <= 0) {
    throw "Incorrect range for mu_q!";
  }
  if (no_of_cnt_unit_cells % 2 == 0)  no_of_cnt_unit_cells ++;

  // calculate distances between atoms in a warped cnt unit cell.
	arma::mat pos_aa = arma::mat(_Nu,2,arma::fill::zeros);
	arma::mat pos_ab = arma::mat(_Nu,2,arma::fill::zeros);
	arma::mat pos_ba = arma::mat(_Nu,2,arma::fill::zeros);
  arma::mat pos_bb = arma::mat(_Nu,2,arma::fill::zeros);

	for (int i=0; i<_Nu; i++)
	{
    pos_aa.row(i) = _pos_a.row(i)-_pos_a.row(0);
    pos_ab.row(i) = _pos_a.row(i)-_pos_b.row(0);
    pos_ba.row(i) = _pos_b.row(i)-_pos_a.row(0);
    pos_bb.row(i) = _pos_b.row(i)-_pos_b.row(0);

		if(pos_aa(i,0) > _ch_vec(0)/2)
      pos_aa(i,0) -= _ch_vec(0);
		if(pos_ab(i,0) > _ch_vec(0)/2)
      pos_ab(i,0) -= _ch_vec(0);
		if(pos_ba(i,0) > _ch_vec(0)/2)
      pos_ba(i,0) -= _ch_vec(0);
		if(pos_bb(i,0) > _ch_vec(0)/2)
      pos_bb(i,0) -= _ch_vec(0);
	}

  arma::cube rel_pos(_Nu*no_of_cnt_unit_cells,2,4,arma::fill::zeros);
  for (int i=-std::floor(double(no_of_cnt_unit_cells)/2.); i<=std::floor(double(no_of_cnt_unit_cells)/2.); i++)
  {
    int idx = (i+std::floor(double(no_of_cnt_unit_cells)/2.))*_Nu;
    for (int j=0; j<_Nu; j++)
    {
      rel_pos.slice(0).row(idx+j) = pos_aa.row(j)+(i*_t_vec.t());
      rel_pos.slice(1).row(idx+j) = pos_ab.row(j)+(i*_t_vec.t());
      rel_pos.slice(2).row(idx+j) = pos_ba.row(j)+(i*_t_vec.t());
      rel_pos.slice(3).row(idx+j) = pos_bb.row(j)+(i*_t_vec.t());
    }
  }

  // calculate vq
  arma::cx_cube vq(nq,n_mu,4,arma::fill::zeros);
  arma::vec q_vec(nq,arma::fill::zeros);

  arma::vec q(2,arma::fill::zeros);
  const double coeff = std::pow(4.*constants::pi*constants::eps0*_Upp/constants::q0/constants::q0,2);
  const std::complex<double> i1(0.,1.);
  auto Uhno = [&](const arma::mat& R){
    return std::exp(i1*arma::dot(q,R))*_Upp/std::sqrt(coeff*(std::pow(R(0),2)+std::pow(R(1),2))+1);
  };

  progress_bar prog(nq, "vq");

  for (int iq=iq_range[0]; iq<iq_range[1]; iq++)
  {
    int iq_idx = iq-iq_range[0];

    prog.step();

    q_vec(iq_idx) = iq*arma::norm(_dk_l,2);
    for (int mu=mu_range[0]; mu<mu_range[1]; mu++)
    {
      int mu_idx = mu - mu_range[0];
      q = iq*_dk_l + mu*_K1;
      // std::cout << "after addition!\n";
      for (int i=0; i<4; i++)
      {
        for (unsigned int k=0; k<_Nu*no_of_cnt_unit_cells; k++)
        {
          vq(iq_idx,mu_idx,i) += Uhno(rel_pos.slice(i).row(k));
        }
      }
    }
  }

  vq = vq/(2*_Nu*no_of_cnt_unit_cells);

  std::cout << "\n...calculated vq\n";

  std::cout << "saved real part of vq\n";
  arma::cube vq_real = arma::real(vq);
  std::string filename = _directory.path()/"vq_real.dat";
  vq_real.save(filename, arma::arma_ascii);

  std::cout << "saved imaginary part of vq\n";
  arma::cube vq_imag = arma::imag(vq);
  filename = _directory.path()/"vq_imag.dat";
  vq_imag.save(filename, arma::arma_ascii);

  std::cout << "saved q_vector for vq\n";
  filename = _directory.path()/"vq_q_vec.dat";
  q_vec.save(filename, arma::arma_ascii);

  // make the vq_struct that is to be returned
  vq_struct vq_s;
  vq_s.data = vq;
  vq_s.iq_range = iq_range;
  vq_s.mu_range = mu_range;
  vq_s.nq = nq;
  vq_s.n_mu = n_mu;

  return vq_s;
}

// polarization of electronic states a.k.a PI(q)
cnt::PI_struct cnt::calculate_polarization(const std::array<int,2> iq_range, const std::array<int,2> mu_range, const cnt::el_energy_struct& elec_struct)
{
  // primary checks for function input
  int nq = iq_range.at(1) - iq_range.at(0);
  if (nq <= 0) {
    throw "Incorrect range for iq in calculate_polarization!";
  }
  int n_mu = mu_range.at(1) - mu_range.at(0);
  if (n_mu <= 0) {
    throw "Incorrect range for mu_q in calculate_polarization!";
  }

  int ikq, mu_kq;
  int ik, mu_k;
  int iq, mu_q;
  // lambda function to wrap iq+ik and mu_k+mu_q inside the K2-extended brillouine zone
  auto get_kq = [&](){
    mu_kq = mu_k+mu_q;
    ikq = ik+iq;
    while (mu_kq >= elec_struct.mu_range[1]) {
      mu_kq -= elec_struct.n_mu;
      ikq += _nk_K1*_M;
    }
    while (mu_kq < elec_struct.mu_range[0]) {
      mu_kq += elec_struct.n_mu;
      ikq -= _nk_K1*_M;
    }
    while (ikq >= elec_struct.ik_range[1]){
      ikq -= elec_struct.nk;
    }
    while (ikq < elec_struct.ik_range[0]){
      ikq += elec_struct.nk;
    }
  };

  arma::mat PI(nq,n_mu,arma::fill::zeros);
  arma::vec q_vec(nq,arma::fill::zeros);

  const int iv = 0;
  const int ic = 1;

  progress_bar prog(nq, "calculate polarization");

  int iq_idx, mu_q_idx;
  int ik_idx, mu_k_idx;
  int i_kq_idx, mu_kq_idx;
  for (iq=iq_range[0]; iq<iq_range[1]; iq++)
  {
    iq_idx = iq - iq_range[0];
    q_vec(iq_idx) = iq*arma::norm(_dk_l);

    prog.step(iq_idx);

    for (mu_q=mu_range[0]; mu_q<mu_range[1]; mu_q++)
    {
      mu_q_idx = mu_q - mu_range[0];
      for (ik=elec_struct.ik_range[0]; ik<elec_struct.ik_range[1]; ik++)
      {
        ik_idx = ik - elec_struct.ik_range[0];
        for (mu_k=elec_struct.mu_range[0]; mu_k<elec_struct.mu_range[1]; mu_k++)
        {
          mu_k_idx = mu_k - elec_struct.mu_range[0];
          get_kq();
          mu_kq_idx = mu_kq - elec_struct.mu_range[0];
          i_kq_idx = ikq - elec_struct.ik_range[0];

          PI(iq_idx,mu_q_idx) += std::pow(std::abs(arma::dot(arma::conj(elec_struct.wavefunc(mu_k_idx).slice(ik_idx).col(iv)),\
                                                             elec_struct.wavefunc(mu_kq_idx).slice(i_kq_idx).col(ic))),2)/ \
                                          (elec_struct.energy(ic,i_kq_idx,mu_kq_idx)-elec_struct.energy(iv,ik_idx,mu_k_idx)) + \
                                 std::pow(std::abs(arma::dot(arma::conj(elec_struct.wavefunc(mu_k_idx).slice(ik_idx).col(ic)), \
                                                             elec_struct.wavefunc(mu_kq_idx).slice(i_kq_idx).col(iv))),2)/ \
                                          (elec_struct.energy(ic,ik_idx,mu_k_idx)-elec_struct.energy(iv,i_kq_idx,mu_kq_idx));
        }
      }
    }
  }

  PI = 2*PI;

  std::cout << "\n...calculated polarization: PI(q)\n";

  std::cout << "saved PI\n";
  std::string filename = _directory.path()/"PI.dat";
  PI.save(filename, arma::arma_ascii);

  std::cout << "saved q_vector for PI\n";
  filename = _directory.path()/"PI_q_vec.dat";
  q_vec.save(filename, arma::arma_ascii);

  // make the vq_struct that is to be returned
  PI_struct PI_s;
  PI_s.data = PI;
  PI_s.iq_range = iq_range;
  PI_s.mu_range = mu_range;
  PI_s.nq = nq;
  PI_s.n_mu = n_mu;

  return PI_s;
}

// dielectric function a.k.a eps(q)
cnt::epsilon_struct cnt::calculate_dielectric(const std::array<int,2> iq_range, const std::array<int,2> mu_range)
{
  // check if vq has been calculated properly before
  if (not (in_range(iq_range,_vq.iq_range) and in_range(mu_range,_vq.mu_range))){
    throw std::logic_error("You need to calculate vq with correct range before \
                            trying to calculate dielectric function");
  }
  // check if PI has been calculated properly before
  if (not (in_range(iq_range,_PI.iq_range) and in_range(mu_range,_PI.mu_range))){
    throw std::logic_error("You need to calculate PI with correct range before \
                            trying to calculate dielectric function");
  }

  int nq = iq_range[1] - iq_range[0];
  int n_mu = mu_range[1] - mu_range[0];
  arma::mat eps = arma::real(arma::mean(_vq.data,2));
  eps = eps.submat(iq_range[0]-_vq.iq_range[0],mu_range[0]-_vq.mu_range[0],arma::size(nq,n_mu));
  eps %= _PI.data.submat(iq_range[0]-_PI.iq_range[0],mu_range[0]-_PI.mu_range[0],arma::size(nq,n_mu));
  std::cout << "size of dielectric function matrix: " << arma::size(eps) << std::endl;
  eps += 1.;

  arma::vec q_vec(nq);
  for (int iq=iq_range[0]; iq<iq_range[1]; iq++)
  {
    int iq_idx = iq - iq_range[0];
    q_vec(iq_idx) = iq*arma::norm(_dk_l);
  }

  std::cout << "\n...calculated dielectric function: epsilon(q)\n";

  std::cout << "saved epsilon\n";
  std::string filename = _directory.path()/"eps.dat";
  eps.save(filename, arma::arma_ascii);

  std::cout << "saved q_vector for epsilon\n";
  filename = _directory.path()/"eps_q_vec.dat";
  q_vec.save(filename, arma::arma_ascii);

  epsilon_struct eps_s;
  eps_s.data = eps;
  eps_s.iq_range = iq_range;
  eps_s.mu_range = mu_range;
  eps_s.nq = nq;
  eps_s.n_mu = n_mu;
  return eps_s;
}

// calculate exciton dispersion
std::vector<cnt::exciton_struct> cnt::calculate_A_excitons(const std::array<int,2> ik_cm_range, const cnt::el_energy_struct& elec_struct) {
  // some utility variables that are going to be used over and over again
  int ik_c, mu_c;
  int ik_v, mu_v;
  int ik_cp, mu_cp;
  int ik_vp, mu_vp;
  int ik_c_diff, mu_c_diff;
  int ik_cm, mu_cm;

  std::complex<double> dir_interaction;
  std::complex<double> xch_interaction;

  const int iv = 0;
  const int ic = 1;

  const int i_valley_1 = 0;
  const int i_valley_2 = 1;

  // lambda function to calculate direct interaction
  auto get_direct_interaction = [&](){
    ik_c_diff = ik_c-ik_cp;
    mu_c_diff = mu_c-mu_cp;
    while(ik_c_diff < elec_struct.ik_range[0]){
      ik_c_diff += elec_struct.nk;
    }
    while(ik_c_diff >= elec_struct.ik_range[1]){
      ik_c_diff -= elec_struct.nk;
    }

    dir_interaction = 0;
    for (int i=0; i<2; i++)
    {
      for (int j=0; j<2; j++)
      {
        dir_interaction += std::conj(elec_struct.wavefunc(mu_c -elec_struct.mu_range[0])(i,ic,ik_c -elec_struct.ik_range[0]))* \
                                     elec_struct.wavefunc(mu_v -elec_struct.mu_range[0])(j,iv,ik_v -elec_struct.ik_range[0]) * \
                                     elec_struct.wavefunc(mu_cp-elec_struct.mu_range[0])(i,ic,ik_cp-elec_struct.ik_range[0]) * \
                           std::conj(elec_struct.wavefunc(mu_vp-elec_struct.mu_range[0])(j,iv,ik_vp-elec_struct.ik_range[0]))* \
                                                          _vq.data(ik_c_diff-_vq.iq_range[0],mu_c_diff-_vq.mu_range[0],2*i+j)/ \
                                                             _eps.data(ik_c_diff-_eps.iq_range[0],mu_c_diff-_eps.mu_range[0]);
      }
    }
    return dir_interaction;
  };

  // lambda function to calculate exchange interaction
  auto get_exchange_interaction = [&](){
    xch_interaction = 0;
    for (int i=0; i<2; i++)
    {
      for (int j=0; j<2; j++)
      {
        xch_interaction += std::conj(elec_struct.wavefunc(mu_c -elec_struct.mu_range[0])(i,ic,ik_c -elec_struct.ik_range[0]))* \
                                     elec_struct.wavefunc(mu_v -elec_struct.mu_range[0])(i,iv,ik_v -elec_struct.ik_range[0]) * \
                                     elec_struct.wavefunc(mu_cp-elec_struct.mu_range[0])(j,ic,ik_cp-elec_struct.ik_range[0]) * \
                           std::conj(elec_struct.wavefunc(mu_vp-elec_struct.mu_range[0])(j,iv,ik_vp-elec_struct.ik_range[0]))* \
                                                                  _vq.data(ik_cm-_vq.iq_range[0],mu_cm-_vq.mu_range[0],2*i+j);
      }
    }
    return xch_interaction;
  };

  // get ik of valence band state by taking care of wrapping around K2-extended zone
  auto get_ikv = [&elec_struct](const int& ik_c, const int& ik_cm){
    int ik_v = ik_c - ik_cm;
    while (ik_v >= elec_struct.ik_range[1]){
      ik_v -= elec_struct.nk;
    }
    while (ik_v < elec_struct.ik_range[0]){
      ik_v += elec_struct.nk;
    }
    return ik_v;
  };
  
  // get ik of conduction band state by taking care of wrapping around K2-extended zone
  auto get_ikc = [&elec_struct](const int& ik_v, const int& ik_cm){
    int ik_c = ik_v + ik_cm;
    while (ik_c >= elec_struct.ik_range[1]){
      ik_c -= elec_struct.nk;
    }
    while (ik_c < elec_struct.ik_range[0]){
      ik_c += elec_struct.nk;
    }
    return ik_c;
  };


  int nk_cm = ik_cm_range[1] - ik_cm_range[0];
  int nk_relev = int(_relev_ik_range[0].size());
  int nk_c = 2*nk_relev;

  arma::mat ex_energy_A1(nk_cm,nk_relev,arma::fill::zeros);
  arma::mat ex_energy_A2_singlet(nk_cm,nk_relev,arma::fill::zeros);
  arma::mat ex_energy_A2_triplet(nk_cm,nk_relev,arma::fill::zeros);
  
  arma::cx_cube ex_psi_A1(nk_c,nk_relev,nk_cm, arma::fill::zeros);
  arma::cx_cube ex_psi_A2_singlet(nk_c,nk_relev,nk_cm, arma::fill::zeros);
  arma::cx_cube ex_psi_A2_triplet(nk_c,nk_relev,nk_cm, arma::fill::zeros);

  arma::ucube ik_idx(4, nk_c, nk_cm);

  arma::cx_mat kernel_11(nk_relev,nk_relev,arma::fill::zeros);
  arma::cx_mat kernel_12(nk_relev,nk_relev,arma::fill::zeros);
  arma::cx_mat kernel_exchange(nk_relev,nk_relev,arma::fill::zeros);
  arma::vec energy;
  arma::cx_mat psi;
  arma::vec k_cm_vec(nk_cm,arma::fill::zeros);

  progress_bar prog(nk_cm, "calculate ex_energy");

  // loop to calculate exciton dispersion
  for (ik_cm=ik_cm_range[0]; ik_cm<ik_cm_range[1]; ik_cm++)
  {
    kernel_11.zeros();
    kernel_12.zeros();
    kernel_exchange.zeros();
    mu_cm = 0;
    int ik_cm_idx = ik_cm - ik_cm_range[0];
    k_cm_vec(ik_cm_idx) = ik_cm*arma::norm(_dk_l);

    prog.step(ik_cm_idx);


    for (int ik_c_idx=0; ik_c_idx<nk_relev; ik_c_idx++)
    {
      ik_c = _relev_ik_range[i_valley_1][ik_c_idx][0];
      mu_c = _relev_ik_range[i_valley_1][ik_c_idx][1];
      ik_v = get_ikv(ik_c,ik_cm);
      mu_v = mu_c;

      kernel_11(ik_c_idx, ik_c_idx) += elec_struct.energy(ic,ik_c-elec_struct.ik_range[0],mu_c-elec_struct.mu_range[0]) - \
                                       elec_struct.energy(iv,ik_v-elec_struct.ik_range[0],mu_c-elec_struct.mu_range[0]);

      // interaction  between valley_1 and valley_1
      for (int ik_cp_idx=0; ik_cp_idx<=ik_c_idx; ik_cp_idx++)
      {
        ik_cp = _relev_ik_range[i_valley_1][ik_cp_idx][0];
        mu_cp = _relev_ik_range[i_valley_1][ik_cp_idx][1];
        ik_vp = get_ikv(ik_cp,ik_cm);
        mu_vp = mu_cp;

        kernel_11(ik_c_idx,ik_cp_idx) -= get_direct_interaction();
        kernel_exchange(ik_c_idx,ik_cp_idx) += std::complex<double>(2,0)*get_exchange_interaction();
      }

      // interaction  between valley_1 and valley_2
      for (int ik_vp_idx=ik_c_idx; ik_vp_idx<nk_relev; ik_vp_idx++)
      {
        ik_vp = _relev_ik_range[i_valley_2][ik_vp_idx][0];
        mu_vp = _relev_ik_range[i_valley_2][ik_vp_idx][1];
        ik_cp = get_ikc(ik_vp,ik_cm);
        mu_cp = mu_vp;

        kernel_12(ik_c_idx,nk_relev-1-ik_vp_idx) -= get_direct_interaction();
      }
    }

    kernel_11 += kernel_11.t();
    kernel_12 += kernel_12.t();
    kernel_exchange += kernel_exchange.t();
    for (int ik_c_idx=0; ik_c_idx<nk_relev; ik_c_idx++)
    {
      kernel_11(ik_c_idx,ik_c_idx) /= std::complex<double>(2,0);
      kernel_12(ik_c_idx,ik_c_idx) /= std::complex<double>(2,0);
      kernel_exchange(ik_c_idx,ik_c_idx) /= std::complex<double>(2,0);
    }

    arma::eig_sym(energy,psi,kernel_11-kernel_12);
    ex_energy_A1.row(ik_cm_idx) = energy.t();
    ex_psi_A1.slice(ik_cm_idx).head_rows(nk_relev) = (+1/std::sqrt(2.))*psi;
    ex_psi_A1.slice(ik_cm_idx).tail_rows(nk_relev) = (-1/std::sqrt(2.))*psi;

    // energy = arma::eig_sym(kernel_11+kernel_12);
    arma::eig_sym(energy,psi,kernel_11+kernel_12);
    ex_energy_A2_triplet.row(ik_cm_idx) = energy.t();
    ex_psi_A2_triplet.slice(ik_cm_idx).head_rows(nk_relev) = (+1/std::sqrt(2.))*psi;
    ex_psi_A2_triplet.slice(ik_cm_idx).tail_rows(nk_relev) = (+1/std::sqrt(2.))*psi;

    // energy = arma::eig_sym(kernel_11+kernel_12+std::complex<double>(2,0)*kernel_exchange);
    arma::eig_sym(energy,psi,kernel_11+kernel_12+std::complex<double>(2,0)*kernel_exchange);
    ex_energy_A2_singlet.row(ik_cm_idx) = energy.t();
    ex_psi_A2_singlet.slice(ik_cm_idx).head_rows(nk_relev) = (+1/std::sqrt(2.))*psi;
    ex_psi_A2_singlet.slice(ik_cm_idx).tail_rows(nk_relev) = (+1/std::sqrt(2.))*psi;


    // save the index of kc and kv states from i_valley_1
    for (int ik_c_idx=0; ik_c_idx<nk_relev; ik_c_idx++)
    {
      ik_c = _relev_ik_range[i_valley_1][ik_c_idx][0];
      mu_c = _relev_ik_range[i_valley_1][ik_c_idx][1];
      ik_v = get_ikv(ik_c,ik_cm);
      mu_v = mu_c;

      ik_idx(0,ik_c_idx,ik_cm_idx) = ik_c-elec_struct.ik_range[0];
      ik_idx(1,ik_c_idx,ik_cm_idx) = mu_c-elec_struct.mu_range[0];
      ik_idx(2,ik_c_idx,ik_cm_idx) = ik_v-elec_struct.ik_range[0];
      ik_idx(3,ik_c_idx,ik_cm_idx) = mu_v-elec_struct.mu_range[0];
    }

    // save the index of kc and kv states from i_valley_2
    for (int ik_v_idx=0; ik_v_idx<nk_relev; ik_v_idx++)
    {
      ik_v = _relev_ik_range[i_valley_2][ik_v_idx][0];
      mu_v = _relev_ik_range[i_valley_2][ik_v_idx][1];
      ik_c = get_ikc(ik_v,ik_cm);
      mu_c = mu_v;

      ik_idx(0,nk_c-1-ik_v_idx,ik_cm_idx) = ik_c-elec_struct.ik_range[0];
      ik_idx(1,nk_c-1-ik_v_idx,ik_cm_idx) = mu_c-elec_struct.mu_range[0];
      ik_idx(2,nk_c-1-ik_v_idx,ik_cm_idx) = ik_v-elec_struct.ik_range[0];
      ik_idx(3,nk_c-1-ik_v_idx,ik_cm_idx) = mu_v-elec_struct.mu_range[0];
    }

  }

  std::cout << "\n...calculated exciton dispersion\n";

  std::cout << "saved exciton dispersion: A2 singlet\n";
  std::string filename = _directory.path()/"ex_energy_A2_singlet.dat";
  ex_energy_A2_singlet.save(filename, arma::arma_ascii);

  std::cout << "saved exciton dispersion: A2 triplet\n";
  filename = _directory.path()/"ex_energy_A2_triplet.dat";
  ex_energy_A2_triplet.save(filename, arma::arma_ascii);

  std::cout << "saved exciton dispersion: A1\n";
  filename = _directory.path()/"ex_energy_A1.dat";
  ex_energy_A1.save(filename, arma::arma_ascii);

  std::cout << "saved k_vector for center of mass\n";
  filename = _directory.path()/"exciton_k_cm_vec.dat";
  k_cm_vec.save(filename, arma::arma_ascii);

  // prepare the values that are to be returned
  std::vector<exciton_struct> excitons(3,exciton_struct(this));

  excitons[0].name = "A1 exciton";
  excitons[0].energy = ex_energy_A1;
  excitons[0].spin = 0;
  excitons[0].mu_cm = 0;
  excitons[0].n_principal = nk_relev;
  excitons[0].nk_c = nk_c;
  excitons[0].nk_cm = nk_cm;
  excitons[0].psi = ex_psi_A1;
  excitons[0].ik_idx = ik_idx;
  excitons[0].ik_cm_range = ik_cm_range;

  excitons[1].name = "A2 triplet exciton";
  excitons[1].energy = ex_energy_A2_triplet;
  excitons[1].spin = 1;
  excitons[1].mu_cm = 0;
  excitons[1].n_principal = nk_relev;
  excitons[1].nk_c = nk_c;
  excitons[1].nk_cm = nk_cm;
  excitons[1].psi = ex_psi_A2_triplet;
  excitons[1].ik_idx = ik_idx;
  excitons[1].ik_cm_range = ik_cm_range;

  excitons[2].name = "A2 singlet exciton";
  excitons[2].energy = ex_energy_A2_singlet;
  excitons[2].spin = 0;
  excitons[2].mu_cm = 0;
  excitons[2].n_principal = nk_relev;
  excitons[2].nk_c = nk_c;
  excitons[2].nk_cm = nk_cm;
  excitons[2].psi = ex_psi_A2_singlet;
  excitons[2].ik_idx = ik_idx;
  excitons[2].ik_cm_range = ik_cm_range;

  return excitons;
}

// call this to do all the calculations at once
void cnt::calculate_exciton_dispersion()
{
  get_parameters();
  get_atom_coordinates();

  // calculate K2-extended representation of electron energy
  std::array<int,2> ik_range_K2 = {0,_Nu/_Q*_nk_K1};
  std::array<int,2> mu_range_K2 = {0,_Q};
  _elec_K2 = electron_energy(ik_range_K2, mu_range_K2, "K2_extended");

  // find valleys and select a range of relevant iks in the focus valleys
  find_valleys(_elec_K2);
  find_relev_ik_range(1.*constants::eV, _elec_K2);

  // calculate vq, and dielectric function for a sufficiently large range of mu and ik.
  std::array<int,2> iq_range = {-(_elec_K2.ik_range[1]-1),_elec_K2.ik_range[1]};
  std::array<int,2> mu_range = {-(_elec_K2.mu_range[1]-1),_elec_K2.mu_range[1]};
  _vq = calculate_vq(iq_range, mu_range, _number_of_cnt_unit_cells);
  _PI = calculate_polarization(iq_range, mu_range, _elec_K2);
  _eps = calculate_dielectric(iq_range, mu_range);

  // calculate exciton dispersions using the information calculated above
  std::array<int,2> ik_cm_range = {-int(_relev_ik_range[0].size()), int(_relev_ik_range[0].size())};
  _excitons = calculate_A_excitons(ik_cm_range, _elec_K2);

}