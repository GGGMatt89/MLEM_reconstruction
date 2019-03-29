#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "image.h"
#include "aqphy.h"
#include "bgeot_ftool.h"

typedef enum { SM_PRECALC, SM_CALC, EV_PRECALC_DISK, EV_CALC_DISK, EV_PRECALC_RAM, EV_CALC_RAM } SM_management; 
typedef enum { ANY, RANGE, KNOWN } Energy_management;
typedef enum { OFF, CONSTANT, LINEAR, FORCE, FORCEINV, ONLYHODO, ALTERNATE } Hodoscope_management;
typedef enum { CV, RTS, RTV } Algorithm;
typedef enum { cos0rho0, cos1rho1, cos1rho2, cos0rho2 } Model;

class SystemMatrix
{
 private :
  int presel;
  int cpp;
  int cps;
  Algorithm algo;
  Model model;
  double width_factor;
  bool spatial_uncertainty;
  // energy
  int flag_energy;
  double  Emin;
  double Emax;
  double Etot;
  //double nb_vox[3];
  //int voxel_length[3];
  //double corner[3];
  //Image line;

 public:

 SystemMatrix() : cpp(0) {};
 SystemMatrix(int pres, int counts_p, int counts_s, Algorithm algorithm, Model flag_mod, double width, bool su) : presel(pres), cpp(counts_p), cps(counts_s), algo(algorithm), model(flag_mod),width_factor(width), spatial_uncertainty(su),flag_energy(ANY) {}

  // energy management strategy
  void set_energy_known(double etot)
  { 
    flag_energy=KNOWN;
    Etot=etot;
  }
  void set_energy_range(double emin, double emax)
  {
    flag_energy=RANGE;
    Emin=emin;
    Emax=emax;
  }

  // calculation
  bool algo_CV_uncertainties(Event &ev, Image &line) ;  // without spatial uncertainty
  bool algo_CV(Event &ev, Image &line) ; // with spatial uncertainty
  bool algo_RTS(Event &ev, Image &line) ;
  bool algo_RTV(Event &ev, Image &line) ;
  bool line_calc(Event &ev, Image &line) ;
  int sample_calc(std::ifstream& data_file, SM_management sm_management, std::vector<double> &ev_RAM, char *SM_filename, char *results_filename, const int nb_vox[], const double voxel_length[], const double corner[], Camera &camera, data_format f);
};


class MLEM
{
 private :
  int nb_iter;
  int first_iter;
  // image dimensions
  int nb_vox[3];
  double voxel_length[3];
  double corner[3];
  // energy
  int flag_energy;
  double  Emin;
  double Emax;
  double Etot;
  // hodoscope
  int flag_hodo;
  double beam_point[3];
  double beam_direction[3];
  int beam_sigma;
  int beam_width_factor;
  Hodoscope_management beam_mng;
  int beam_first_iter;
  int beam_nb_iter;
  // system matrix
  SM_management mng;
  Algorithm algo;
  Model model;
  double width_factor;
  bool spatial_uncertainty;
  // sensitivity 
  std::string file_sens;

 public:
 MLEM(): nb_iter(0) {}
 MLEM(int iter, int first, int nbvox[], double voxellength[], double corn[], double width) : nb_iter(iter), first_iter(first), flag_energy(ANY), flag_hodo(OFF), mng(SM_CALC),  width_factor(width)
  {
    nb_vox[0]=nbvox[0];
    nb_vox[1]=nbvox[1];
    nb_vox[2]=nbvox[2];
    voxel_length[0]=voxellength[0];
    voxel_length[1]=voxellength[1];
    voxel_length[2]=voxellength[2];
    corner[0]=corn[0];
    corner[1]=corn[1];
    corner[2]=corn[2];
  }

  void set_energy_range(double  emin, double emax)
  { 
    flag_energy=RANGE;
    Emin=emin;
    Emax=emax;
  }
  void set_energy_known(double etot)
  { 
    flag_energy=KNOWN;
    Etot=etot;
  }
  void set_hodoscope(int flag_hodoscope, double beam_pt[], double beam_dir[], int beam_sig, int beam_width, Hodoscope_management beam_incl_type, int bfiter, int bnbiter);
  void init_SM(SM_management m, Algorithm alg, Model mod, double width, bool su)
    { 
      mng=m;
      algo=alg;
      model=mod;
      width_factor=width;
      spatial_uncertainty=su;
    }
  const void hodoscope_matrix(Image &BM);

  //bool run_STORE(const char *results_file,  const char *SM_file,const std::vector<double> &sens, Image &BM, Camera &camera);
  bool run(const char *results_file,  const char *SM_file, const std::vector<double> &ev_RAM, const std::vector<double> &sens, Image &BM, Camera &camera, data_format f);

};
