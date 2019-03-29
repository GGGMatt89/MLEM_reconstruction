#include <iostream>
#include <fstream>
#include <string>


//using namespace std;

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <unistd.h>

#include "bgeot_ftool.h"
#include "gmm_std.h"
//#include "File_examinator.h"
//#include "traitement.h"
#include "image.h"
#include "aqphy.h"
#include "mlem.h"

//#include "gmm/gmm.h"

//#define seuil_pds 0.0001
//#define pi 3.14159265


//int nb_tra(std::string);

int main(int argc, char *argv[])
{

  if (argc<2)
    {
      std::cout<<"missing argument, needs parameters .m file.\n";
      return 1;
    }

  char *extension = strrchr( argv[1], '.');
  char * pwd = strrchr( argv[1], '/');
  std::string RunName( argv[1]);
  
  if (extension ==NULL || strncmp(extension, ".m",2)!=0)
    {
      std::cout<<"extension error, expected .m file. \n";
      return 2;
    }
  if (pwd!=NULL)
    {
      RunName.erase(0,pwd-argv[1]+1);
      RunName=RunName.substr(0,extension-pwd-1);
    }
  else
    {
      RunName=RunName.substr(0,extension-argv[1]);
    }
  std::cout <<"\n\nRun name: " <<RunName<< '\n';

  bgeot::md_param PARAM;
  PARAM.read_command_line(argc, argv);
  std::string data_file = PARAM.string_value("DATA_FILE","Data file ");
  std::cout << "Data from " << data_file << '\n';
  std::string data_type = PARAM.string_value("DATA_TYPE","Data type ");
  data_format df;
  if (data_type.compare("MEGALIB") ==0)
    df=MEGALIB;
  if (data_type.compare("IPNL") ==0)
    df=IPNL;
  if (data_type.compare("WP3") ==0)
    df=WP3;

  std::string camera_file(data_file);
  if (data_type.compare("WP3") ==0)
    {
      if (data_file.find( "_data.txt") ==std::string::npos)
	{
	  std::cout<<"extension error, for WP3 data expected _data.txt file. \n";
	  return 2;
	}
      camera_file.replace(data_file.find( "_data.txt"),5,"_camera_parameters");
      std::cout <<"Parameters of the camera in: " <<camera_file<< '\n';
    }

  std::string sens_file = PARAM.string_value("SENSITIVITY_FILE","file with the sensitivity matrix ");
  if (sens_file.empty()) 
    std::cout<< "No sensitivity file given\n";
  else 
    std::cout<< "Sensitivity matrix from "<< sens_file<<"\n\n";
  
  std::string results_path = PARAM.string_value("RESULTS_DIR","directory with the results ");
  std::cout<< "The results will be put in "<<results_path <<'\n';
  std::string matrix_path = PARAM.string_value("SM_DIR","directory with the system matrices ");
  struct stat buf;
  int ret;

  ret = stat (matrix_path.c_str(), &buf);
  if (ret == -1 && errno == ENOENT) 
    {
      ret = mkdir(matrix_path.c_str(), 0755);
      if (ret == -1) 
	{
	  std::cout<< "Could't create "<<  matrix_path << '\n';
	}
    } 
  ret = stat (results_path.c_str(), &buf);
  if (ret == -1 && errno == ENOENT) 
    {
      ret = mkdir(results_path.c_str(), 0755);
      if (ret == -1) 
	{
	  std::cout<< "Could't create "<<  results_path << '\n';
	}
    } 

  int nb_samples=PARAM.int_value("SAMPLES", "number of samples from data file ");
  std::string results_file(results_path+RunName);
  std::cout<< "Results files: "<< results_file <<".id*.r*.bin\n\n";
  std::string matrix_file(matrix_path+RunName);
  std::cout<< "The system matrices will be put in "<<matrix_path <<'\n';
  std::cout<< "System matrix files: "<< matrix_file << ".id*.SM*.bin\n";



  std::cout<< "\n~~~~~~ * COUNTS *~~~~~~~\n";
  int cps = PARAM.int_value("COUNTS_PER_SAMPLE", "number of events per sample ");
  int presel=PARAM.int_value("PRESELECT", "event selection ");
  if (presel>0)
    std::cout<< "Number of USEFUL events per sample (cones intersecting the volume): " << cps<<'\n';
  else 
    std::cout<< "Total number of DETECTED events per sample (not only cones intersecting the volume): " << cps<<'\n';
  int first = PARAM.int_value("FIRST", "event to begin with ");
  std::cout<< "Event to begin with: "<< first << "\n";
  std::cout<< "Number of samples from the data: " << nb_samples<< '\n';



  std::cout<< "\n~~~~~~ * VOLUME *~~~~~~~\n";
  double volume_length[3];
  std::vector<bgeot::md_param::param_value> tmp= PARAM.array_value("VOLUME_DIMENSIONS", "dimensions of the volume");
  std::cout << "Volume dimensions: [ " ;
  for (int i=0; i<3; i++) 
    {
      volume_length[i]=tmp[i].real();
      std::cout<<volume_length[i] << " ";
    }
  std::cout << "]\n";
  int nb_vox[3];
  tmp= PARAM.array_value("VOXELS", "nb voxels of the volume");
  std::cout << "Voxels in the volume: [ " ;
  for (int i=0; i<3; i++) 
    {
      nb_vox[i]=int(tmp[i].real());
      std::cout<<nb_vox[i] << " ";
    }
  std::cout << "]\n";
  double voxel_length[3] ;
  voxel_length[0]=volume_length[0]/nb_vox[0];
  voxel_length[1]=volume_length[1]/nb_vox[1];
  voxel_length[2]=volume_length[2]/nb_vox[2];
  int total_nb_voxels=nb_vox[0]*nb_vox[1]*nb_vox[2]; 
  double frame_center[3];
  tmp= PARAM.array_value("VOLUME_CENTRE", "centre of the volume");
  std::cout << "Volume centre: [ " ;
  for (int i=0; i<3; i++) 
    {
      frame_center[i]=tmp[i].real();
      std::cout<<frame_center[i] << " ";
    }
  std::cout << "]\n";
  double corner[3];
  corner[0]=frame_center[0]-volume_length[0]/2;
  corner[1]=frame_center[1]-volume_length[1]/2;
  corner[2]=frame_center[2]-volume_length[2]/2;



  std::cout<< "\n~~~~~~ * ITERATIONS *~~~~~~~\n";
  SM_management sm_management;
  int cpp=0;
  std::string sa=PARAM.string_value("TYPE_SM", "System matrix management ");;
  bool flag_test=(sa=="SM_PRECALC")||(sa=="SM_CALC")||(sa=="EV_PRECALC_DISK")||(sa== "EV_CALC_DISK")||(sa== "EV_PRECALC_RAM")||(sa=="EV_CALC_RAM"); 
  if (!flag_test)  
    { 
      std::cout<< "System matrix management: error\n";
      return 1;
    }
  if (sa=="SM_CALC")
    {
      std::cout<< "The system matrix will be calculated and stored\n";
      sm_management=SM_CALC;
    };
  if (sa=="SM_PRECALC")
    {
      std::cout<< "The system matrices are supposed to be precalculated \n";
      sm_management=SM_PRECALC;
    };
  if (sa=="EV_PRECALC_DISK")
    {
      std::cout<< "The preprocessed events file is supposed to be stored; it will be read line by line \n";
      sm_management=EV_PRECALC_DISK;
    };
  if (sa=="EV_PRECALC_RAM")
    {
      std::cout<< "The preprocessed events file is supposed to be stored; it will be load in the RAM \n";
      sm_management=EV_PRECALC_RAM;
    };
  if (sa=="EV_CALC_DISK")
    {
      std::cout<< "The preprocessed events file will be calculated and stored; it will be read line by line \n";
      sm_management=EV_CALC_DISK;
    };
  if (sa=="EV_CALC_RAM")
    {
      std::cout<< "The preprocessed events will be calculated and kept in the RAM \n";
      sm_management=EV_CALC_RAM;
    };
  bool clear_SM=PARAM.int_value("CLEAR_SM", "System matrix or preprocessed events file will be erased or not 1|0") ;
  if ((sm_management==SM_PRECALC)||(sm_management==SM_CALC))
    cpp=PARAM.int_value("COUNTS_PER_PARTITION","Counts per partition");
  else 
    cpp=cps;
  
  int iter_first = PARAM.int_value("FIRST_ITERATION", "first iteration");
  if (iter_first<1) 
    {
      std::cout << "FIRST_ITERATION should not be < 1; set to 1\n";
      iter_first=1;
    }
  int iter_nb = PARAM.int_value("ITERATIONS", "number of iterations");
  std::cout << "Number of iterations: " << iter_nb << "\n";


  //int algorithm = PARAM.string_value("ALGORITHM", "algorithm");
  sa=PARAM.string_value("ALGORITHM", "algorithm");
  Algorithm algorithm;
  flag_test=(sa=="CV")||(sa=="RTS")||(sa== "RTV"); 
  if (!flag_test)  
    { 
      std::cout<< "Algorithm: error\n";
      return 1;
    };
  if (sa=="CV")
    {
      std::cout<< "Algorithm: CV \n";
      algorithm=CV;
    };
  if (sa=="RTS")
    {
      std::cout<< "Algorithm: RayTracing Surface \n";
      algorithm=RTS;
    };
  if (sa== "RTV")
    {
      std::cout<< "Algorithm: RayTracing Volume \n";
      algorithm=RTV;
    };
  if (sa== "")    
    { std::cout<< "Unknown algorithm\n";
      return 1;
    };


  sa =  PARAM.string_value("MODEL", "model"); 
  Model model;
  flag_test=(sa=="cos0rho0")||(sa=="cos1rho1")||(sa== "cos1rho2")||(sa== "cos0rho2"); 
  if (!flag_test)  
    { 
      std::cout<< "Unknown model\n";
      return 1;
    };
  if (sa=="cos0rho0")
    {
      std::cout<< "Model: cos0rho0 \n";
      model=cos0rho0;
    };
  if (sa=="cos1rho1")
    {
      std::cout<< "Model: cos1rho1 \n";
      model=cos1rho1;
    };
  if (sa== "cos1rho2")
    {
      std::cout<< "Model: cos1rho2 \n";
      model=cos1rho2;
    };
  if (sa== "cos0rho2")
    {
      std::cout<< "Model: cos0rho2 \n";
      model=cos0rho2;
    };

  double width_factor = PARAM.real_value("WIDTH_FACTOR", " thicker cones ");



  std::cout<< "\n~~~~~~ * ENERGY SELECTION *~~~~~~~\n";
  int energy_selection = PARAM.int_value("ENERGY_FLAG", "energy flag ");
  Energy_management flag_energy;
  double  Emin=-1, Emax=-1, Etot=-1;
  switch (energy_selection)
    {
    case 0:
      std::cout<< "Energy selection OFF \n" ;
      flag_energy=ANY;
      break;
    case 1:
      std::cout<< "Energy selection ON: acceptance interval \n" ;
      flag_energy=RANGE;
      Emin = PARAM.real_value("ENERGY_MIN", "minimum allowed total energy"); 
      Emax = PARAM.real_value("ENERGY_MAX", "maximum allowed total energy"); 
      std::cout<< "Minimal allowed total energy: " << Emin <<" keV\n";
      std::cout<< "Maximal allowed total energy: " << Emax <<" keV\n";
      break;
    case 2:
      std::cout<< "Energy selection ON: known total energy \n" ;
      flag_energy=KNOWN;
      Etot = PARAM.real_value("ENERGY_TOTAL", "known initial energy"); 
      std::cout<< "Initial energy of the photons: " << Etot <<" keV \n";
      break;
    default :
      std::cout<< "Unknown energy selection parameter \n" ;
    }



  std::cout<< "\n~~~~~~ * HODOSCOPE *~~~~~~~\n";
  double beam_point[3];
  double beam_direction[3];
  double beam_sigma=0, beam_width_factor=0;
  int beam_first_iter, beam_nb_iter;
  std::string beam_inclusion;
  Hodoscope_management beam_inclusion_type=OFF;

  int flag_hodoscope = PARAM.int_value("HODOSCOPE_FLAG", "hodoscope ");
  if (flag_hodoscope)
    {
      std::cout<< "Hodoscope ON \n" ;
      tmp= PARAM.array_value("BEAM_ENTRY_POINT", "point on the beam line");
      std::cout << "Beam entry point: [ " ;
      for (int i=0; i<3; i++) 
	{
	  beam_point[i]=tmp[i].real();
	  std::cout<<beam_point[i] << " ";
	}
      std::cout << "]\n";
      tmp= PARAM.array_value("BEAM_DIRECTION", "direction of the beam line");
      std::cout << "Beam direction: [ " ;
      for (int i=0; i<3; i++) 
	{
	  beam_direction[i]=tmp[i].real();
	  std::cout<<beam_direction[i] << " ";
	}
      std::cout << "]\n"; 
      beam_sigma = PARAM.real_value("BEAM_SIGMA", "sigma for beam "); 
      beam_width_factor = PARAM.real_value("BEAM_WIDTH_FACTOR", "nb sigma for beam ");
      std::cout << "sigma for beam = "<< beam_sigma<< "\n"; 
      std::cout << "Number of sigma for beam = "<< beam_width_factor<< "\n\n\n"; 
      beam_first_iter = PARAM.int_value("BEAM_FIRST_ITERATION", "Inclusion of the beam a priori from iter "); 
      std::cout << "Inclusion of the beam a priori from iter "<< beam_first_iter<< "\n"; 
      beam_nb_iter = PARAM.int_value("BEAM_ITERATIONS", "Inclusion of the beam a priori for ? iterations "); 
      std::cout << "Inclusion of the beam a priori on "<< beam_nb_iter<< " iterations\n"; 
      beam_inclusion = PARAM.string_value("BEAM_INCLUSION", "Inclusion of the beam a priori type "); 
      if (beam_inclusion=="CONSTANT")
	{
	  std::cout << "Inclusion type CONSTANT\n"; 
	  beam_inclusion_type=CONSTANT;
	}
      if (beam_inclusion=="LINEAR")
	{
	  std::cout << "Inclusion type LINEAR\n"; 
	  beam_inclusion_type=LINEAR;
	}
      if (beam_inclusion=="FORCE")
	{
	  std::cout << "Inclusion type FORCE line\n"; 
	  beam_inclusion_type=FORCE;
	}
      if (beam_inclusion=="FORCEINV")
	{
	  std::cout << "Inclusion type FORCEINV line\n"; 
	  beam_inclusion_type=FORCEINV;
	}
      if (beam_inclusion=="ONLYHODO")
	{
	  std::cout << "Inclusion type ONLYHODO line\n"; 
	  beam_inclusion_type=ONLYHODO;
	}
      if (beam_inclusion=="ALTERNATE")
	{
	  std::cout << "Inclusion type ALTERNATE\n"; 
	  beam_inclusion_type=ALTERNATE;
	}
    }
  else 
    std::cout<< "Hodoscope OFF \n" ;


 

  // Config camera --------------------------------------------

  std::cout<<"\n\n---*Configuration of the Compton Camera*---\n";
  Camera camera(df);
  int nb_layers=PARAM.int_value("NB_LAYERS", "Number of layesr in scatterer ");
  char name[30], label[100] ;
  Point centre;
  std::string empty("");
  for (int i=1; i<=nb_layers; i++)
    {
      sprintf(name, "LAY_CENTRE_%d",i);
      sprintf(label, "centre of layer %d ",i);
      tmp= PARAM.array_value(name, label);
      centre[0]=tmp[0].real();
      centre[1]=tmp[1].real();
      centre[2]=tmp[2].real();
      camera.add_sca_layer(centre,empty);
    }
  tmp= PARAM.array_value("ABS_CENTRE", " centre of the absorber ");
  centre[0]=tmp[0].real();
  centre[1]=tmp[1].real();
  centre[2]=tmp[2].real();
  camera.add_abs(centre,empty);
  tmp= PARAM.array_value("LAY_SIZE", " dimensions of a layer ");
  camera.set_layer_dim(tmp[0].real(),tmp[1].real(),tmp[2].real());
  tmp= PARAM.array_value("ABS_SIZE", " dimensions of the absorber ");
  camera.set_abs_dim(tmp[0].real(),tmp[1].real(),tmp[2].real());
  tmp= PARAM.array_value("LAY_VOXELS", " nb detector units ");
  camera.set_layer_vox(int(tmp[0].real()),int(tmp[1].real()),int(tmp[2].real()));
  tmp= PARAM.array_value("ABS_VOXELS", " nb detector units ");
  camera.set_abs_vox(int(tmp[0].real()),int(tmp[1].real()),int(tmp[2].real()));
  bool spatial_uncertainty= (PARAM.int_value("SPATIAL_UNCERTAINTY", " flag for the spatial uncertainty consideration ")==1);
  if (spatial_uncertainty)
    {
      tmp= PARAM.array_value("LAY_VOXEL_SAMPLING", " nb points in a voxel ");
      camera.set_layer_vox_sampling(int(tmp[0].real()),int(tmp[1].real()),int(tmp[2].real()));
      tmp= PARAM.array_value("ABS_VOXEL_SAMPLING", " nb points in a voxel ");
      camera.set_abs_vox_sampling(int(tmp[0].real()),int(tmp[1].real()),int(tmp[2].real()));
    }
  Point Ox, Oy, Oz;
  tmp= PARAM.array_value("Ox", "Local frame (camera related), vector parallel to an edge ");
  for (int i=0; i<3; i++) 
    Ox[i]=tmp[i].real();
  tmp= PARAM.array_value("Oy", "Local frame (camera related), vector parallel to an edge ");
  for (int i=0; i<3; i++) 
    Oy[i]=tmp[i].real();
  tmp= PARAM.array_value("Oz", "Local frame (camera related), vector parallel to an edge ");
  for (int i=0; i<3; i++) 
    Oz[i]=tmp[i].real();
  if (!camera.set_frame(Ox,Oy,Oz))
    {
      std::cout<< "Error in camera.set_frame: axis vectors cannot be null\n";
      return 1;
    }
  if(!camera.set_sca_box()) 
    {
      std::cout<< "Error in camera.set_box: the frame is not yet defined\n";
      return 1;
    }

  int ok;
  /*
  if (camera.get_format()==WP3)
    ok=camera.read_attributes(camera_file);
  else 
    ok=camera.read_attributes(data_file);
  if (ok>0)
    {
      std::cout << "Error in camera file \n";
      return 1;
    }
  */
  
  std::cout<<"Number of layers: "<<camera.get_layer_tab_size()<<'\n';
  std::vector<double> dim=camera.get_layer_dim();
  std::cout<<"Size of scatterer layers: "<<dim[0]<<"x"<<dim[1]<<"x"<<dim[2]<<'\n';
  dim=camera.get_sca_dim();
  std::cout<<"Size of the scatterer: "<<dim[0]<<"x"<<dim[1]<<"x"<<dim[2]<<'\n';
  dim=camera.get_sca_centre();
  std::cout<<"Centre of the scatterer: ("<<dim[0]<<","<<dim[1]<<","<<dim[2]<<")\n";
  dim=camera.get_abs_dim();
  std::cout<<"Size of the absorber: "<<dim[0]<<"x"<<dim[1]<<"x"<<dim[2]<<'\n';
  dim=camera.get_abs();
  std::cout<<"Centre of the absorber: ("<<dim[0]<<","<<dim[1]<<","<<dim[2]<<")\n";

  for (int i=0;i<camera.get_layer_tab_size(); i++)
    {
      std::cout<<"Scatterer " << i<< " " << camera.get_layer_name(i) << " at : ["<<camera.get_sca_layer(i)[0]<<","<<camera.get_sca_layer(i)[1]<<","<<camera.get_sca_layer(i)[2]<<"]\n";
    }


  // Sensitivity  ----------------------------------------------

  std::vector<double> sens(nb_vox[0]*nb_vox[1]*nb_vox[2],1);
  if(!sens_file.empty())
    { 
      std::ifstream FS;
      FS.open(sens_file.c_str(), std::ios::in|std::ios::binary);
      if (!FS.is_open())
	{
	  std::cout<<"Error while opening the sensitivity matrix file\n";
	  return -1;
	}
      FS.read((char *)(&(sens[0])), nb_vox[0]*nb_vox[1]*nb_vox[2]*sizeof(double));
      if (!FS) 
	{
	  std::cout << "Error in " << sens_file << ": only "<< FS.gcount()<< " items could be read.\n" ;
	  return false;
	}
      FS.close();
    }


  // Prepare MLEM ------------------------------------------------------
  // open the data file, only whem the system matrix is not precalculated
  std::ifstream dfp;
  bool end_of_series=false;
  Event ev(df); 
  if ((sm_management==SM_CALC)||(sm_management==EV_CALC_DISK)||(sm_management==EV_CALC_RAM))
    {
      dfp.open(data_file.c_str());
      if (!dfp)
	{
	  std::cout<<"Error while opening the data file\n";
	  return 1;
	}  
      int total_read=1;
      if (df==WP3)  
	{
	  std::string line;  
	  getline(dfp,line);  // names of the columns
	}
      while (!dfp.eof() && (total_read<first))
	{
	  ok=ev.read_event(dfp);
	  total_read++;
	}
      if (dfp.eof()) 
	{
	  std::cout<< "The data file contains less than first=" << first << " elements. Nothing to calculate.\n";
	  return 1;
	}
    }

  // initialisation of the MLEM object
  MLEM mlem(iter_nb,iter_first,nb_vox,voxel_length,corner, width_factor);
  if ((sm_management==EV_PRECALC_DISK)||(sm_management==EV_PRECALC_RAM)||(sm_management==EV_CALC_DISK)||(sm_management==EV_CALC_RAM))
    mlem.init_SM(sm_management,algorithm,model, width_factor,spatial_uncertainty);
  if (flag_hodoscope)
    mlem.set_hodoscope(flag_hodoscope,beam_point,beam_direction,beam_sigma,beam_width_factor, beam_inclusion_type, beam_first_iter, beam_nb_iter);
  
  // calculate the hodoscope matrix
  Image beam_matrix;
  if (flag_hodoscope) 
    {
      beam_matrix=Image(nb_vox,voxel_length, corner);
      mlem.hodoscope_matrix(beam_matrix);
    }

  // Run over series --------------------------------------------------
  bool end_of_data_file=false, end_of_events=false;
  bool flag_mlem;
  int sample_id=0;
  char sample_results_file[100];
  std::vector<double> ev_RAM;
  std::cout << "\n\n Run over samples ......\n";
  while ((!end_of_series)&&(sample_id<nb_samples))
    {
      std::cout << "\n Sample "<< sample_id<<'\n';
      // Define filenames --------------------------------------------
      sprintf(sample_results_file, "%s.id%d",results_file.c_str(),sample_id);
      char work_file[100];
      if ((sm_management==SM_PRECALC) || (sm_management==SM_CALC)) 
	sprintf(work_file, "%s.id%d",matrix_file.c_str(),sample_id);
      if ((sm_management==EV_PRECALC_DISK) || (sm_management==EV_PRECALC_RAM) || (sm_management==EV_CALC_DISK))
	sprintf(work_file,"%s.id%d.ev.txt", matrix_file.c_str(), sample_id);
      // If necessary, calculate SM or the preprocessed events ------
      clock_t t=clock();
      if (sm_management==EV_PRECALC_RAM)
	{
	  end_of_events=!ev.load_post_treated_all(work_file,ev_RAM);
	  //std::cout << "vector size" <<ev_RAM.size()<<'\n';
	}
      if ((sm_management==SM_CALC)||(sm_management==EV_CALC_DISK)||(sm_management==EV_CALC_RAM))
	{
	  SystemMatrix SM(presel, cpp,cps,algorithm, model, width_factor,spatial_uncertainty);
       //std::cout <<"dfp "<<dfp<<'\n';
         std::cout <<"test "<<'\n';
	  if (flag_energy==RANGE) SM.set_energy_range(Emin,Emax);
	  if (flag_energy==KNOWN) SM.set_energy_known(Etot);
        
        //<<" sm_management "<<sm_management<<" ev_RAM "<<ev_RAM<<" work_file "<<work_file<<" sample_results_file "<<sample_results_file<<" nb_vox "<<nb_vox<<" voxel_length "<<voxel_length<<" corner "<<corner<<" camera "<<camera<<"\n";
	  end_of_data_file=(SM.sample_calc(dfp, sm_management, ev_RAM, work_file,sample_results_file, nb_vox,voxel_length,corner,camera,df)==1);
        
	}
      t=clock()-t;
      std::cout << "Time spent in calculation of the system matrix or events: "<< ((float)t)/CLOCKS_PER_SEC<< " seconds\n";
      // Run MLEM on this subset ------------------------------
      flag_mlem=mlem.run(sample_results_file, work_file,ev_RAM,sens,beam_matrix,camera,df);
      end_of_series=end_of_data_file || end_of_events || !flag_mlem;
      // Delete files containing the system matrix for the sample
      if (clear_SM && (sm_management!=EV_CALC_RAM))
	{
	  int part=0;
	  std::ifstream FSM;
	  char filename[100];
	  while (true)
	    {
	      if ((sm_management==SM_PRECALC) || (sm_management==SM_CALC)) 
		sprintf(filename,"%s.SM%d.bin",work_file, part);
	      if ((sm_management==EV_CALC_DISK)||(sm_management==EV_PRECALC_RAM) || (sm_management==EV_PRECALC_DISK))
		sprintf(filename,"%s",work_file);
	      FSM.open(filename,std::ios::in| std::ios::binary);
	      if (FSM)
		{
		  FSM.close();
		  unlink(filename);
		}
	      else break;
	      part++;
	    }
	}
      sample_id++;
    }

  if (dfp.is_open()) dfp.close();

  return 0;




}

