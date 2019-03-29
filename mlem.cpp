
#include <math.h>

#include "mlem.h"
#include "aqphy.h"



bool SystemMatrix::algo_CV_uncertainties(Event &ev, Image &line) 
{
  // updated version of algo_CV 
  // that accounts for spatial uncertainties
  bool non_zero=false;
  double taille_max=std::max(std::max(line.VoxelSize[0],line.VoxelSize[1]),line.VoxelSize[2]) ;   
  Point V1=ev.get_V1(), V2=ev.get_V2();
  std::vector<Point> S1,S2;
  Camera Cc=ev.get_camera();
  if (spatial_uncertainty)
    {
      if (!Cc.layer_sample(V1,S1)) return false;
      if (!Cc.abs_sample(V2,S2)) 
	{
	  if (!Cc.layer_sample(V2,S2))
	    return false;
	}
    }
  else 
    {
      S1.push_back(V1);
      S2.push_back(V2);
    }
  Point OV1;
  OV1[0] = line.Corner[0]+line.DimInCm[0]/2-V1[0];
  OV1[1] = line.Corner[1]+line.DimInCm[1]/2-V1[1];
  OV1[2] = line.Corner[2]+line.DimInCm[2]/2-V1[2];
  double seuil = 2*taille_max/OV1.norm2();
    //sqrt(OV1[0]*OV1[0]+OV1[1]*OV1[1]+OV1[2]*OV1[2]);       //valeur du seuil
  double sigma_beta=width_factor*ev.get_sigma_beta();
  if(sigma_beta<seuil) 
    {
      sigma_beta=seuil;
    }
  //int nb_sigma=3;
  double fact_gauss=1./(sigma_beta*sqrt(2*M_PI));
  double Oj[3];
  double rho_m,delta_m,theta_m, gauss, kbl_m, norm, alpha;
  std::vector<Point> cone_axis(S1.size()*S2.size());
  if (spatial_uncertainty)
    {
      for (int n=0; n<S2.size();n++)
	for (int m=0; m<S1.size(); m++)
	  {
	    cone_axis[m+n*S1.size()][0]=S1[m][0]-S2[n][0];
	    cone_axis[m+n*S1.size()][1]=S1[m][1]-S2[n][1];
	    cone_axis[m+n*S1.size()][2]=S1[m][2]-S2[n][2];
	    norm=(cone_axis[m+n*S1.size()]).norm2();
	    cone_axis[m+n*S1.size()][0] /=norm;
	    cone_axis[m+n*S1.size()][1] /=norm;
	    cone_axis[m+n*S1.size()][2] /=norm;
	  }
    }
  else 
    {
      ev.get_cone_axis(cone_axis[0]);
    }
  Point camera_axis=ev.get_normal_to_camera();
  std::vector<double> dv_1(S1.size());
  std::vector<double> dv_2(S1.size());
  std::vector<double> delta_1(S1.size()*S2.size());
  std::vector<double> delta_2(S1.size()*S2.size());
  std::vector<double> theta_1(S1.size());
  std::vector<double> theta_2(S1.size());

    
  /*std::cout << "beta="<< ev.get_beta()*180/M_PI << " , sigma_beta=" << sigma_beta*180/M_PI << '\n';
  std::cout << "cone axis "<< cone_axis[0]<< ' ' << cone_axis[1]<< ' '<<cone_axis[2]<< " alpha="<< ev.get_alpha()*180/M_PI<<'\n';
  std::cout << "camera axis "<< camera_axis[0]<< ' ' << camera_axis[1]<< ' '<<camera_axis[2]<< '\n';
  */
  for(int k=0; k<line.DimInVoxels[2]; k++)
    {
      Oj[2]=line.Corner[2]+(k+0.5)*line.VoxelSize[2];
      for (int m=0; m<S1.size(); m++)
	{
	  dv_2[m]=(Oj[2]-S1[m][2])*(Oj[2]-S1[m][2]);
	  theta_2[m]=camera_axis[2]*(Oj[2]-S1[m][2]);
	  for (int n=0; n<S2.size(); n++)
	    delta_2[m+n*S1.size()]=cone_axis[m+n*S1.size()][2]*(Oj[2]-S1[m][2]);
	}
      for(int j=0; j<line.DimInVoxels[1]; j++)
        {
	  Oj[1]=line.Corner[1]+(j+0.5)*line.VoxelSize[1];
	  for (int m=0; m<S1.size(); m++)
	    {
	      dv_1[m]=(Oj[1]-S1[m][1])*(Oj[1]-S1[m][1]);
	      theta_1[m]=camera_axis[1]*(Oj[1]-S1[m][1]);
	      for (int n=0; n<S2.size(); n++)
		delta_1[m+n*S1.size()]=cone_axis[m+n*S1.size()][1]*(Oj[1]-S1[m][1]);
	    }
	  for(int i=0; i<line.DimInVoxels[0]; i++)
            {
	      Oj[0]=line.Corner[0]+(i+0.5)*line.VoxelSize[0];
	      for (int m=0; m<S1.size(); m++)
		{
		  rho_m=sqrt((Oj[0]-S1[m][0])*(Oj[0]-S1[m][0])+dv_1[m]+dv_2[m]);                          // ρ=||V₁Oⱼ||
		  theta_m=acos((camera_axis[0]*(Oj[0]-S1[m][0])+theta_1[m]+theta_2[m])/rho_m);  // θ=arcos((V₁Oⱼ.z)/ρ)
		  for (int n=0; n<S2.size(); n++)
		    {
		      delta_m=acos((cone_axis[m+n*S1.size()][0]*(Oj[0]-S1[m][0])+delta_1[m+n*S1.size()]+delta_2[m+n*S1.size()])/rho_m);                     // δ=arcos((N.V₁Oⱼ)/ρ)	      
	      //if (fabs(delta_m-ev.get_beta())<=nb_sigma*sigma_beta)
	      //{                 
		      gauss=fact_gauss*exp(-(delta_m-ev.get_beta())*(delta_m-ev.get_beta())/(2*sigma_beta*sigma_beta));
		      kbl_m=ev.get_P()*ev.get_P()*0.5*(ev.get_P()+1/ev.get_P()-1+cos(delta_m)*cos(delta_m));
		      switch(model)
			{
			case cos1rho1:
			  alpha= acos((cone_axis[m+n*S1.size()]).dot_product(Cc.get_Oz()));
			  kbl_m *=fabs(cos(theta_m)*cos(alpha))/rho_m;
			  break;
			case cos1rho2:
			  alpha= acos((cone_axis[m+n*S1.size()]).dot_product(Cc.get_Oz()));
			  kbl_m *=(fabs(cos(theta_m)*cos(alpha)))/(rho_m*rho_m);                       // Model_rho2
			  break;
			case cos0rho2:
			  kbl_m *=1/(rho_m*rho_m);      // Model_rho2 sans les cosinus
			  break;
			}
		      line.Value[line.index_3Dto1D(i,j,k)] +=kbl_m*gauss/(S1.size()*S2.size());
		    }
		}
	      //non_zero=true;
		  //};
            }
        }
    }
  return true;
}


bool SystemMatrix::algo_CV(Event &ev, Image &line) 
{
  //do not accounts for spatial uncertainties;
  //see new version algo_CV
  bool non_zero=false;
  double taille_max=std::max(std::max(line.VoxelSize[0],line.VoxelSize[1]),line.VoxelSize[2]) ;   
  Point V1=ev.get_V1(), V2=ev.get_V2();
  Point OV1;
  OV1[0] = line.Corner[0]+line.DimInCm[0]/2-V1[0];
  OV1[1] = line.Corner[1]+line.DimInCm[1]/2-V1[1];
  OV1[2] = line.Corner[2]+line.DimInCm[2]/2-V1[2];
  double seuil = 2*taille_max/OV1.norm2();
    //sqrt(OV1[0]*OV1[0]+OV1[1]*OV1[1]+OV1[2]*OV1[2]);       //valeur du seuil
  double sigma_beta=width_factor*ev.get_sigma_beta();
  if(sigma_beta<seuil) 
    {
      sigma_beta=seuil;
    }
  int nb_sigma=3;
  double fact_gauss=1./(sigma_beta*sqrt(2*M_PI));
  double inv_var=0.5*1./(sigma_beta*sigma_beta);
  double Oj[3];
  double rho_m, cos_delta_m,delta_m, theta_m, gauss, kbl_m;
  double ddelta, cos_theta_m;
  Point cone_axis;
  ev.get_cone_axis(cone_axis);
  double P2=ev.get_P()*ev.get_P();
  double K_part=0.5*ev.get_P()*(1-ev.get_P()+P2);
  double cos_alpha=cos(ev.get_alpha());

  Point camera_axis=ev.get_normal_to_camera();
  double dv_0,dv_1,dv_2,delta_2,delta_1, theta_1,theta_2;
    
  /*std::cout << "beta="<< ev.get_beta()*180/M_PI << " , sigma_beta=" << sigma_beta*180/M_PI << '\n';
  std::cout << "cone axis "<< cone_axis[0]<< ' ' << cone_axis[1]<< ' '<<cone_axis[2]<< " alpha="<< ev.get_alpha()*180/M_PI<<'\n';
  std::cout << "camera axis "<< camera_axis[0]<< ' ' << camera_axis[1]<< ' '<<camera_axis[2]<< '\n';
  */

  Oj[2]=line.Corner[2]+0.5*line.VoxelSize[2];
  for(int k=0; k<line.DimInVoxels[2]; k++)
    {
      //Oj[2]=line.Corner[2]+(k+0.5)*line.VoxelSize[2];
      dv_2=Oj[2]-V1[2];
      delta_2=dv_2*cone_axis[2]; //cone_axis[2]*(Oj[2]-V1[2])
      theta_2=dv_2*camera_axis[2]; // camera_axis[2]*(Oj[2]-V1[2]);
      dv_2 *=dv_2;  // (Oj[2]-V1[2])*(Oj[2]-V1[2])
      Oj[1]=line.Corner[1]+0.5*line.VoxelSize[1];
      for(int j=0; j<line.DimInVoxels[1]; j++)
        {
	  //Oj[1]=line.Corner[1]+(j+0.5)*line.VoxelSize[1];
	  dv_1=Oj[1]-V1[1];
	  delta_1=dv_1*cone_axis[1];  // cone_axis[1]*(Oj[1]-V1[1]);
	  theta_1=dv_1*camera_axis[1]; // camera_axis[1]*(Oj[1]-V1[1]);
	  dv_1 *=dv_1; // (Oj[1]-V1[1])*(Oj[1]-V1[1]);
	  Oj[0]=line.Corner[0]+0.5*line.VoxelSize[0];
	  for(int i=0; i<line.DimInVoxels[0]; i++)
            {
	      //Oj[0]=line.Corner[0]+(i+0.5)*line.VoxelSize[0];
	      dv_0=Oj[0]-V1[0];
	      rho_m=sqrt(dv_0*dv_0+dv_1+dv_2);                                    // ρ=||V₁Oⱼ||
	      cos_delta_m=(cone_axis[0]*dv_0+delta_1+delta_2)/rho_m;
	      delta_m=acos(cos_delta_m);                     // δ=arcos((N.V₁Oⱼ)/ρ)	      
	      ddelta=fabs(delta_m-ev.get_beta());
	      if (ddelta<=nb_sigma*sigma_beta)
		{
		  cos_theta_m=(camera_axis[0]*dv_0+theta_1+theta_2)/rho_m;
		  theta_m=acos(cos_theta_m);                    // θ=arcos((V₁Oⱼ.z)/ρ)
		  gauss=fact_gauss*exp(-ddelta*ddelta*inv_var);
		  kbl_m=K_part+0.5*P2*cos_delta_m*cos_delta_m;
		  switch(model)
		    {
		    case cos1rho1:
		      kbl_m *=fabs(cos_theta_m*cos_alpha)/rho_m;
		      break;
		    case cos1rho2:
		      kbl_m *=fabs(cos_theta_m*cos_alpha)/(rho_m*rho_m);                       // Model_rho2
		      break;
		    case cos0rho2:
		      kbl_m *=1/(rho_m*rho_m);      // Model_rho2 sans les cosinus
		      break;
		    }
		  line.Value[line.index_3Dto1D(i,j,k)]=kbl_m*gauss;
		  non_zero=true;
		}
	      Oj[0] += line.VoxelSize[0];
            }
	  Oj[1] += line.VoxelSize[1];
        }
      Oj[2] += line.VoxelSize[2];
    }
  return non_zero;
}


bool SystemMatrix::algo_RTS(Event &ev, Image &line) 
{
  
}

bool SystemMatrix::algo_RTV(Event &ev, Image &line) 
{

}

bool SystemMatrix::line_calc(Event &ev, Image &line)
{
  // returns true if the line is not full of zeroes
  line.initialize(0);
  switch (algo)
    {
    case CV :
      if (spatial_uncertainty)
	return algo_CV_uncertainties( ev, line);
      else
	return algo_CV( ev, line);
      break;
    case RTS:
      return algo_RTS( ev, line);
      break;
    case RTV:
      return algo_RTV( ev, line);
      break;
    }
}


int SystemMatrix::sample_calc(std::ifstream& data_file, SM_management sm_management, std::vector<double> &ev_RAM, char *SM_filename, char* results_filename, const int nb_vox[], const double voxel_length[], const double corner[], Camera &camera, data_format f)
// calculation of the system matrix for one sample, with all its partitions
{
  int ok=1;
  Event ev(f);
  int count_p=0, part_id=0;
  int total_read=0, good_events=0;
  int camera_fail=0, energy_fail=0, volume_fail=0, not_compton=0, detection_problem=0;
  std::ofstream FSM;
  char filename[100];
  //int volume_size=nb_vox[0]*nb_vox[1]*nb_vox[2];
 
  if ((sm_management==EV_CALC_DISK)||(sm_management==EV_CALC_RAM))
    std::cout << "... Calculation of the preprocessed events file ...\n";
  if (sm_management==SM_CALC)
    std::cout << "... Calculation of the system matrices ...\n";

  Image line(nb_vox,voxel_length, corner);
  Image lambda(nb_vox,voxel_length, corner);
  while (!data_file.eof() && ((presel==0 && total_read<cps) || (presel==1 && good_events<cps)))
    {
      if (sm_management==SM_CALC)
	{
	  sprintf(filename,"%s.SM%d.bin",SM_filename,part_id);
	  FSM.open(filename,std::ios::out| std::ios::binary);
	  if (!FSM.is_open())
	    {
	      std::cout << "File " << filename << " could not be opened at writting\n" ;
	      return -1;
	    }
	}
      if (sm_management==EV_CALC_DISK)
	{
	  sprintf(filename, "%s",SM_filename);
	  FSM.open(filename);
	  if (!FSM)
	    {
	      std::cout << "File " << filename << " could not be opened\n" ;
	      return -1;
	    } 
	}
      count_p=0;
      ok=0;

      while (ok>=0  && count_p<cpp && ((presel==0 && total_read<cps)||(presel==1 && good_events<cps)))
	{ 
	  ok=ev.read_event(data_file);
	  if (ok<0) break;
	  /*std::cout << " New event\n";
	    std::cout << "ID " << ev.get_id()<< '\n';
	    std::cout << "Energies: Ee= " << ev.get_Ee() << " dEe=" << ev.get_dEe()<< " Eg=" << ev.get_Eg()<< " dEg=" << ev.get_dEg() << '\n';
	    std::cout << "V1(" << ev.get_V1()[0] << ',' << ev.get_V1()[1] << ','<<ev.get_V1()[2] << ")\n";
	    std::cout << "V2(" << ev.get_V2()[0] << ',' << ev.get_V2()[1] << ','<<ev.get_V2()[2] << ")\n";
	    std::cout << "ok "<< ok << '\n';
	  */

	  total_read++;
	  if (ok>0)
	    {
	      not_compton++;
	      continue;
	    }
	  ev.set_camera(camera);
	  if (flag_energy==KNOWN)
	    ev.set_energy_parameters(Etot);

	  ok=ev.is_valid();

	  if (ok==3) 
	    {
	      camera_fail++;
	      continue;
	    }
	  if (ok>0)
	    {
	      detection_problem++;
	      continue;
	    }
	  if ((flag_energy==RANGE) && ((ev.get_Eg()+ev.get_Ee()<Emin)||(ev.get_Eg()+ev.get_Ee()>Emax))) 
	    {
	      energy_fail++;
	      continue;
	    }
	  if (!line_calc( ev, line))
	    {
	      volume_fail++;
	      continue;
	    }
	  //store line, calcul lambda0

	  switch (sm_management)
	    {
	    case SM_CALC: 
	      line.writeToFile(FSM);
	      break;
	    case EV_CALC_DISK:
	      ev.store_post_treated(FSM);
	      break;
	    case EV_CALC_RAM:
	      ev.push_post_treated(ev_RAM);
	    }
	  for (int i=0; i<line.NbVoxels; i++)
	    lambda.Value[i] += line.Value[i];
	  count_p++;
	  good_events++;
	}
      if (sm_management==SM_CALC) 
	std::cout << "Partition " << part_id << ": wrote " <<count_p<<" lines\n";
      part_id++;	  
      FSM.close();
    }
  std::cout << "Total read: " << total_read <<'\n';
  std::cout << "Not Compton: " << not_compton <<'\n';  
  std::cout << "Incorrect detection or cos(beta)<-1: " << detection_problem <<'\n';  
  std::cout << "Not in camera: " << camera_fail <<'\n';
  std::cout << "Not in energy range: " << energy_fail <<'\n';
  std::cout << "Not in volume: " << volume_fail <<'\n';
  std::cout << "Total stored: "<< good_events <<'\n';
  /*if (sm_management==EV_CALC_RAM) 
    {
      int length=ev_RAM.size();
      std::cout << "Size of the vector of events: "<< length << ", divided by 12: "<< double(length)/12. <<'\n';
      if (length!=12*good_events) 
	{
	  std::cout << "Error in the calculation of the vector of events\n";
	  return -1;
	}
    }
  */
  char iter0_results_file[100];
  sprintf(iter0_results_file, "%s.r0.bin",results_filename);
  if (lambda.writeFile(iter0_results_file)!=0) return -2;
  return data_file.eof() ? 1:0;
}



void MLEM::set_hodoscope(int flag_hodoscope, double beam_pt[], double beam_dir[], int beam_sig, int beam_width,Hodoscope_management beam_incl_type, int bfiter, int bnbiter)
{
  flag_hodo=flag_hodoscope;
  for (int i=0; i<3; i++)
    {
      beam_point[i]=beam_pt[i];
      beam_direction[i]=beam_dir[i];
    }
  beam_sigma=beam_sig;
  beam_width_factor=beam_width;
  beam_mng=beam_incl_type;
  beam_first_iter=bfiter;
  beam_nb_iter=bnbiter;
}

const void MLEM::hodoscope_matrix(Image &BM)
{
  Point Oj;
  double vec[3];
  double norm2=beam_direction[0]*beam_direction[0]+beam_direction[1]*beam_direction[1]+beam_direction[2]*beam_direction[2];
  double facteur_gauss=1./(sqrt(2*M_PI)*beam_sigma);
  double distance;
  int i,j,k;
  
  for(k=0; k<BM.DimInVoxels[2]; k++)
    {
      Oj[2]=BM.Corner[2]+(k+0.5)*BM.VoxelSize[2];
      for(j=0; j<BM.DimInVoxels[1]; j++)
        {
	  Oj[1]=BM.Corner[1]+(j+0.5)*BM.VoxelSize[1];
	  
	  vec[0]=(Oj[1]-beam_point[1])*beam_direction[2]-(Oj[2]-beam_point[2])*beam_direction[1]; 	  
	  for(i=0; i<BM.DimInVoxels[0]; i++)
            {
	      //int temp=(i+j*nb_vox[0]+k*nb_vox[0]*nb_vox[1]);
	      Oj[0]=BM.Corner[0]+(i+0.5)*BM.VoxelSize[0];
	      vec[1]=(Oj[2]-beam_point[2])*beam_direction[0]-(Oj[0]-beam_point[0])*beam_direction[2];
	      vec[2]=(Oj[0]-beam_point[0])*beam_direction[1]-(Oj[1]-beam_point[1])*beam_direction[0];
	      distance=sqrt((vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])/norm2);
	      if (distance<=(beam_width_factor*beam_sigma))
		{
		  BM.Value[BM.index_3Dto1D(i,j,k)]=facteur_gauss*exp(-0.5*distance*distance/(beam_sigma*beam_sigma));
		}
            }
        }
    }
}


/*bool MLEM::run_STORE( const char *results_file,  const char *SM_file,const std::vector<double> &sens, Image &BM, Camera &camera)
{
  //Lit une partition à la fois. Gain en temps sur le cluster ? 
  //  Cette fonction n'est pas appelée ni entretenue 
  // System matrix 
  std::ifstream FSM;
  char filename[100];
  std::vector<double> SM_vector;
  int length;
  int nb_lines=0;

  // Iterations
  int iter=first_iter; 
  bool fin=false;
  int part=0;
  Image BackProj(nb_vox,voxel_length, corner);
  double ForwardProj=0;
  int count=0, i, j, k;
  double mu;

  // read file  iter first_iter-1 ---------------------------------
  Image lambda(nb_vox,voxel_length, corner);
  char iter_results_file[100];
  sprintf(iter_results_file, "%s.r%d.bin",results_file, first_iter-1);
  if (lambda.readFile(iter_results_file)!=0) return false;



  while (iter-first_iter<nb_iter)
    {
      std::cout << "\n Iteration " << iter << '\n';
      part=0;
      fin=false;


      while (true)
	{
	  // Load the SM for the partition
	  sprintf(filename,"%s.SM%d.bin",SM_file,part);
	  FSM.open(filename,std::ios::in| std::ios::binary);  
	  if (!FSM.is_open()) break;
	  std::cout <<"Processing matrix "<<filename<< "\n";
	  FSM.seekg (0, FSM.end);
	  length = FSM.tellg();
	  SM_vector.resize(double(length)/sizeof(double),0);
	  nb_lines=SM_vector.size()/lambda.NbVoxels;
	  //std::cout << "Size of SM_value "<< SM_vector.size()<< '\n';
	  //std::cout << "Number of lines: "<< nb_lines<<'\n';
	  FSM.seekg (0, FSM.beg);
	  FSM.read((char *)(&(SM_vector[0])), length);
	  if (!FSM) 
	    {
	      std::cout << "Error in " << filename << ": only "<< FSM.gcount()<< " items could be read.\n" ;
	      return 1;
	    }
	  FSM.close();

	  for(i=0; i<nb_lines; i++)
	    {
	      // calcul projection
	      ForwardProj=0;
	      if (flag_hodo && iter>=beam_first_iter && iter-beam_first_iter<beam_nb_iter)
		{
		  mu=(double(iter-beam_first_iter))/(beam_nb_iter-1.);
		  j=i*lambda.NbVoxels;
		  for ( k=0; k<lambda.NbVoxels; k++)
		    {
		      switch (beam_mng) 
			{
			case CONSTANT :
			  SM_vector[j] *=BM.Value[k];
			  break;
			case LINEAR :
			  SM_vector[j] *=(1.-mu+mu*BM.Value[k]);
			  //if (k==0) std::cout<< "mu=" << mu<< '\n';
			  break;
			case FORCE :
			  SM_vector[j] =(1.-mu)*SM_vector[j]+mu*BM.Value[k];
			  break;
			case FORCEINV :
			  SM_vector[j] =mu*SM_vector[j]+(1.-mu)*BM.Value[k];
			  break;
			case ONLYHODO :
			  SM_vector[j] =BM.Value[k];
			}
		      ForwardProj += SM_vector[j]*lambda.Value[k];
		      j++;
		    }
		}
	      else
		{
		  j=i*lambda.NbVoxels;    
		  for (int k=0; k<lambda.NbVoxels; k++)
		    {
		      ForwardProj += SM_vector[j]*lambda.Value[k];
		      j++;
		    }
		}
	      // ajout terme correctif
	      if ( ForwardProj==0)
		{ 
		  //std::cout << "Error: ForwardProj is zero \n";
		  continue;
		}
	      j=i*lambda.NbVoxels;  
	      for (k=0; k<lambda.NbVoxels; k++)
		{
		  BackProj.Value[k] += SM_vector[j]/ForwardProj;
		  j++;
		}
	    }
	  part++;
	}
      // multipl par lambda et sj
      for (k=0; k<lambda.NbVoxels; k++)
	lambda.Value[k] =lambda.Value[k]*BackProj.Value[k]/sens[k];
      BackProj.initialize(0);
      sprintf(iter_results_file, "%s.r%d.bin",results_file, iter);
      lambda.writeFile(iter_results_file);
      iter++;

    }

}
*/

bool MLEM::run( const char *results_file,  const char *SM_file, const std::vector<double> &ev_RAM, const std::vector<double> &sens, Image &BM, Camera &camera, data_format f)
{
  // read file  iter first_iter-1 ---------------------------------

  Image lambda(nb_vox,voxel_length, corner);
  char iter_results_file[100];
  sprintf(iter_results_file, "%s.r%d.bin",results_file, first_iter-1);
  if (lambda.readFile(iter_results_file)!=0) return false;

  SystemMatrix SM(0,0,0,algo, model, width_factor,spatial_uncertainty);
  Image line(nb_vox,voxel_length,corner);
  //Event ev(camera);
  Event ev(f);
  std::ifstream FSM;
  char filename[100];
  int iter=first_iter, part=0, k; 
  bool fin=false;  // weight_calculated=false;
  Image BackProj(nb_vox,voxel_length, corner);
  double ForwardProj=0;
  //std::vector<double> weight;
  clock_t t;
  while (iter-first_iter<nb_iter)
    {
      t=clock();
      std::cout << "Iteration " << iter << '\n';
      part=0;
      // open events file or SM file --------------------------------
      if ((mng==EV_PRECALC_DISK) || (mng==EV_CALC_DISK)) 
	{
	  FSM.open(SM_file);
	  if (!FSM)
	    {
	      std::cout << "File " << SM_file << " could not be opened\n" ;
	      return false;
	    }
	} 
      if ((mng==SM_PRECALC) || (mng==SM_CALC)) 
	{
	  sprintf(filename,"%s.SM0.bin",SM_file);
	  FSM.open(filename,std::ios::in| std::ios::binary);
	  if (!FSM)
	    {
	      std::cout << "File " << filename << " could not be opened\n" ;
	      return false;
	    } 
	}
      // run iteration ---------------------------------
      fin=false;
      int count =0;
      while (!fin)
	{
	  if ((mng==EV_PRECALC_DISK) || (mng==EV_CALC_DISK)) 
	    {
	      //lecture d'un event
	      fin=ev.load_post_treated(FSM);
	      //calcul d'une ligne
	      if (fin) continue;
	      if (SM.line_calc(ev, line)==0)
		std::cout << "Empty SM line - not  usual\n";
	      count++;
	    }
	  if ((mng==SM_PRECALC) || (mng==SM_CALC)) 
	    {
	      // si fin fichier, close et open suivant
	      if (FSM.eof())
		{
		  FSM.close();
		  part++;
		  sprintf(filename,"%s.SM%d.bin",SM_file,part);
		  FSM.open(filename,std::ios::in| std::ios::binary);
		}
	      //lecture d'une ligne
	      if (FSM) line.readFromFile(FSM);
	      else fin=true;
	      if (FSM.eof()) continue;
	      count++;
	    } 
	  if ((mng==EV_PRECALC_RAM) || (mng==EV_CALC_RAM)) 
	    {
	      //lecture d'un event
	      //std::cout << "Reads event " <<ev_RAM.size()<<'\n' ;
	      fin=!ev.pop_post_treated(ev_RAM,count);
	      if (fin) continue;
	      //std::cout << " New event\n";
	      ev.set_camera(camera);
	      //ok=ev.is_valid();
	      if (SM.line_calc(ev, line)==0)
		std::cout << "Empty SM line - not  usual\n";
	      count++;
	    }
	  // calcul projection
	  ForwardProj=0; 
	  /*if (!weight_calculated)
	    {
	      double sum=0;
	      for (k=0; k<lambda.NbVoxels; k++)
		sum +=line.Value[k];
	      //if ((mng==SM_PRECALC) || (mng==SM_CALC))
	      weight.push_back(sum);
	      //else
	      //weight.push_back(sum/tan(ev.get_beta()/2));
	      //std::cout<< sum << "  " << weight[count-1]<< '\n';
	      } 
	  */
	  if (flag_hodo && iter>=beam_first_iter && iter-beam_first_iter<beam_nb_iter)
	    {
	      double mu=(double(iter-beam_first_iter))/(beam_nb_iter-1.);
	      switch (beam_mng) 
		{
		case CONSTANT :
		  for (k=0; k<lambda.NbVoxels; k++)
		    line.Value[k] *=BM.Value[k];
		  break;
		case LINEAR :
		  for (k=0; k<lambda.NbVoxels; k++)
		    line.Value[k] *=(1.-mu+mu*BM.Value[k]);
		  break;
		case FORCE :
		  for (k=0; k<lambda.NbVoxels; k++)
		    line.Value[k] =(1.-mu)*line.Value[k]+mu*BM.Value[k];
		  break;
		case FORCEINV :
		  for (k=0; k<lambda.NbVoxels; k++)
		    line.Value[k] =mu*line.Value[k]+(1.-mu)*BM.Value[k];
		  break;
		case ONLYHODO :
		  for (k=0; k<lambda.NbVoxels; k++)
		    line.Value[k] =BM.Value[k];
		  break;
		case ALTERNATE :
		  if ((iter-first_iter)/2!=double(iter-first_iter)/2)
		    for (k=0; k<lambda.NbVoxels; k++)
		      line.Value[k] =BM.Value[k];
		}
	    }
	  for (k=0; k<lambda.NbVoxels; k++)
	    ForwardProj += line.Value[k]*lambda.Value[k];
	  // ajout terme correctif
	  if ( ForwardProj==0)
	    { 
	      //std::cout << "Error: ForwardProj is zero \n";
	      continue;
	    }
	  for (k=0; k<lambda.NbVoxels; k++)
	    BackProj.Value[k] += line.Value[k]/ForwardProj;   //*weight[count-1];
	}
      //std::cout << "count= "<< count<< '\n';
      if (mng!= EV_CALC_RAM) FSM.close();
      // multipl par lambda et sj
      for (k=0; k<lambda.NbVoxels; k++)
	lambda.Value[k] =lambda.Value[k]*BackProj.Value[k]/sens[k];
      BackProj.initialize(0);
      t =clock()-t;
      std::cout<< "Time spent in iteration: "<< ((float)t)/CLOCKS_PER_SEC<<" seconds \n";
      sprintf(iter_results_file, "%s.r%d.bin",results_file, iter);
      lambda.writeFile(iter_results_file);
      // weight_calculated=true;
      iter++;
    }
  return true;
}
