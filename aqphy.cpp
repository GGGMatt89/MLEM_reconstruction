#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "aqphy.h"


int line_content::read_line(std::ifstream &file, std::string prefix, std::string postfix, std::string comment)
{
  std::string word;
  bool found_prefix=false, found=false; 
  int i, pos;
  std::istringstream line, line1, line2;
  double val;

  clear_value();
  while(!file.eof() && !found)
    {
      file >>word;
      if (found_prefix && (word.find(postfix)!=std::string::npos))
	found=true;
      found_prefix=(word.compare(prefix) ==0);
      /*
      if (found_prefix)
	std::cout << word << " : prefix\n";
      else if (word.find(postfix)!=std::string::npos)
      	std::cout << word << " : postfix\n";
      else
        std::cout << word << " : rien\n"; 
      */
    }
  if (file.eof())
    {
      std::cout<< comment << " not found\n";
      return 1;
    }
  if (lt==STRING_VALUE)
    {
      string_value=word;
    }
  else if (lt==REAL_VALUE)
    {
      if( (pos=word.find_first_of("+-0123456789."))==std::string::npos)
	getline(file,word);
      if( (pos=word.find_first_of("+-0123456789."))==std::string::npos)
	{
	  std::cout<< comment << "real value not found\n";
	  return 1;
	}
      word.erase(0,pos);
      line.str(word);
      if (!( line>> real_value))
	{
	  std::cout<< comment << " real value not found\n";
	  return 1;
	} 
    }
  else
    {
      i=-1;
      if( (pos=word.find_first_of("+-0123456789."))!=std::string::npos)
	{
	  word.erase(0,pos);
	  line1.str(word);
	  if ( line1>> val)
	    {
	      i++;
	      array_value.push_back (val);
	    };
	}
      getline(file,word);
      if(( (pos=word.find_first_of("+-0123456789."))==std::string::npos) && (i<0))
	{
	  std::cout<< comment << " array value not found\n";
	  return 1;
	} 
      word.erase(0,pos);
      line2.str(word);
      //std::cout  << "ligne "<< line2.str()  << '\n';
      while (line2 >> val)
	{
	  array_value.push_back (val);
	  i++;
	  //std::cout <<line2.str() << "val " << i << " " << val << '\n';
	}
      if (i<0) 
	{
	  std::cout<< comment << " array value not found\n";
	  return 1;
	}
    }
  return 0;
}


void Camera::add_sca_layer(Point &P, std::string & s)
{
  layer_tab.push_back(P);
  layer_name.push_back(s);
}

void Camera::add_abs(Point & P, std::string & s)
{
  abs_centre[0]=P[0];
  abs_centre[1]=P[1];
  abs_centre[2]=P[2];
  abs_name=s;
}

bool Camera::set_frame(Point &ox, Point &oy, Point &oz) 
  {
    double n=ox.norm2();
    if (n==0) return false;
    Ox[0]=ox[0]/n; 
    Ox[1]=ox[1]/n; 
    Ox[2]=ox[2]/n;  
    n=oy.norm2(); 
    if (n==0) return false;
    Oy[0]=oy[0]/n; 
    Oy[1]=oy[1]/n; 
    Oy[2]=oy[2]/n;  
    n=oz.norm2(); 
    if (n==0) return false;
    Oz[0]=oz[0]/n; 
    Oz[1]=oz[1]/n; 
    Oz[2]=oz[2]/n;  
    return true;   
  }

bool Camera::set_sca_box() 
{ 
  int i,i_min=0, i_max=0;
  double mu, mu_min=0, mu_max=0;
  if (layer_tab.size()==0) return false;
  for (i=1; i<layer_tab.size();i++)
    {
      if (Oz[0]!=0)
	mu=(layer_tab[i][0]-layer_tab[0][0])/Oz[0];
      else
	if (Oz[1]!=0)
	  mu=(layer_tab[i][1]-layer_tab[0][1])/Oz[1];
	else
	  mu=(layer_tab[i][2]-layer_tab[0][2])/Oz[2];
      if (mu<mu_min)
	{
	  mu_min=mu;
	  i_min=i;
	}
      if (mu>mu_max)
	{
	  mu_max=mu;
	  i_max=i;
	}
    }
  sca_centre[0]=(layer_tab[i_min][0]+layer_tab[i_max][0])/2.0;
  sca_centre[1]=(layer_tab[i_min][1]+layer_tab[i_max][1])/2.0;
  sca_centre[2]=(layer_tab[i_min][2]+layer_tab[i_max][2])/2.0;
  sca_dim[0]=layer_dim[0];
  sca_dim[1]=layer_dim[1];
  sca_dim[2]=mu_max-mu_min+layer_dim[2];
  /*sca_centre[0]=sca_tab[0][0]; 
    sca_centre[1]=sca_tab[0][1]; 
    sca_centre[2]=(sca_tab[0][2]+sca_tab[sca_tab.size()-1][2])/2; */
  return true;
}

int Camera::read_attributes(std::string filename)
{
  int ok=0;
  switch (format)
    {
    case MEGALIB:
      ok=read_MEGA(filename);
      break;
    case IPNL:
      ok=read_IPNL(filename);
      break;
    case WP3:
      ok=read_WP3(filename);
    }
  return ok ;
}

int Camera::read_MEGA(std::string filename)
{
  std::string line(""), geo_file, geo_path(""), sca_file , abs_file;
  std::string sca_path, dim_char;
  std::ifstream sca,abs;
  int pos;
  std::vector<double> dimension;

  // Scan for the name of .setup.geo files
  std::ifstream F;
  line_content PARAM;
  F.open(filename.c_str());

  if(!F)
    {
      std::cout<<"Error, could not open " << filename <<'\n';
      return 1;
    } 
  PARAM.lt=STRING_VALUE;
  if (PARAM.read_line(F, "Geometry",".geo.setup", "Geometry name "))
    {
      return 1;
    }
  geo_file=PARAM.string_value;
 
  pos=geo_file.rfind("/");
  geo_path=geo_file.substr(0,pos+1);
  // geo_file.erase(0,pos+1);
  std::cout<<"Geometry file: "<<geo_file<<'\n';
  std::cout<<"Geometry path: "<<geo_path<<'\n';

  // Scan for the name of .geo files

  std::ifstream Fg;
  Fg.open(geo_file.c_str());

  if(!Fg)
    {
      std::cout<<"Error, could not open " << geo_file <<'\n';
      return 1;
    } 
  PARAM.lt=STRING_VALUE;
  if (PARAM.read_line(Fg, "Include","Detector.geo", "Detector name "))
    {
      return 1;
    }
  std::string sca_layer=PARAM.string_value;
  Fg.seekg(0);
  PARAM.lt=STRING_VALUE;
  if (PARAM.read_line(Fg, "Include", "Absorber.geo", "Absorber name "))
    {
      return 1;
    }
  std::string abs_layer=PARAM.string_value;
  sca_file = geo_path+sca_layer;
  sca_layer.erase(sca_layer.find(".geo"), std::string::npos);
  std::cout<<"Scatterer file: "<<sca_file<<'\n';
  abs_file=geo_path+abs_layer;
  abs_layer.erase(abs_layer.find(".geo"),std::string::npos);
  std::cout<<"Absorber file: "<<abs_file<<"\n";
  Fg.close();

  // Scan for dimensions of a layer from the scatterer 

  sca.open(sca_file.c_str());
  if(!sca)
    {
      std::cout<<"Error, could not open " << sca_file <<'\n';
      return 1;
    } 

  /* PARAM.lt=ARRAY_VALUE;
  F.seekg(0);
  read_line(F, "Shape", "", PARAM, "real size 1 ");
  for (int i=0 ;i<PARAM.array_value.size(); i++)
    std::cout <<  "real size 1 "<< PARAM.array_value[i] << '\n';
  */

  PARAM.lt=ARRAY_VALUE;
  //std::cout << "Sca layer : " << sca_layer<< '\n';
  if (PARAM.read_line(sca, sca_layer+".Shape", "", "Scatterer size "))
    {
      return 1;
    }
  if (PARAM.array_value.size()<3)
    {
      std::cout << "Not enough dimensions for the scatterer\n";
      return 1;
    }
  set_layer_dim ( 2*PARAM.array_value[0],  2*PARAM.array_value[1], 2*PARAM.array_value[2]);
  //std::cout << "Sca dim " << scatterer_length<< "  " << scatterer_depth<<'\n';
  sca.close();

  // Scan for dimensions of the absorber

  abs.open(abs_file.c_str());
  if(!abs)
    {
      std::cout<<"Error, could not open " << abs_file <<'\n';
      return 1;
    } 
  PARAM.lt=ARRAY_VALUE;
  //std::cout << "Abs layer : " << abs_layer<< '\n';
  if (PARAM.read_line(abs, abs_layer+".Shape", "", "Scatterer size "))
    {
      return 1;
    }
  if (PARAM.array_value.size()<3)
    {
      std::cout << "Not enough dimensions for the absorber\n";
      return 1;
    }
  set_abs_dim (2*PARAM.array_value[0], 2*PARAM.array_value[1],2*PARAM.array_value[2]);
  //std::cout << "Abs dim " << absorber_length<< "  " << absorber_depth<<'\n';
  abs.close();

  // Seek for the number of scatterers
  //std::cout << "Seek for the number of scatterers \n";
  //Fg.seekg(0);
  Fg.open(geo_file.c_str());
  Point P;
  std::string layer_name;
  while (!Fg.eof())
    {
      PARAM.lt=STRING_VALUE;
      if (PARAM.read_line(Fg, sca_layer+".Copy" ,"", "End scan: other scatterer layer "))
	break;
      PARAM.lt=ARRAY_VALUE;
      layer_name=PARAM.string_value;
      //Fg.seekg(0);
      //std::cout << layer_name << '\n';
      if (PARAM.read_line(Fg, layer_name+".Position" ,"", "Layer coord "))
	{
	  std::cout << "Error in " << geo_file << " at "<< layer_name+".Position" << '\n';
	  return 1;
	}
      if (PARAM.array_value.size()<3)
	{
	  std::cout << "Not enough dim for layer " << layer_name << '\n';
	  return 1;
	}
      P[0]=PARAM.array_value[0];
      P[1]=PARAM.array_value[1];
      P[2]=PARAM.array_value[2];
      add_sca_layer(P,layer_name);
    }
  Fg.close();

  // Seek for the absorber

  
  Fg.open(geo_file.c_str());  //Fg.seekg(0);
  PARAM.lt=STRING_VALUE;
  //std::cout<< "abs " <<abs_layer+".Copy\n";
  if (PARAM.read_line(Fg, abs_layer+".Copy" ,"", " Absorber "))
    {
      std::cout << "Absorber copy " << abs_layer << " not found." << '\n';
      return 1;
    }
  PARAM.lt=ARRAY_VALUE;
  layer_name=PARAM.string_value;
  if (PARAM.read_line(Fg, layer_name+".Position" ,"", "Absorber coord "))
    {
      std::cout << "Error in " << geo_file << " at "<< layer_name+".Position" << '\n';
      return 1;
    }
  if (PARAM.array_value.size()<3)
    {
      std::cout << "Not enough dim for absorber " << layer_name << '\n';
      return 1;
    }
  P[0]=PARAM.array_value[0];
  P[1]=PARAM.array_value[1];
  P[2]=PARAM.array_value[2];
  add_abs (P,layer_name);
  //std::cout << "abs " << abs_layer << " " << abs_tab[0][2];
  Fg.close();
  //max_min();
  set_sca_box();
  //set_sca_centre();
  return 0;

}

int Camera::read_IPNL(std::string filename)
{
  std::string line, scatterer_geo,absorber_geo, path;
  std::ifstream txt,geo,scatterer,absorber;
  std::string useless;
  double source_first_scatterer_distance,first_to_last_scatterer_distance,last_scatterer_absorber_distance, x,y,z;
  int nb_layer, nb_abs;

  //std::cout<<"nom: "<<filename<<std::endl;
  
  txt.open(filename.c_str());
  if(!txt)
    {
      std::cout << "File " << filename << " not found\n";
      return 1;
    }
  int i=0;
  while (!txt.eof() && (i<4))
    {
      getline(txt, line);
      i++;
    }
  if (txt.eof())
    {
      std::cout << "Error in " << filename << '\n';
      return 1;
    }
  i=0;
  while (!txt.eof() && (i<5))
    {
      txt>>useless;
      i++;
    }
  if (txt.eof())
    {
      std::cout <<"Error at source_first_scatterer_distance\n";
      return 1;
    }
  txt>>source_first_scatterer_distance;
  getline(txt,line);
  i=0;
  while (!txt.eof() && (i<7))
    {
      txt>>useless;
      i++;
    }
  if (txt.eof())
    {
      std::cout <<"Error at first_to_last_scatterer_distance\n";
      return 1;
    }
  txt>>first_to_last_scatterer_distance;
  getline(txt,line);
  i=0;
  while (!txt.eof() && (i<6))
    {
      txt>>useless;
      i++;
    }
  if (txt.eof())
    {
      std::cout <<"Error at last_scatterer_absorber_distance\n";
      return 1;
    }
  txt>>last_scatterer_absorber_distance;
  getline(txt,line);
  txt>>useless;txt>>useless;
  txt>>nb_layer;
  txt>>useless;txt>>useless;
  txt>>x; txt>>useless;
  txt>>y;txt>>useless;
  txt>>z;
  set_layer_dim(x,y,z);
  getline(txt,line);
  txt>>useless;txt>>useless;
  txt>>nb_abs;
  txt>>useless;txt>>useless;
  txt>>x ; txt>>useless;
  txt>>y;txt>>useless;
  txt>> z;
  set_abs_dim(x,y,z);

  //std::cout<<"Configuration of the camera: \n";
  Point P;
  std::string empty;

  for(i=0;i<nb_layer;i++)
    {
      P[0]=0;
      P[1]=0;
      P[2]=-(source_first_scatterer_distance+i*(first_to_last_scatterer_distance/(nb_layer-1)));
      add_sca_layer(P,empty);
     // sca_tab.push_back (P);
      //std::cout<<"layer "<<i<<" : " << sca_tab[i][2];
    }
  //std::cout<<"number of scattering layers: "<<sca_tab.size() <<'\n';

  P[0]=0;
  P[1]=0;
  P[2]=-(source_first_scatterer_distance+first_to_last_scatterer_distance+last_scatterer_absorber_distance);
  add_abs(P,empty);
  //std::cout<<"absorber : "<<abs_coord[2] <<'\n';
  //max_min();
  set_sca_box();
  //set_sca_centre();
  return 0;
}

int Camera::read_WP3(std::string filename)
{
  std::string line, scatterer_geo,absorber_geo, path;
  double Lx,Ly,Lz,Dx,Dy,Dz;
  Point P;

  std::ifstream txt ;
  txt.open(filename.c_str(), std::ifstream::in);
  if(!txt.is_open())
    {
      std::cout<<"camera not found."<<'\n';
      return 1;
    }
  getline(txt, line);
  txt >> Lx>>Ly>>Lz>>P[0]>>P[1]>>P[2];

  //cout<<layers[0]<<" "<<layers[1]<<" "<<layers[2]<<" "<<endl;
  //std::cout<<"position diffuseur: "<<P[0]<<" "<<P[1]<<" "<<P[2]<<" "<<'\n';

  set_layer_dim(Lx,Ly,Lz);
  //sca_tab.push_back(P);
  std::string empty("");
  add_sca_layer(P,empty);
  //set_min(P[2]);
  //set_max(P[2]);

  txt >> Lx>>Ly>> Lz>>P[0]>>P[1]>>P[2];
  add_abs(P,empty);
  //cout<<layers[6]<<" "<<layers[7]<<" "<<layers[8]<<" "<<endl;
  //std::cout<<"position absorbeur: "<<abs_coord[0]<<" "<<abs_coord[1]<<" "<<abs_coord[2]<<" "<<'\n';
  set_abs_dim(Lx,Ly,Lz);
  set_sca_box();
  //set_sca_centre();
  txt.close();
  return 0;
}


/*
void Camera::max_min()
{

  double m=10000,M=-10000;

  for(int i=0;i<sca_tab.size() ;i++)
    {
      m=(m<sca_tab[i][2]) ? m : sca_tab[i][2];
      M=(M>sca_tab[i][2]) ? M : sca_tab[i][2];
    }
  set_max(M);
  set_min(m);
}
*/

const int Camera::layer_coord2index_3D(Point &V1, int ind_V1[])
{
  bool inside=false;
  int i=0;
  Point v,centre;
  //std::cout<< "Coordinates for V1 ("<< V1[0] << ',' << V1[1] << ',' << V1[2] << ")\n";
  while (!inside && i<get_layer_tab_size())
    {
      centre=get_sca_layer(i);
      V1.get_local_coord(v,centre,Ox,Oy,Oz);
      inside=(fabs(v[0])<=layer_dim[0]) && (fabs(v[1])<=layer_dim[1]) && (fabs(v[2])<=layer_dim[2]);
      i++;
    }
  if (!inside) return -1;
  //std::cout<< "New coordinates for V1 ("<< v[0] << ',' << v[1] << ',' << v[2] << ")\n";
  ind_V1[0]=int(std::min(int(floor((v[0]/layer_dim[0]+0.5)*layer_vox[0])), layer_vox[0]));
  ind_V1[1]=int(std::min(int(floor((v[1]/layer_dim[1]+0.5)*layer_vox[1])), layer_vox[1]));
  ind_V1[2]=int(std::min(int(floor((v[2]/layer_dim[2]+0.5)*layer_vox[2])), layer_vox[2]));
  return i-1;
}

const bool Camera::layer_sample(Point &V1, std::vector<Point> &S)
{
  int ind_V1[3];
  //std::cout<< "Coordinates for V1 ("<< V1[0] << ',' << V1[1] << ',' << V1[2] << ")\n";
  S.clear();
  int l=layer_coord2index_3D(V1,ind_V1 );
  //std::cout<< "Index for V1 ("<< ind_V1[0] << ',' << ind_V1[1] << ',' << ind_V1[2] << ")\n";
  if (l<0) return false;
  Point corner_vox;
  double coord;
  double xi,yi,zi,xj,yj,zj;
  corner_vox[0]=(double(ind_V1[0])/layer_vox[0]-0.5)*layer_dim[0];
  corner_vox[1]=(double(ind_V1[1])/layer_vox[1]-0.5)*layer_dim[1];
  corner_vox[2]=(double(ind_V1[2])/layer_vox[2]-0.5)*layer_dim[2];
  //std::cout << "Corner ("<< corner_vox[0]<<','<< corner_vox[1]<<',' <<corner_vox[2]<<")\n";
  for (int i=0; i<layer_vox_sampling[0]; i++)
    {
      coord=corner_vox[0]+(i+1)*layer_dim[0]/layer_vox[0]/(layer_vox_sampling[0]+1);
      xi = coord*Ox[0];
      yi = coord*Oy[0];
      zi = coord*Oz[0];
      for (int j=0; j<layer_vox_sampling[1]; j++)
	{
	  coord=corner_vox[1]+(j+1)*layer_dim[1]/layer_vox[1]/(layer_vox_sampling[1]+1);
	  xj = coord*Ox[1];
	  yj = coord*Oy[1];
	  zj = coord*Oz[1];
	  for (int k=0; k<layer_vox_sampling[2]; k++)
	    {
	      coord=corner_vox[2]+(k+1)*layer_dim[2]/double(layer_vox[2])/(layer_vox_sampling[2]+1);
	      Point P(get_sca_layer(l));
	      P[0] += xi+xj+coord*Ox[2];
	      P[1] += yi+yj+coord*Oy[2];
	      P[2] += zi+zj+coord*Oz[2];
	      //std::cout << "New V1 ("<< P[0]<<','<< P[1]<<',' <<P[2]<<")\n";
	      S.push_back(P);
	    }
	}
    }
  return true;
}

const int Camera::abs_coord2index_3D(Point &V2, int ind_V2[])
{
  bool inside=false;
  int i=0;
  Point v,centre;
  //std::cout<< "Coordinates for V2 ("<< V2[0] << ',' << V2[1] << ',' << V2[2] << ")\n";
  centre=get_abs();
  V2.get_local_coord(v,centre,Ox,Oy,Oz);
  inside=(fabs(v[0])<=abs_dim[0]) && (fabs(v[1])<=abs_dim[1]) && (fabs(v[2])<=abs_dim[2]);
  if (!inside) return -1;
  //std::cout<< "New coordinates for V2 ("<< v[0] << ',' << v[1] << ',' << v[2] << ")\n";
  ind_V2[0]=int(std::min(int(floor((v[0]/abs_dim[0]+0.5)*abs_vox[0])), abs_vox[0]));
  ind_V2[1]=int(std::min(int(floor((v[1]/abs_dim[1]+0.5)*abs_vox[1])), abs_vox[1]));
  ind_V2[2]=int(std::min(int(floor((v[2]/abs_dim[2]+0.5)*abs_vox[2])), abs_vox[2]));
  return 0;
}


const bool Camera::abs_sample(Point V2, std::vector<Point> &S)
{
  int ind_V2[3];
  //std::cout<< "Coordinates for V2 ("<< V2[0] << ',' << V2[1] << ',' << V2[2] << ")\n";
  S.clear();
  int l=abs_coord2index_3D(V2,ind_V2 );
  //std::cout<< "Index for V2 ("<< ind_V2[0] << ',' << ind_V2[1] << ',' << ind_V2[2] << ") l=" << l<< "\n";
  if (l<0) return false;
  Point corner_vox;
  double coord;
  double xi,yi,zi,xj,yj,zj;
  corner_vox[0]=(double(ind_V2[0])/abs_vox[0]-0.5)*abs_dim[0];
  corner_vox[1]=(double(ind_V2[1])/abs_vox[1]-0.5)*abs_dim[1];
  corner_vox[2]=(double(ind_V2[2])/abs_vox[2]-0.5)*abs_dim[2];
  //std::cout << "Corner ("<< corner_vox[0]<<','<< corner_vox[1]<<',' <<corner_vox[2]<<")\n";
  for (int i=0; i<abs_vox_sampling[0]; i++)
    {
      //std::cout << abs_dim[2]/double(abs_vox[2])/(abs_vox_sampling[2]+1)<<'\n';
      coord=corner_vox[0]+(i+1)*abs_dim[0]/double(abs_vox[0])/(abs_vox_sampling[0]+1);
      xi = coord*Ox[0];
      yi = coord*Oy[0];
      zi = coord*Oz[0];
      for (int j=0; j<abs_vox_sampling[1]; j++)
	{
	  coord=corner_vox[1]+(j+1)*abs_dim[1]/double(abs_vox[1])/(abs_vox_sampling[1]+1);
	  xj = coord*Ox[1];
	  yj = coord*Oy[1];
	  zj = coord*Oz[1];
	  for (int k=0; k<abs_vox_sampling[2]; k++)
	    {
	      coord=corner_vox[2]+(k+1)*abs_dim[2]/double(abs_vox[2])/(abs_vox_sampling[2]+1);
	      Point P(get_abs());
	      P[0] += xi+xj+coord*Ox[2];
	      P[1] += yi+yj+coord*Oy[2];
	      P[2] += zi+zj+coord*Oz[2];
	      //std::cout << "New V2 ("<< P[0]<<','<< P[1]<<',' <<P[2]<<")\n";
	      S.push_back(P);
	    }
	}
    }
  return true;
}




void Event::set_energy_parameters(double Etot)
// default Etot=-1; known initial energy value when Etot<>Ee+Eg
{
  if (Etot<0)
    {
      cosbeta=1-(double)((Ee*511)/(Eg*(Eg+Ee)));
      if (cosbeta<-1) 
	{
	  beta=-M_PI;
	  P=0;
	  K=0;
	  sigma_beta=0;
	  return;
	}
      beta=acos(cosbeta); 
      P=Eg/(Eg+Ee);
      K=P*P*0.5*(P+(1./P)-1+cosbeta*cosbeta);	  
      /*if (format==WP3)  
	{
	double dee=(6+0.45*sqrt(Ee))/2.355;  // Xavier :(6+0.15*sqrt(Ee))/2.355;
	double deg=1.75*sqrt(eg);      // Xavier : 0.5*1.75*sqrt(Eg);
	set_energies(Ee,eg,dee,deg);
	} */
      double Tg=Ee*(2*Eg+Ee)/(Eg*Eg)*dEg;
      sigma_beta=511./((Eg+Ee)*(Eg+Ee)*sin(beta))*sqrt(dEe*dEe+Tg*Tg);
    }
  else
    {
      double eg=Etot-Ee;
      cosbeta=1-(double)((Ee*511)/(eg*Etot));
      if (cosbeta<-1) 
	{
	  beta=-M_PI;
	  P=0;
	  K=0;
	  sigma_beta=0;
	  return;
	}
      beta=acos(cosbeta); 
      P=eg/Etot;
      K=P*P*0.5*(P+(1./P)-1+cosbeta*cosbeta);	
      double Tg=(1./(eg*eg)-1./(Etot*Etot))*dEg;
      double Te=dEe/(Etot*Etot);
      sigma_beta=511*dEe/(sin(beta)*eg*eg);
    }
}

void Event::get_cone_axis(Point & cone_axis) 
{
  cone_axis[0]=(V1[0]-V2[0]); 
  cone_axis[1]=(V1[1]-V2[1]);
  cone_axis[2]=(V1[2]-V2[2]);
  double n=cone_axis.norm2();
    //norm=sqrt((V1[0]-V2[0])*(V1[0]-V2[0])+(V1[1]-V2[1])*(V1[1]-V2[1])+(V1[2]-V2[2])*(V1[2]-V2[2]));
  cone_axis[0] /=n; 
  cone_axis[1] /=n;
  cone_axis[2] /=n;
}

double Event::get_alpha()
{
  // calculation of the parameters of the axis
  Point cone_axis;
  get_cone_axis(cone_axis);
  //Point ca=Cc->get_Oz();   //Cc->get_axis();
  return acos(cone_axis.dot_product(Cc->get_Oz()));
  //cone_axis[0]*ca[0]+cone_axis[1]*ca[1]+cone_axis[2]*ca[2]);

}


int Event::read_event(std::ifstream& file)
{
  int ok=0;
  switch (format)
    {
    case MEGALIB:
      ok=read_MEGA(file);
      break;
    case IPNL:
      ok=read_IPNL(file);
      break;
    case WP3:
      ok=read_WP3(file);
    }
  // calculation of energy-related parameters
  set_energy_parameters(-1);
  return ok ;
}


int Event::read_MEGA(std::ifstream& file)
{
  /*
    return 0 if a Compton event was read
    return 1 if the event is not Compton
    return -1 if end of file reached
  */

  std::string line;
  bool new_event=false;
  double x=0;
  double eg,deg,ee,dee;
  Point V;
  do
    {
      getline(file,line);
      new_event = (line=="SE");
    }
  while (!file.eof() && !new_event);
  if (file.eof()) return -1;
  getline(file,line);
  if (line!="ET CO")  return 1;
  file>>line;
  file>> id;
//pour nouveaux fichiers .tra
  for (int i=0;i<6; i++)
    getline(file,line);
  file >> line ; // >> eg >> deg>> ee >> dee;
  if (line!="CE") 
    {
      getline(file,line);
      file >> line >> eg >> deg>> ee >> dee;
    }
  else 
    file >>  eg >> deg>> ee >> dee;
  set_energies(ee , eg, dee, deg);
  file >> line >> V[0]>> V[1]>> V[2]>> x>> x>> x;
  set_V1(V);
  file >> V[0]>> V[1]>> V[2] ;
  set_V2(V);
  getline(file,line);
  getline(file,line);

// pour anciens fichiers .tra
/*  for (int i=0;i<8; i++)
    getline(file,line);
  file >> line ; // >> eg >> deg>> ee >> dee;
  if (line!="CE") 
    {
      getline(file,line);
      file >> line >> eg >> deg>> ee >> dee;
    }
  else 
    file >>  eg >> deg>> ee >> dee;
  set_energies(ee , eg, dee, deg);
  file >> line >> V[0]>> V[1]>> V[2]>> x>> x>> x;
  set_V1(V);
std::cout<<"V1 = "<<V[0]<<" "<<V[1]<<" "<<V[2]<<std::endl;
  file >> V[0]>> V[1]>> V[2] ;
  set_V2(V);
std::cout<<"V2 = "<<V[0]<<" "<<V[1]<<" "<<V[2]<<std::endl<<std::endl;
  getline(file,line);
  getline(file,line);
  getline(file,line);*/

  return 0;
}

int Event::read_IPNL(std::ifstream& file)
{
  /*
    return 0 if a Compton event was read
    return -1 if end of file reached
  */

  std::string line;
  bool new_event=false;
  double x=0,eg,ee,deg,dee;
  Point V;
  do
    {
      getline(file,line);
      new_event = (line.compare(0,19,"---Nouvel evenement")==0);
    }
  while (!file.eof() && !new_event);
  if (file.eof()) { return -1;}
  file >> line ; file >> line ; file >> V[0] >> V[1] >> V[2];
  V[0] /=10; V[1] /=10; V[2] /=10;   // initial values en mm
  set_V1(V);
  file >> line ; file >> line ; file >> ee;    // initial values in MeV
  ee *=1000;
  file >> line ; file >> line ; file >> V[0]>> V[1]>> V[2];
  V[0] /=10; V[1] /=10; V[2] /=10;  
  set_V2(V);  
  file >> line ; file >> line ; file >> eg;
  eg *=1000;

  // ERREURS ENERGIES IPNL:
  double pce=3.65/1000;           // pair creation energy in kev
  double F=0.115;                 // Fano factor
  double ENC=600;                 // Equivalent Noise Charge

  dee=pce*sqrt(ENC*ENC+F*ee/pce);     // These MarieHelene: Eg_FWHM=2.355*formuleFano
  deg=(8./100.)/2.35*sqrt(eg);
  set_energies(ee,eg,dee,deg);
  return 0;
}

int Event::read_WP3(std::ifstream& file)
{
  /*
    return 0 if a Compton event was read
    return -1 if end of file reached
  */
  if (file.eof()) return -1;
  Point V;
  double eg,ee,deg,dee;
  std::string line;
  file>> id;
  //std::cout << "id "<< id <<'\n';
  file >> V[0]>> V[1]>> V[2] >> ee;
  set_V1(V);
  file >> V[0]>> V[1]>> V[2] >> eg;
  set_V2(V);
  dee=(6+0.45*sqrt(ee))/2.355;  // Xavier :(6+0.15*sqrt(Ee))/2.355;
  deg=1.75*sqrt(eg);      // Xavier : 0.5*1.75*sqrt(Eg);
  set_energies(ee,eg,dee,deg);
  return 0;
}

int Event::is_valid()
{
  if(cosbeta<-1)   // cos(beta)<-1;
    return 4;
  Point Ox, Oy, Oz;
  get_camera_frame(Ox,Oy,Oz);
  if ((V1[0]-V2[0])*Oz[0]+(V1[1]-V2[1])*Oz[1]+(V1[2]-V2[2])*Oz[2]==0)  // not scattered
    return 1;
  if(Eg==0)   // no second interaction
    return 2;
  Point cs=Cc->get_sca_centre();
  Point v1;
  V1.get_local_coord(v1,cs,Ox,Oy,Oz);
  Point ca=Cc->get_abs_centre();
  Point v2_s, v2_a;
  Point ds = Cc->get_sca_dim();
  Point da = Cc->get_abs_dim();
  V2.get_local_coord(v2_s,cs,Ox,Oy,Oz);
  V2.get_local_coord(v2_a,ca,Ox,Oy,Oz); 

  //if ((fabs(v1[0])>cs[0]/2) || (fabs(v1[1])>cs[1]/2) || (fabs(v1[2])>cs[2]/2) || (fabs(v2[0])>ca[0]/2) || (fabs(v2[1])>ca[1]/2) || (fabs(v2[2])>ca[2]/2))  // out of camera
  if ( !(((fabs(v1[0])<=ds[0]/2) && (fabs(v1[1])<=ds[1]/2) && (fabs(v1[2])<=ds[2]/2) && (fabs(v2_s[0])<=ds[0]/2) && (fabs(v2_s[1])<=ds[1]/2) && (fabs(v2_s[2])<=ds[2]/2)) /* V1 pas dans le diffuseur et V2 pas dans le diffuseur */
	 ||
	 ((fabs(v1[0])<=ds[0]/2) && (fabs(v1[1])<=ds[1]/2) && (fabs(v1[2])<=ds[2]/2) && (fabs(v2_a[0])<=da[0]/2) && (fabs(v2_a[1])<=da[1]/2) && (fabs(v2_a[2])<=da[2]/2)) /* V1 pas dans le diffuseur et V2 pas dans l'absorbeur */
	 ) )  // out of camera 
    return 3;
  //if(!(V1.inDetector(Cc->get_sca_dim(),Cc->get_sca_centre() ) && V2.inDetector(Cc->get_abs_dim(),Cc->get_abs()))) 
  return 0;
}

/*
int Event::coord2index_3D(int ind_V1[],int ind_V2[])
{
  bool inside=false;
  int i=0;
  Point Ox, Oy, Oz;
  get_camera_frame(Ox,Oy,Oz);
  Point centre, dim;
  int nb_vox[3];
  Point v;
  dim=Cc->get_layer_dim();
  //std::cout<< "Coordinates for V1 ("<< V1[0] << ',' << V1[1] << ',' << V1[2] << ")\n";
  while (!inside && i<Cc->get_layer_tab_size())
    {
      centre=Cc->get_sca_layer(i);
      V1.get_local_coord(v,centre,Ox,Oy,Oz);
      inside=(fabs(v[0])<=dim[0]) && (fabs(v[1])<=dim[1]) && (fabs(v[2])<=dim[2]);
      i++;
    }
  if (!inside) return -1;
  //std::cout<< "New coordinates for V1 ("<< v[0] << ',' << v[1] << ',' << v[2] << ")\n";
  Cc->get_layer_vox(nb_vox);
  ind_V1[0]=int(std::min(int(floor((v[0]/dim[0]+0.5)*nb_vox[0])), nb_vox[0]));
  ind_V1[1]=int(std::min(int(floor((v[1]/dim[1]+0.5)*nb_vox[1])), nb_vox[1]));
  ind_V1[2]=int(std::min(int(floor((v[2]/dim[2]+0.5)*nb_vox[2])), nb_vox[2]));

  dim=Cc->get_abs_dim();
  //std::cout<< "Coordinates for V2 ("<< V2[0] << ',' << V2[1] << ',' << V2[2] << ")\n";
  centre=Cc->get_abs();
  V2.get_local_coord(v,centre,Ox,Oy,Oz);
  //std::cout<< "New coordinates for V1 ("<< v[0] << ',' << v[1] << ',' << v[2] << ")\n";
  inside=(fabs(v[0])<=dim[0]) && (fabs(v[1])<=dim[1]) && (fabs(v[2])<=dim[2]);
  if (!inside) return -2;
  Cc->get_abs_vox(nb_vox);
  ind_V2[0]=int(std::min(int(floor((v[0]/dim[0]+0.5)*nb_vox[0])), nb_vox[0]));
  ind_V2[1]=int(std::min(int(floor((v[1]/dim[1]+0.5)*nb_vox[1])), nb_vox[1]));
  ind_V2[2]=int(std::min(int(floor((v[2]/dim[2]+0.5)*nb_vox[2])), nb_vox[2]));
  return(i-1);
}
*/

void Event::store_post_treated(std::ofstream& file)
{
  file << id << " : " << V1[0] << ", "<< V1[1] << ", "<< V1[2] << "; " ;
  file << V2[0] << ", "<< V2[1] << ", "<< V2[2] << "; " ;
  //file << alpha << ", " << alpha*180/M_PI<< "; ";
  file << beta << ", " << beta*180/M_PI<< "; ";
  file << sigma_beta<< ", " << sigma_beta*180/M_PI<<"; " ;
  file << P<< ", " << K <<'\n'; 

}

bool Event::load_post_treated(std::ifstream& file)
{
  double tmp;
  std::string line;
  file >> id ;
  if (file.eof()) return true;
  file>> line;
  // std::cout << "ID " << id<< '\n';
  file >> V1[0] >> line >> V1[1] >> line >> V1[2] >> line ;
  // std::cout << "V1(" << V1[0] << ',' << V1[1] << ','<<V1[2] << ")\n";
  file >> V2[0]  >> line>> V2[1] >> line>> V2[2] >> line ;
  // std::cout << "V2(" << V2[0] << ',' << V2[1] << ','<<V2[2] << ")\n";
  //file >> alpha >> line >> tmp >> line ;
  // std::cout <<"alpha="<< alpha << '\n';
  file >> beta >> line >> tmp >> line ;
  // std::cout <<"beta="<< beta << '\n';
  file >> sigma_beta >> line >> tmp >> line ;
  // std::cout <<"sigma_beta="<< sigma_beta << '\n';
  file >> P>> line >> K ;
  //std::cout <<"K="<< K << '\n';
  getline(file,line);
  return false;
}


bool Event::load_post_treated_all(char filename[], std::vector<double>& ev_RAM)
{
  std::ifstream fp;
  fp.open(filename);
  if (!fp)
    {
      std::cout<<"Error at opening of the post-treated events file " << filename<< "\n";
      return false;
    } 
  Event ev;
  bool fin=false;
  while (!fin)
    { 
      fin=ev.load_post_treated(fp);
      if (!fin)	ev.push_post_treated(ev_RAM);
    }
  fp.close();
  return true;
}


void Event::push_post_treated(std::vector<double>& ev_RAM)
{
  ev_RAM.push_back(double(id));
  ev_RAM.push_back(V1[0]);
  ev_RAM.push_back(V1[1]);
  ev_RAM.push_back(V1[2]);
  ev_RAM.push_back(V2[0]);
  ev_RAM.push_back(V2[1]);
  ev_RAM.push_back(V2[2]);
  //ev_RAM.push_back(alpha);
  ev_RAM.push_back(Eg);
  ev_RAM.push_back(Ee);
  ev_RAM.push_back(beta);
  ev_RAM.push_back(sigma_beta);
  ev_RAM.push_back(P);
  ev_RAM.push_back(K);
}

bool Event::pop_post_treated(const std::vector<double>& ev_RAM, int count)
{
  int nb_champs=13;
  if (ev_RAM.size()<count*nb_champs+nb_champs-1) 
    {
      //std::cout << "Error, not enough element in the array of events\n";
      return false;
    }
  id=ev_RAM[count*nb_champs];
  V1[0]=ev_RAM[count*nb_champs+1];
  V1[1]=ev_RAM[count*nb_champs+2];
  V1[2]=ev_RAM[count*nb_champs+3];
  V2[0]=ev_RAM[count*nb_champs+4];
  V2[1]=ev_RAM[count*nb_champs+5];
  V2[2]=ev_RAM[count*nb_champs+6];
  //alpha=ev_RAM[count*nb_champs+7];
  Eg=ev_RAM[count*nb_champs+7];
  Ee=ev_RAM[count*nb_champs+8];
  beta=ev_RAM[count*nb_champs+9];
  sigma_beta=ev_RAM[count*nb_champs+10];
  P=ev_RAM[count*nb_champs+11];
  K=ev_RAM[count*nb_champs+12];
  return true;
}
