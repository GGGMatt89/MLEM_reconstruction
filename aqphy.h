#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "image.h"
#include "bgeot_ftool.h"

/* Classe qui contient les paramètres de la caméra et des évenements */

typedef enum { MEGALIB, IPNL, WP3 } data_format;  
typedef enum { REAL_VALUE, STRING_VALUE, ARRAY_VALUE } line_type;



class line_content {
 public:
  line_type lt;
  double real_value;
  std::string string_value;
  std::vector<double> array_value;
  
 line_content() : real_value(0), string_value(""), array_value(0) {}
  void clear_value()
  {
    real_value=0;
    string_value="";
    array_value.clear();
  };
  int read_line(std::ifstream &, std::string , std::string , std::string);
};

class Camera
{

 private:
  data_format format;
  Point layer_dim; // layer as a box
  std::vector<std::string> layer_name;  // names of the layers
  std::vector<Point> layer_tab;    // vector containing the centres
  Point sca_dim;   // scatterer as a box
  Point sca_centre;             // for is_valid(), hit in scatterer
  Point abs_dim;  // abs as a box
  Point abs_centre;             // for is_valid(), hit in absorber
  std::string abs_name;         // name of the absorber
  Point Ox;                   // norm=1, || to a layer edge
  Point Oy;                   // norm=1, || to another layer edge
  Point Oz;                   // norm=1, orthogonal to layers
  int layer_vox[3];
  int abs_vox[3];               // norm=1, orthogonal to layers
  int layer_vox_sampling[3];
  int abs_vox_sampling[3];

  //double max;
  //double min;

  int read_MEGA(std::string);
  int read_IPNL(std::string);
  int read_WP3(std::string);
  //void max_min();

 public:

  Camera() {}
 Camera(data_format df) : format(df) {}
  void set_format(data_format df) {format=df; }
  void set_layer_dim(double length, double width, double thick) { layer_dim[0]=length; layer_dim[1]=width;  layer_dim[2]=thick;}
  void set_layer_vox(double length, double width, double thick) { layer_vox[0]=length; layer_vox[1]=width;  layer_vox[2]=thick;}
  void set_layer_vox_sampling(double length, double width, double thick) { layer_vox_sampling[0]=length; layer_vox_sampling[1]=width;  layer_vox_sampling[2]=thick;}
  //void set_sca_dim() { sca_dim[0]=layer_dim[0]; sca_dim[1]=layer_dim[1];  sca_dim[2]=max-min+layer_dim[2];}
  void set_abs_dim(double length, double width, double thick) { abs_dim[0]=length; abs_dim[1]=width;  abs_dim[2]=thick;}
  void set_abs_vox(double length, double width, double thick) { abs_vox[0]=length; abs_vox[1]=width;  abs_vox[2]=thick;}
  void set_abs_vox_sampling(double length, double width, double thick) { abs_vox_sampling[0]=length; abs_vox_sampling[1]=width;  abs_vox_sampling[2]=thick;}
  bool set_sca_box() ;
  //void set_min(double val) { min=val; }
  //void set_max(double val) { max=val; }
  void add_sca_layer(Point &, std::string &);
  void add_abs(Point &, std::string &);
  bool set_frame(Point &ox, Point &oy, Point &oz);

  data_format get_format(void ) { return format;}
  const Point& get_layer_dim() { return layer_dim;}
  const Point& get_sca_dim() {return sca_dim; }
  const Point& get_abs_dim() { return abs_dim;}
  const Point& get_sca_centre() {return sca_centre;}
  const Point& get_abs_centre() {return abs_centre;}
  const void get_layer_vox(int a[]) 
  { 
    a[0]=layer_vox[0]; a[1]=layer_vox[1]; a[2]=layer_vox[2];
  }
  const void get_abs_vox(int a[]) 
  { 
    a[0]=abs_vox[0]; a[1]=abs_vox[1]; a[2]=abs_vox[2];
  }
  const std::string& get_layer_name( int i) { return layer_name[i] ;}
  const std::string& get_abs_name( void) { return abs_name; }
  int get_layer_tab_size() { return layer_tab.size(); }
  const Point& get_sca_layer(int i) {return layer_tab[i];}
  //const Point& get_sca_layer(const int i) {return sca_tab[i];}
  const Point& get_abs(void) {return abs_centre ; } 
  //const Point& get_abs(void) {return abs_centre ; } 
  const Point& get_Ox(void) {return Ox; }
  const Point& get_Oy(void) {return Oy; }
  const Point& get_Oz(void) {return Oz; }
  const int layer_coord2index_3D(Point &V1, int ind_V1[]);
  const bool layer_sample(Point &V1, std::vector<Point> &S);
  const int abs_coord2index_3D(Point &V2, int ind_V2[]);
  const bool abs_sample(Point V2, std::vector<Point> &S);

  int read_attributes(std::string);
};

class Event
{
 private:

  data_format format;
  int id;
  double Eg;
  double Ee;
  double dEg;
  double dEe;
  double cosbeta;
  double beta;    // Compton angle
  //double alpha;   // axis polar angle
  //double delta;   // axis azimuth
  double sigma_beta;
  double P;
  double K;
  Point V1;
  Point V2;

  Camera *Cc;

 public:

 Event() : id(-1), Eg(0), Ee(0), dEg(0), dEe(0), Cc(NULL) {};
  //Event(Camera &cp) : id(-1), Eg(0), Ee(0), dEg(0), dEe(0), Cc(&cp), format(cp.get_format()) { };
 Event(data_format f) : id(-1), Eg(0), Ee(0), dEg(0), dEe(0), format(f) { };

  //acces to attributes
  void set_format(data_format f) {  format = f; }
  data_format get_format() {return format; }
  void set_id(int i) { id=i; }
  int get_id(void) {return id; }
  void set_energies(double ee, double eg, double dee, double deg) { Eg=eg; Ee=ee; dEg=deg; dEe=dee;}
  void set_energy_parameters(double Etot);
  double get_Eg(void) {return Eg;}
  double get_Ee(void) {return Ee;}
  double get_dEg(void) {return dEg;}
  double get_dEe(void) {return dEe;}
  void set_V1(Point &V) {V1[0]=V[0]; V1[1]=V[1]; V1[2]=V[2]; }
  Point & get_V1(void) {return V1; }
  void set_V2(Point &V) {V2[0]=V[0]; V2[1]=V[1]; V2[2]=V[2]; }
  Point & get_V2(void) {return V2; }
  void get_cone_axis(Point & cone_axis);
  double get_sigma_beta() { return sigma_beta; }
  double get_beta() { return beta; }
  double get_alpha();
  double get_K() { return K; }
  double get_P() { return P; }

  void set_camera(Camera& cp) { Cc=&cp; }
  //Camera& get_camera(void) { return *Cc ;}
  const Camera& get_camera(void) { return *Cc ;}
  const void get_camera_frame( Point &Ox, Point &Oy, Point &Oz) 
    { 
      Ox = Cc->get_Ox() ;
      Oy = Cc->get_Oy() ;
      Oz = Cc->get_Oz() ;
    }
  const Point & get_normal_to_camera() {return Cc->get_Oz() ;}

  // initialisation
  int read_event(std::ifstream& file);
  int read_MEGA(std::ifstream& file);
  int read_IPNL(std::ifstream& file);
  int read_WP3(std::ifstream& file);

  // analyse and storage
  int is_valid();
  //int coord2index_3D(int ind_V1[],int ind_V2[]);
  void store_post_treated(std::ofstream& file);
  bool load_post_treated(std::ifstream& file);
  bool load_post_treated_all(char filename[], std::vector<double>& ev_RAM);
  void push_post_treated(std::vector<double>& ev_RAM);
  bool pop_post_treated(const std::vector<double>& ev_RAM, int count);
};

   
