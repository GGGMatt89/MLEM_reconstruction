#pragma once

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>


class Point : public std::vector<double>
{
 public:

  Point() : std::vector<double>(3) {}
 Point(int val) : std::vector<double>(3,val) {}
 Point(const std::vector<double> &my_point) : std::vector<double>(my_point) 
    {
      assert(size() == 3);
    }

  double norm2() {return sqrt((*this)[0]*(*this)[0]+(*this)[1]*(*this)[1]+(*this)[2]*(*this)[2]);}
  const double dot_product(const Point &P) { return (*this)[0]*P[0]+(*this)[1]*P[1]+(*this)[2]*P[2]; }
  void get_local_coord(Point &P_loc, const Point &O, const Point &Ox, const Point &Oy, const Point &Oz)
{
  P_loc[0]=((*this)[0]-O[0])*Ox[0]+((*this)[1]-O[1])*Ox[1]+((*this)[2]-O[2])*Ox[2] ;
  P_loc[1]=((*this)[0]-O[0])*Oy[0]+((*this)[1]-O[1])*Oy[1]+((*this)[2]-O[2])*Oy[2] ; 
  P_loc[2]=((*this)[0]-O[0])*Oz[0]+((*this)[1]-O[1])*Oz[1]+((*this)[2]-O[2])*Oz[2] ;
}

  bool inVolume(const std::vector<double>& corner, const std::vector<double>& volume_dim);

  bool inDetector(const std::vector<double>& dim_detector, const std::vector<double>& center_detector);
};



class Image
{

 public:

  int DimInVoxels[3]; // image dimensions in voxels
  double DimInCm[3]; // image dimensions in cm
  double VoxelSize[3]; // voxel dimensions in cm
  double Corner[3]; // coordinates of the corner bottom left, front
  int NbVoxels;
  std::vector<double> Value;


  //Les constructeurs
  Image() : NbVoxels(0) {}
  Image(const int* image_size, const double* voxel_size, const double* corner);

  // Methodes qui retournent des indices 1D ou 3D dans l'image
  bool index_1Dto3D(int index_voxel, int& i, int& j, int& k) const;
  int index_3Dto1D(int i,int j, int k) const;
  bool coord2index_3D(double x, double y, double z, int& i, int& j, int& k) const;
  int coord2index_1D(double x, double y, double z) const;
	
  // Calculs sur l'image	
  void initialize(double value) ;
  void multiply(const std::vector<double> mult, const std::vector<double> divisor);
  void setIntensity(int index_voxel, double value);//Intensité du voxel d'indice_voxel = value
  double Maximum() const;//Retourne la valeur maximale des intensitées	

  // lecture/ecriture dans un fichier
  int readFile(char *image_file);
  const int writeFile(char *image_file);
  bool readFromFile(std::ifstream&  image_file);  // returns eof
  const void writeToFile(std::ofstream& image_file);


};
