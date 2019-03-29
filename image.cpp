#include <iostream>
#include <fstream>
#include <math.h>
#include "image.h"



bool Point::inVolume(const std::vector<double>& corner, const std::vector<double>& volume_dim)
{
  bool x=((*this)[0]>=corner[0]) && ((*this)[0]<=corner[0]+volume_dim[0]);
  bool y=((*this)[1]>=corner[1]) && ((*this)[1]<=corner[1]+volume_dim[1]);
  bool z=((*this)[2]>=corner[2]) && ((*this)[2]<=corner[2]+volume_dim[2]);
  return (x && y && z);
}

bool Point::inDetector(const std::vector<double>& dim_detector, const std::vector<double>& center_detector)
{
  bool x=fabs((*this)[0]-center_detector[0])<=dim_detector[0];
  bool y=fabs((*this)[1]-center_detector[1])<=dim_detector[1];
  bool z=fabs((*this)[2]-center_detector[2])<=dim_detector[2];
  return (x && y && z);
}


Image::Image(const int* image_size, const double* voxel_size, const double* corner)
{
  for (int i=0;i<3; i++)
    {
      DimInVoxels[i]=image_size[i];
      VoxelSize[i]=voxel_size[i];
      Corner[i]=corner[i];
      DimInCm[i]=DimInVoxels[i]*VoxelSize[i];
    }
  NbVoxels=DimInVoxels[0]*DimInVoxels[1]*DimInVoxels[2];
  Value.assign( NbVoxels,0.0);
}

/*
void Image::getVoxel_size( std::vector<double>& voxel_size) const
{  
  for (int i=0;i<3; i++)  voxel_size[i]=m_voxel_size[i];
}

void Image::getCorner( std::vector<double>& corner) const
{
  for (int i=0;i<3; i++)  corner[i]=m_corner[i];
}

void Image::getImage_size(std::vector<int>& image_size) const
{
  for (int i=0;i<3; i++)  image_size[i]=m_image_size[i];
}

double Image::getIntensity(int index_voxel) const
{
  return (m_volume[index_voxel]);
}

void Image::getVolume(std::vector<double>& volume ) const
{
   volume=m_volume;
}

*/

bool Image::index_1Dto3D(int index_voxel, int& i, int& j, int& k) const
{
  if ((index_voxel<0) || (index_voxel>=NbVoxels)) return false;
  k=index_voxel / (DimInVoxels[0]*DimInVoxels[1]);
  j=index_voxel % (DimInVoxels[0]*DimInVoxels[1]);
  i= j % DimInVoxels[0];
  j = j/DimInVoxels[0];
  return true;
}

int Image::index_3Dto1D(int i,int j, int k) const
{
  if ((i<0)||(i>=DimInVoxels[0])||(j<0)||(j>=DimInVoxels[1])||(k<0)||(k>=DimInVoxels[2])) return -1;
  return (i+j*DimInVoxels[0]+k*DimInVoxels[0]*DimInVoxels[1]);
}

bool Image::coord2index_3D(double x, double y, double z, int& i, int& j, int& k) const
{
  bool inside=((x>=Corner[0]) && (x<=Corner[0]+DimInCm[0]) && (y>=Corner[1]) && (y<=Corner[1]+DimInCm[1]) && (z>=Corner[2]) && (z<=Corner[2]+DimInCm[2]));
  i=int(std::min(int(floor((x-Corner[0])/VoxelSize[0])), DimInVoxels[0]));
  j=int(std::min(int(floor((y-Corner[1])/VoxelSize[1])), DimInVoxels[1]));
  k=int(std::min(int(floor((z-Corner[2])/VoxelSize[2])), DimInVoxels[2]));
  // std::cout<< i<< ' ' << j<< ' ' << k << ' ' << '\n';
  return(inside);
}

int Image::coord2index_1D(double x, double y, double z) const
{
  int i,j,k;
  bool inside=Image::coord2index_3D(x,y,z,i,j,k);
  if (inside) 
    return Image::index_3Dto1D(i,j,k);
  else
    return -1;
}
	
void  Image::initialize(double value) 
{
  Value.assign( NbVoxels, value);
}

void Image::multiply(const std::vector<double> mult, const std::vector<double> divisor)
{
  for(int i=0; i<NbVoxels; i++)
    Value[i] *= mult[i]/divisor[i];

}

void Image::setIntensity(int index_voxel, double value)
{
  Value[index_voxel] = value;
}

double Image::Maximum() const
{
  double max = 0.0;
  for (int i = 0; i<NbVoxels;i++)
    if (Value[i] > max) max = Value[i];
  return max;
}

int Image::readFile(char *image_file)
{
  std::ifstream fp;
  fp.open(image_file, std::ios::in|std::ios::binary);
  if (!fp.is_open())
    {
      std::cout<<"Unknown file "<<image_file<< '\n';
      return 1;
    }
  fp.read((char *)(&(Value[0])), NbVoxels*sizeof(double));
  if (!fp) 
    {
      std::cout << "Error in " << image_file << ": only "<< fp.gcount()<< " items could be read.\n" ;
      return 1;
    }
  fp.close();
  return 0;
}

const int Image::writeFile(char *image_file)
{
  std::ofstream fp;
  fp.open(image_file, std::ios::out|std::ios::binary);
  if (!fp.is_open())
    {
      std::cout<<"Unknown file "<<image_file<< '\n';
      return 1;
    }
  fp.write((char*)(&(Value[0])), NbVoxels*sizeof(double));
  fp.close();
  return 0;
}


bool Image::readFromFile(std::ifstream& fp)
{
  fp.read((char *)(&(Value[0])), NbVoxels*sizeof(double));
  return fp.eof();
  /*
  if (!fp) 
    {
      std::cout << "Error in " << image_file << ": only "<< fp.gcount()<< " items could be read.\n" ;
      return 1;
    }
  */
}

const void Image::writeToFile(std::ofstream& fp)
{
  fp.write((char*)(&(Value[0])), NbVoxels*sizeof(double));
}
