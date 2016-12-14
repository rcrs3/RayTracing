#include <math.h>

struct PointVector{
  double x,y,z;
  PointVector(){};
  PointVector(double x, double y, double z):x(x),y(y),z(z){};
  PointVector operator+(PointVector p){
    return PointVector(x+p.x, y+p.y, z+p.z);
  }
  PointVector operator-(PointVector p){
    return PointVector(x-p.x, y-p.y, z-p.z);
  }
  PointVector operator*(PointVector p){
    return PointVector(x*p.x, y*p.y, z*p.z);
  }
  PointVector operator*(double p){
    return PointVector(x*p, y*p, z*p);
  }
};

struct Triangle{
    PointVector x,y,z;
    Triangle(PointVector a, PointVector b, PointVector c): x(a), y(b), z(c){};
};
