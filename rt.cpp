#include <bits/stdc++.h>
#include "Includes/header.h"

#define INF 1000000000

using namespace std;

Size size;
Ortho ortho;
vector<Light> lights;
vector<Object> objects;
double ambient;
bool superSampling;
int depth;
T3 background;
double w, h;

T3 getDirection(int i, int j) {
    T3 dir;
    dir.x = (ortho.x0 + w/2) + w*j;
    dir.y = (ortho.y1 - h/2) + h*i;
    dir.z = 0;

    return dir;
}

T3 normalize(T3 vec) {
	double normal;
	normal = sqrt((vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z));
	vec.x /= normal;
	vec.y /= normal;
	vec.z /= normal;

	return vec;
}

double intersect(Ray ray, Quad *obj) {
   double  a, b, c, d, e;// Coefficents of equation of..
   double  f, g, h, j, k;// ..quadric surface
   double  acoef, bcoef, ccoef;// Intersection coefficents
   double  dx, dy, dz;// Direction - origin coordinates
   double  disc;// Distance to intersection
   double  root;// Root of distance to intersection
   double  t;// Distance along ray to intersection
   double  x0, y0, z0;// Origin coordinates


   a = obj->a;
   b = obj->b;
   c = obj->c;
   d = obj->d;
   e = obj->e;
   f = obj->f;
   g = obj->g;
   h = obj->h;
   j = obj->j;
   k = obj->k;

   dx = ray.dir.x - ray.org.x;
   dy = ray.dir.y - ray.org.y;
   dz = ray.dir.z - ray.org.z;

   x0 = ray.org.x;
   y0 = ray.org.y;
   z0 = ray.org.z;

   acoef = 2 * f * dx * dz + 2 * e * dy * dz + c * dz * dz + b * dy * dy +
           a * dx * dx + 2 * d * dx * dy;

   bcoef = 2 * b * y0 * dy + 2 * a * x0 * dx + 2 * c * z0 * dz +
           2 * h * dy + 2 * g * dx + 2 * j * dz + 2 * d * x0 * dy +
           2 * e * y0 * dz + 2 * e * z0 * dy + 2 * d * y0 * dx +
           2 * f * x0 * dz + 2 * f * z0 * dx;

   ccoef = a * x0 * x0 + 2 * g * x0 + 2 * f * x0 * z0 + b * y0 * y0 +
           2 * e * y0 * z0 + 2 * d * x0 * y0 + c * z0 * z0 + 2 * h * y0 +
           2 * j * z0 + k;


   if ( 1.0 + acoef == 1.0 ) {
      if ( 1.0 + bcoef == 1.0 ) {
         return -1.0;
      }

      t = ( -ccoef ) / ( bcoef );

   } else {
      disc = bcoef * bcoef - 4 * acoef * ccoef;
      if ( 1.0 + disc < 1.0 ) {
         return -1.0;
      }

      root = sqrt( disc );
      t = ( -bcoef - root ) / ( acoef + acoef );
      if ( t < 0.0 ) {
         t = ( -bcoef + root ) / ( acoef + acoef );
      }
   }

   if ( t < 0.001 )
      return -1.0;

   return t;
}


int objIndex() {
	int ret = -1;
	double dist = INF, aux = -1.0;

	for(int i = 0; i < objects.size(); i++) {
		aux = 

		if(aux > 0.0 && aux < dist) {
			dist = aux;
			ret = i;
		}
	}
	if(dist != INF) return ret;
	return -1;
}

int main() {
	SDL* sdl = new SDL("onesphere.sdl");
	cout << sdl->getOutput() << endl;

	size = sdl.getSize();
	ortho = sdl.getOrtho();
	lights = sdl.getLights();
	objects = sdl.getObjects();
	background = sdl.getBackground();
	ambient = sdl.getAmbient();
	superSampling = sdl.getSuperSampling();
	depth = sdl.getDepth();
	w = fabs(ortho.x1 - ortho.x0)/size.w;
	h = fabs(orthor.y1 - ortho.y0)/size.h;


	for(int i = 0; i < size.h; i++) {
		for(int j = 0; i < size.w; i++) {
			T3 dir = (getDiretion(i, j) - sdl.getEye());
			dir = normalize(dir);
		}
	}


	delete sdl;
	return 0;
}