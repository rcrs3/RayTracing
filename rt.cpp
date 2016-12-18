#include <bits/stdc++.h>
#include "Includes/header.h"

#define INF 1000000000
#define EPS 1e-12
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
    dir.y = (ortho.y1 - h/2) - h*i;
    dir.z = 0;

    return dir;
}

T3 normalize(T3 v) {
   double denom;// Temporary denominator
   //  Absolute value of vector's coordinates

   double x = ( v.x > 0.0 ) ? v.x : - v.x;
   double y = ( v.y > 0.0 ) ? v.y : - v.y;
   double z = ( v.z > 0.0 ) ? v.z : - v.z;


   if ( x > y ) {
      if ( x > z ) {
         y = y / x;
         z = z / x;
         denom = 1.0 / ( x * sqrt( 1.0 + y * y + z * z ) );

      } else { // z > x > y
         if ( 1.0 + z > 1.0 ) {
            y = y / z;
            x = x / z;
            denom = 1.0 / ( z * sqrt( 1.0 + y * y + x * x ) );
         }
      }

   } else {
      if ( y > z ) {
         z = z / y;
         x = x / y;
         denom = 1.0 / ( y * sqrt( 1.0 + z * z + x * x ) );

      } else { // x < y < z
         if ( 1.0 + z > 1.0 ) {
            y = y / z;
            x = x / z;
            denom = 1.0 / ( z * sqrt( 1.0 + y * y + x * x ) );
         }
      }
   }

   if ( 1.0 + x + y + z > 1.0 ) {
      v = v * denom;
   }
   return v;
}// End procedure normalize

double intersect(Ray ray, Object *obj) {
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

   ray.dir = normalize(ray.dir);

   dx = ray.dir.x;
   dy = ray.dir.y;
   dz = ray.dir.z;

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


int objIndex(Ray ray) {
	int ret = -1;
	double dist = INF, aux = -1.0;

	for(int i = 0; i < objects.size(); i++) {
		aux = intersect(ray, &objects[i]);

		if(aux > 0.0 && aux < dist) {
			dist = aux;
			ret = i;
		}
	}
	return ret;
}

T3 intersectPoint(Ray ray, Object obj) {
    T3 ret = ray.dir;
    ret = normalize(ret);
    double dist = intersect(ray, &obj);
    ret = ret * dist;
    ret = ret + ray.org;
    return ret;
}

T3 normalQuadric(Object obj, T3 p) {
    double x = p.x;
    double y = p.y;
    double z = p.z;
    double a = obj.a;
    double b = obj.b;
    double c = obj.c;
    double d = obj.d;
    double e = obj.e;
    double f = obj.f;
    double g = obj.g;
    double h = obj.h;
    double j = obj.j;
    p.x = (2*a*x) + (2*d*y) + (2*f*z) + (2*g);
    p.y = (2*b*y) + (2*d*x) + (2*e*z) + (2*h);
    p.z = (2*c*z) + (2*e*y) + (2*f*x) + (2*j);

    return p;
}

bool isShadow(Ray ray) {
	for(int i = 0; i < objects.size(); i++) {
		if(intersect(ray, &objects[i]) != -1.0) return true;
	}
	return false;
}

double intensityS(Ray ray, Object obj, Light light) {
	T3 iP = intersectPoint(ray, obj);
    T3 normalObj = normalQuadric(obj, iP);
    normalObj = normalize(normalObj);

    T3 VLight = light.coords;
    VLight = normalize(VLight);

    double c = VLight & normalObj;
    c = max(c, 0.0);
    T3 r = (normalObj*(2*c)) - VLight;
    r = normalize(r);

    Ray RLight = Ray(iP, VLight, 0);
    if(!isShadow(RLight)){
        T3 v = ray.dir;
        v = v * -1;
        v = normalize(v);
        double rv = r & v;
        rv = max(rv, 0.0);

        return light.intensity * obj.ks * pow(rv, obj.n);
    }
    return 0.0;
}

double intensityD(Ray ray, Object obj, Light light) {
	T3 iP = intersectPoint(ray, obj);
	T3 normalObj = normalQuadric(obj, iP);
	normalObj = normalize(normalObj);

	T3 VLight = light.coords;
	VLight = normalize(VLight);

	Ray RLight = Ray(iP, VLight, 0);
	if(!isShadow(RLight)) {
		double nl = VLight & normalObj;
		nl = max(nl, 0.0);
		double ret = light.intensity * obj.kd * nl;
		return ret;
	}
	return 0.0; 
}

void normalizeIlumination(T3 *v) {
	v->x = min(v->x, 1.0);
	v->y = min(v->y, 1.0);
	v->z = min(v->z, 1.0);
}

T3 shadow(Ray ray, Object obj) {
	T3 ret;
	double Ia, Id, Is;
	Ia = obj.ka*ambient;
	Id = Is = 0.0;
	
	for(int i = 0; i < lights.size(); i++) {
		Id += intensityD(ray, obj, lights[i]);
		Is += intensityS(ray, obj, lights[i]);
	}

	double aux = Ia + Id;
	ret.x = obj.color.r*aux + Is;
	ret.y = obj.color.g*aux + Is;
	ret.z = obj.color.b*aux + Is;

	normalizeIlumination(&ret);

	return ret;
}

Ray reflection(Ray ray, Object obj) {
	T3 iP = intersectPoint(ray, obj);
	T3 normalObj = normalQuadric(obj, iP);
	normalObj = normalize(normalObj);

	T3 RLight = ray.dir * -1.0;
	RLight = normalize(RLight);

	double nl = RLight & normalObj;

	return Ray(iP, (normalObj * (2.0 * nl)) - RLight, ray.depth-1);
}

T3 calcColor(Ray ray) {
	T3 color;
	int ind = objIndex(ray);
	if(ind < 0) {
		if(ray.depth == depth) {
			color = background;
			return color;
		}
		color.x = color.y = color.z = 0.0;
		return color;
	}
	color = shadow(ray, objects[ind]);

	T3 refC = T3(0.0, 0.0, 0.0);

	if(ray.depth > 0) {
		if(objects[ind].KS > 0.0) {
			Ray ref = reflection(ray, objects[ind]);
			refC = calcColor(ref);
		}
	}

	color.x = ((1 - objects[ind].KS - objects[ind].KT) * color.x) + (objects[ind].KS * refC.x);
	color.y = ((1 - objects[ind].KS - objects[ind].KT) * color.y) + (objects[ind].KS * refC.y);
	color.z = ((1 - objects[ind].KS - objects[ind].KT) * color.z) + (objects[ind].KS * refC.z);

	normalizeIlumination(&color);

	return color;
}

int main() {
	SDL sdl = SDL("twoplanesphere.sdl");

	size = sdl.getSize();
	ortho = sdl.getOrtho();
	lights = sdl.getLights();
	objects = sdl.getObjects();
	background = sdl.getBackground();
	ambient = sdl.getAmbient();
	superSampling = sdl.getSuperSampling();
	depth = sdl.getDepth();
	w = fabs(ortho.x1 - ortho.x0)/size.w;
	h = fabs(ortho.y1 - ortho.y0)/size.h;

	ofstream ofs(sdl.getOutput(), ios::out | ios::binary); 
    ofs << "P6\n" << size.w << " " << size.h << "\n255\n";

	for(int i = 0; i < size.h; i++) {
		for(int j = 0; j < size.w; j++) {
			T3 color = T3(0.0, 0.0, 0.0);
			T3 dir = T3();
			Ray ray = Ray();

			if(superSampling) {
				for(double subPi = i; subPi <= i + 1.0; subPi += 0.5) {
					for(double subPj = j; subPj <= j + 1.0; subPj += 0.5) {
						dir = getDirection(subPi, subPj) - sdl.getEye();
						ray = Ray(sdl.getEye(), dir, depth);
						color = color + calcColor(ray);
					}
				}
			
				color = color/9.0;
			} else {
				dir = getDirection(i, j) - sdl.getEye();
				ray = Ray(sdl.getEye(), dir, depth);
				color = calcColor(ray);
			}
			
			ofs << (unsigned char)(color.x * 255)
                << (unsigned char)(color.y * 255)
                << (unsigned char)(color.z * 255);
		}
	}
	ofs.close();
	return 0;
}
