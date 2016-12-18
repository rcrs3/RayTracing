#ifndef RT_HEADER_H
#define RT_HEADER_H

#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>

using namespace std;

typedef struct T3 {
    double x, y, z;
    T3() {}
    T3(double x, double y, double z) : x(x), y(y), z(z) {}
    
     T3 operator +(const T3 &b) const {
        T3 ret;
        ret.x = x + b.x;
        ret.y = y + b.y;
        ret.z = z + b.z;
        return ret;
    }
    T3 operator -(const T3 &b) const {
        T3 ret;
        ret.x = x - b.x;
        ret.y = y - b.y;
        ret.z = z - b.z;
        return ret;
    }
    double operator &(const T3 &b) const {
        T3 ret;
        ret.x = x * b.x;
        ret.y = y * b.y;
        ret.z = z * b.z;
        return ret.x + ret.y + ret.z;
    }

    T3 operator *(const double &b) const {
        T3 ret;
        ret.x = x * b;
        ret.y = y * b;
        ret.z = z * b;
        return ret;
    }
    T3 operator /(const double &b) const {
        T3 ret;
        ret.x = x / b;
        ret.y = y / b;
        ret.z = z / b;
        return ret;
    }

} T3;

typedef struct Color {
    double r, g, b;
    Color() {}
} Color;

typedef struct Ortho {
    double x0, y0, x1, y1;
} Ortho;

typedef struct Size {
    int w, h;
} Size;

typedef struct Light {
    T3 coords;
    double intensity;

    Light(){
        this->coords = T3();
    }
} Light;

typedef struct Object {
    double a, b, c, d, e ,f, g, h, j, k, ka, kd, ks, n, KS, KT, ir;
    Color color;

    Object(){
        this->color = Color();
    }
} Object;

typedef struct Obj {
    std::string name;
    double ka, kd,ks, n, KS, KT, ir;
    Color color;

    Obj() {
        this->color = Color();
    }
} Obj;

typedef struct Ray {
    T3  org;// Origin of ray
   T3  dir;// Direction of ray
   int  depth;// Depth (or length) of ray
   Ray() {}
   Ray(T3 org, T3 dir, int depth) : org(org), dir(dir), depth(depth) {}
} Ray;

vector<T3> V[3];
vector<T3> F[3];

void readObj(string file) {
    ifstream ifs;
    ifs.open (file, ifstream::in);
    int cont = 0;
    while(!ifs.eof()) {
        if(ifs.peek() == '\n' || ifs.peek() == ' ') {
            ifs.get();
        } else if(ifs.peek() == '#'){
            ifs.ignore(numeric_limits<streamsize>::max(), '\n');
        } else {
            string tag;
            ifs >> tag;
            if(!tag.compare("v")){
                T3 aux = T3();
                ifs >> aux.x
                    >> aux.y
                    >> aux.z;
                V[0].push_back(aux);
            } else if(!tag.compare("f")) {
                T3 aux[3];
                string s[3];
                
                ifs >> s[0]
                    >> s[1]
                    >> s[2];
                
                if(s[0].find("//") != -1) {
                   int ind = s[0].find_first_of("/");
                   aux[0].x = stoi(s[0].substr(0, ind-1));
                   aux[2].x = stoi(s[0].substr(ind+2, s[0].size()-1));

                   ind = s[1].find_first_of("/");
                   aux[0].y = stoi(s[1].substr(0, ind-1));
                   aux[2].y = stoi(s[1].substr(ind+2, s[1].size()-1));

                   ind = s[2].find_first_of("/");
                   aux[0].z = stoi(s[2].substr(0, ind-1));
                   aux[2].z = stoi(s[2].substr(ind+2, s[2].size()-1));

                    F[0].push_back(aux[0]);
                    F[2].push_back(aux[2]);

                } else if(s[0].find("/") != -1) {
                    int ind = s[0].find_first_of("/");
                    int ind2 = s[0].find_last_of("/");
                    
                    aux[0].x = stoi(s[0].substr(0, ind-1));
                    aux[1].x = stoi(s[0].substr(ind+1, ind2-1));
                    aux[2].x = stoi(s[0].substr(ind2+1, s[0].size()-1));

                    ind = s[1].find_first_of("/");
                    ind2 = s[1].find_last_of("/");
                    
                    aux[0].y = stoi(s[1].substr(0, ind-1));
                    aux[1].y = stoi(s[1].substr(ind+1, ind2-1));
                    aux[2].y = stoi(s[1].substr(ind2+1, s[1].size()-1));

                    ind = s[2].find_first_of("/");
                    ind2 = s[2].find_last_of("/");
                    
                    aux[0].z = stoi(s[2].substr(0, ind-1));
                    aux[1].z = stoi(s[2].substr(ind+1, ind2-1));
                    aux[2].z = stoi(s[2].substr(ind2+1, s[2].size()-1));

                    F[0].push_back(aux[0]);
                    F[1].push_back(aux[1]);
                    F[2].push_back(aux[2]);

                } else {
                    aux[0].x = stoi(s[0]);
                    aux[0].y = stoi(s[1]);
                    aux[0].z = stoi(s[2]);

                    F[0].push_back(aux[0]);
                }

            } else if(!tag.compare("vt")) {
                T3 aux = T3();
                ifs >> aux.x
                    >> aux.y
                    >> aux.z;
                V[1].push_back(aux);

            } else if(!tag.compare("vn")) {
                T3 aux = T3();
                ifs >> aux.x
                    >> aux.y
                    >> aux.z;
                V[2].push_back(aux);

            }

        }

    }
}

#include "SDL.cpp"
#endif
