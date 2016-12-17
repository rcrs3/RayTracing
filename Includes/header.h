#ifndef RT_HEADER_H
#define RT_HEADER_H

#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

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
    T3 operator *(const T3 &b) const {
        T3 ret;
        ret.x = x * b.x;
        ret.y = y * b.y;
        ret.z = z * b.z;
        return ret;
    }

    T3 operator *(const double &b) const {
        T3 ret;
        ret.x = x * b;
        ret.y = y * b;
        ret.z = z * b;
        return ret;
    }

} T3;

typedef struct Color {
    double r, g, b;
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

    light(){
        this->coords = new T3();
    }
    ~light(){
        delete this->coords;
    }
} Light;

typedef struct Object {
    double a, b, c, d, e ,f, g, h, j, k, ka, kd, ks, n, KS, KT, ir;
    Color* color;

    object(){
        this->color = new Color();
    }
    ~object(){
        delete this->color;
    }
} Object;

typedef struct Ray {
    T3  org;// Origin of ray
   T3  dir;// Direction of ray
   int  depth;// Depth (or length) of ray
   Ray() {}
   Ray(T3 org, T3 dir, int depth) : org(org), dir(dir), depth(depth) {}
} Ray;

#include "SDL.cpp"
#endif