#ifndef RT_HEADER_H
#define RT_HEADER_H

#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#ifndef _GLIBCXX_NO_ASSERT
#include <cassert>
#endif
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#if __cplusplus >= 201103L
#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cwchar>
#include <cwctype>
#endif

// C++
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

#if __cplusplus >= 201103L
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <forward_list>
#include <future>
#include <initializer_list>
#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <scoped_allocator>
#include <system_error>
#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#endif

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
    //string name;
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

#include "SDL.cpp"
#endif
