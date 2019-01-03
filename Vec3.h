//
// Copyright (C) 2019 by Yuri Victorovich. All rights reserved.
//

#pragma once

#include <boost/array.hpp>
#include <ostream>
#include <math.h>

typedef double Float;

/// Vector helper classes
enum VecCoords {X=1, Y=2, Z=3};

class Vec3 : public boost::array<Float, 3> {
public:
  typedef Float val_t;
  Vec3() { }
  Vec3(Float x, Float y, Float z) : boost::array<Float, 3>({{x, y, z}}) { }
  Float& operator()(unsigned idx) { // 1-based index access
    return (*this)[idx-1];
  }
  Float operator()(unsigned idx) const { // 1-based index access
    return (*this)[idx-1];
  }
  static Vec3 one(unsigned idx, Float val) {Vec3 vec(0,0,0); vec(idx) = val; return vec;}
  Float len2() const {return sq((*this)(X)) + sq((*this)(Y)) + sq((*this)(Z));}
  Float len() const {return sqrt(len2());}
  Vec3 normalizeZ() const {
    auto l = len();
    if (l != 0)
      return (*this)/l;
    else
      return Vec3(0,0,0);
  }
  friend std::ostream& operator<<(std::ostream &os, const Vec3 &v) {
    os << "{" << v(X) << "," << v(Y) << "," << v(Z) << "}";
    return os;
  }
  Vec3 operator-() const {
    const Vec3 &i = *this;
    return Vec3(-i(X), -i(Y), -i(Z));
  }
  Vec3 operator*(Float m) const {
    const Vec3 &i = *this;
    return Vec3(i(X)*m, i(Y)*m, i(Z)*m);
  }
  Vec3& operator*=(Float m) {
    Vec3 &i = *this;
    i(X) *= m;
    i(Y) *= m;
    i(Z) *= m;
    return *this;
  }
  Float operator*(const Vec3 &v) const { // scalar multiplication
    const Vec3 &i = *this;
    return i(X)*v(X) + i(Y)*v(Y) + i(Z)*v(Z);
  }
  Vec3 cross(const Vec3 &v) const { // cross (vector) multiplication
    const Vec3 &i = *this;
    return Vec3(-i(Z)*v(Y) + i(Y)*v(Z), i(Z)*v(X) - i(X)*v(Z), -i(Y)*v(X) + i(X)*v(Y)); // see https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
  }
  Vec3 operator/(Float m) const {
    const Vec3 &i = *this;
    return Vec3(i(X)/m, i(Y)/m, i(Z)/m);
  }
  Vec3 operator+(const Vec3 &v) const {
    const Vec3 &i = *this;
    return Vec3(i(X)+v(X), i(Y)+v(Y), i(Z)+v(Z));
  }
  Vec3& operator+=(const Vec3 &v) {
    Vec3 &i = *this;
    i(X) += v(X);
    i(Y) += v(Y);
    i(Z) += v(Z);
    return *this;
  }
  Vec3 operator-(const Vec3 &v) const {
    const Vec3 &i = *this;
    return Vec3(i(X)-v(X), i(Y)-v(Y), i(Z)-v(Z));
  }
  Vec3& operator-=(const Vec3 &v) {
    Vec3 &i = *this;
    i(X) -= v(X);
    i(Y) -= v(Y);
    i(Z) -= v(Z);
    return *this;
  }
  Vec3 project(const Vec3 &dir) const { // dir is assimed to be a normal vector or zero vector
    return dir*((*this)*dir);
  }
  Vec3 divOneByOne(const Vec3 &d) const { // dir is assimed to be a normal vector or zero vector
    const Vec3 &i = *this;
    return Vec3(i(X)/d(X), i(Y)/d(Y), i(Z)/d(Z));
  }
private:
  static Float sq(Float f) {return f*f;}
};
