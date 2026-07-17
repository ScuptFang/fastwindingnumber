#pragma once

#include <cmath>
#include <ostream>

namespace fwn {

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    double operator[](int axis) const {
        return axis == 0 ? x : (axis == 1 ? y : z);
    }

    double& operator[](int axis) {
        return axis == 0 ? x : (axis == 1 ? y : z);
    }
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vec3 operator-(const Vec3& v) {
    return Vec3(-v.x, -v.y, -v.z);
}

inline Vec3 operator*(const Vec3& v, double s) {
    return Vec3(v.x * s, v.y * s, v.z * s);
}

inline Vec3 operator*(double s, const Vec3& v) {
    return v * s;
}

inline Vec3 operator/(const Vec3& v, double s) {
    return Vec3(v.x / s, v.y / s, v.z / s);
}

inline Vec3& operator+=(Vec3& a, const Vec3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline Vec3& operator-=(Vec3& a, const Vec3& b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

inline Vec3& operator*=(Vec3& v, double s) {
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

inline double squaredNorm(const Vec3& v) {
    return dot(v, v);
}

inline double norm(const Vec3& v) {
    return std::sqrt(squaredNorm(v));
}

inline Vec3 normalized(const Vec3& v) {
    const double n = norm(v);
    return n > 0.0 ? v / n : Vec3();
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    os << v.x << ' ' << v.y << ' ' << v.z;
    return os;
}

}  // namespace fwn
