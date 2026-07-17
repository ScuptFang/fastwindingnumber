#pragma once

#include <array>
#include <vector>

#include "fwn/Vec3.h"

namespace fwn {

struct Triangle {
    int v0;
    int v1;
    int v2;

    Triangle() : v0(0), v1(0), v2(0) {}
    Triangle(int a, int b, int c) : v0(a), v1(b), v2(c) {}
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Triangle> triangles;

    bool empty() const {
        return vertices.empty() || triangles.empty();
    }
};

}  // namespace fwn
