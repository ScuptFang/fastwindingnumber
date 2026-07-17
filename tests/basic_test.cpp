#include "fwn/FastWindingNumber.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void addOrientedFace(fwn::Mesh& mesh, int a, int b, int c, const fwn::Vec3& inside_hint) {
    const fwn::Vec3& va = mesh.vertices[static_cast<std::size_t>(a)];
    const fwn::Vec3& vb = mesh.vertices[static_cast<std::size_t>(b)];
    const fwn::Vec3& vc = mesh.vertices[static_cast<std::size_t>(c)];
    const fwn::Vec3 centroid = (va + vb + vc) / 3.0;
    const fwn::Vec3 normal = fwn::cross(vb - va, vc - va);

    if (fwn::dot(normal, inside_hint - centroid) > 0.0) {
        mesh.triangles.emplace_back(a, c, b);
    } else {
        mesh.triangles.emplace_back(a, b, c);
    }
}

fwn::Mesh makeTetrahedron() {
    fwn::Mesh mesh;
    mesh.vertices = {
        fwn::Vec3(1.0, 1.0, 1.0),
        fwn::Vec3(-1.0, -1.0, 1.0),
        fwn::Vec3(-1.0, 1.0, -1.0),
        fwn::Vec3(1.0, -1.0, -1.0)
    };

    const fwn::Vec3 inside_hint(0.0, 0.0, 0.0);
    addOrientedFace(mesh, 0, 2, 1, inside_hint);
    addOrientedFace(mesh, 0, 1, 3, inside_hint);
    addOrientedFace(mesh, 0, 3, 2, inside_hint);
    addOrientedFace(mesh, 1, 2, 3, inside_hint);
    return mesh;
}

bool near(double a, double b, double eps) {
    return std::abs(a - b) <= eps;
}

fwn::Vec3 torusPoint(
    double major_radius,
    double minor_radius,
    int u_index,
    int v_index,
    int u_segments,
    int v_segments
) {
    const double u = 2.0 * kPi * static_cast<double>(u_index) / static_cast<double>(u_segments);
    const double v = 2.0 * kPi * static_cast<double>(v_index) / static_cast<double>(v_segments);
    const double ring_radius = major_radius + minor_radius * std::cos(v);
    return fwn::Vec3(
        ring_radius * std::cos(u),
        ring_radius * std::sin(u),
        minor_radius * std::sin(v)
    );
}

fwn::Mesh makeTorus(double major_radius, double minor_radius, int u_segments, int v_segments) {
    fwn::Mesh mesh;
    mesh.vertices.reserve(static_cast<std::size_t>(u_segments * v_segments));
    mesh.triangles.reserve(static_cast<std::size_t>(u_segments * v_segments * 2));

    for (int u = 0; u < u_segments; ++u) {
        for (int v = 0; v < v_segments; ++v) {
            mesh.vertices.push_back(torusPoint(major_radius, minor_radius, u, v, u_segments, v_segments));
        }
    }

    const auto index = [v_segments](int u, int v) {
        return u * v_segments + v;
    };

    for (int u = 0; u < u_segments; ++u) {
        const int next_u = (u + 1) % u_segments;
        for (int v = 0; v < v_segments; ++v) {
            const int next_v = (v + 1) % v_segments;
            const int p00 = index(u, v);
            const int p10 = index(next_u, v);
            const int p01 = index(u, next_v);
            const int p11 = index(next_u, next_v);
            mesh.triangles.emplace_back(p00, p10, p11);
            mesh.triangles.emplace_back(p00, p11, p01);
        }
    }
    return mesh;
}

}  // namespace

int main() {
    const fwn::Mesh mesh = makeTetrahedron();
    const fwn::FastWindingNumber fwn(mesh, 2);

    const double inside = fwn.exactWindingNumber(fwn::Vec3(0.0, 0.0, 0.0));
    const double outside = fwn.exactWindingNumber(fwn::Vec3(8.0, 0.0, 0.0));

    if (!near(std::abs(inside), 1.0, 1e-12)) {
        std::cerr << "inside exact winding number expected abs 1, got " << inside << '\n';
        return 1;
    }
    if (!near(outside, 0.0, 1e-12)) {
        std::cerr << "outside exact winding number expected 0, got " << outside << '\n';
        return 1;
    }

    fwn::QueryOptions exact;
    exact.fast = false;
    const double via_options = fwn.windingNumber(fwn::Vec3(0.0, 0.0, 0.0), exact);
    if (!near(via_options, inside, 1e-12)) {
        std::cerr << "exact option mismatch: " << via_options << " vs " << inside << '\n';
        return 1;
    }

    fwn::QueryOptions fast;
    fast.fast = true;
    fast.beta = 2.0;
    const double far_fast = fwn.windingNumber(fwn::Vec3(12.0, 3.0, -2.0), fast);
    const double far_exact = fwn.exactWindingNumber(fwn::Vec3(12.0, 3.0, -2.0));
    if (std::abs(far_fast - far_exact) > 1e-3) {
        std::cerr << "far-field approximation drift too large: "
                  << far_fast << " vs " << far_exact << '\n';
        return 1;
    }

    const double major_radius = 1.0;
    const double minor_radius = 0.35;
    const fwn::Mesh torus = makeTorus(major_radius, minor_radius, 64, 24);
    const fwn::FastWindingNumber torus_fwn(torus, 8);
    fwn::QueryOptions accurate_fast;
    accurate_fast.fast = true;
    accurate_fast.beta = 5.0;
    accurate_fast.expansion_order = 2;

    std::mt19937 rng(1337);
    std::uniform_real_distribution<double> xy_dist(-(major_radius + minor_radius), major_radius + minor_radius);
    std::uniform_real_distribution<double> z_dist(-minor_radius, minor_radius);

    double worst_error = 0.0;
    double total_error = 0.0;
    const int query_count = 128;
    for (int i = 0; i < query_count; ++i) {
        const fwn::Vec3 q(xy_dist(rng), xy_dist(rng), z_dist(rng));
        const double fast_value = torus_fwn.windingNumber(q, accurate_fast);
        const double exact_value = torus_fwn.exactWindingNumber(q);
        const double error = std::abs(fast_value - exact_value);
        total_error += error;
        worst_error = std::max(worst_error, error);
    }

    const double mean_error = total_error / static_cast<double>(query_count);
    if (worst_error > 1e-3 || mean_error > 1e-4) {
        std::cerr << "torus random-query approximation error too large: worst "
                  << worst_error << ", mean " << mean_error << '\n';
        return 1;
    }

    return 0;
}
