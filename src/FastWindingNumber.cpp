#include "fwn/FastWindingNumber.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace fwn {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kFourPi = 4.0 * kPi;

std::size_t moment1Index(int area_axis, int offset_axis) {
    return static_cast<std::size_t>(area_axis * 3 + offset_axis);
}

std::size_t moment2Index(int area_axis, int offset_axis_a, int offset_axis_b) {
    return static_cast<std::size_t>(area_axis * 9 + offset_axis_a * 3 + offset_axis_b);
}

double delta(int a, int b) {
    return a == b ? 1.0 : 0.0;
}

}  // namespace

void FastWindingNumber::Aabb::expand(const Vec3& p) {
    if (!valid) {
        min = p;
        max = p;
        valid = true;
        return;
    }
    min.x = std::min(min.x, p.x);
    min.y = std::min(min.y, p.y);
    min.z = std::min(min.z, p.z);
    max.x = std::max(max.x, p.x);
    max.y = std::max(max.y, p.y);
    max.z = std::max(max.z, p.z);
}

void FastWindingNumber::Aabb::expand(const Aabb& box) {
    if (!box.valid) {
        return;
    }
    expand(box.min);
    expand(box.max);
}

Vec3 FastWindingNumber::Aabb::center() const {
    return (min + max) * 0.5;
}

Vec3 FastWindingNumber::Aabb::extent() const {
    return max - min;
}

FastWindingNumber::FastWindingNumber(const Mesh& mesh, std::size_t leaf_size)
    : leaf_size_(std::max<std::size_t>(1, leaf_size)) {
    if (mesh.vertices.empty()) {
        throw std::invalid_argument("mesh has no vertices");
    }
    if (mesh.triangles.empty()) {
        throw std::invalid_argument("mesh has no triangles");
    }

    triangles_.reserve(mesh.triangles.size());
    for (const Triangle& t : mesh.triangles) {
        if (t.v0 < 0 || t.v1 < 0 || t.v2 < 0 ||
            static_cast<std::size_t>(t.v0) >= mesh.vertices.size() ||
            static_cast<std::size_t>(t.v1) >= mesh.vertices.size() ||
            static_cast<std::size_t>(t.v2) >= mesh.vertices.size()) {
            throw std::invalid_argument("triangle index is outside the vertex array");
        }

        TriangleData tri;
        tri.v0 = mesh.vertices[static_cast<std::size_t>(t.v0)];
        tri.v1 = mesh.vertices[static_cast<std::size_t>(t.v1)];
        tri.v2 = mesh.vertices[static_cast<std::size_t>(t.v2)];
        tri.centroid = (tri.v0 + tri.v1 + tri.v2) / 3.0;
        tri.area_vector = cross(tri.v1 - tri.v0, tri.v2 - tri.v0) * 0.5;
        tri.area = norm(tri.area_vector);
        tri.bounds.expand(tri.v0);
        tri.bounds.expand(tri.v1);
        tri.bounds.expand(tri.v2);
        triangles_.push_back(tri);
    }

    order_.resize(triangles_.size());
    std::iota(order_.begin(), order_.end(), 0);
    nodes_.reserve(triangles_.size() * 2);
    root_ = buildNode(0, order_.size());
}

double FastWindingNumber::windingNumber(const Vec3& query, const QueryOptions& options) const {
    if (nodes_.empty()) {
        return 0.0;
    }
    const double solid_angle = options.fast
        ? solidAngleRecursive(root_, query, options)
        : exactWindingNumber(query, options.singular_epsilon) * kFourPi;
    return solid_angle / kFourPi;
}

double FastWindingNumber::exactWindingNumber(const Vec3& query, double singular_epsilon) const {
    double solid_angle = 0.0;
    for (const TriangleData& tri : triangles_) {
        solid_angle += triangleSolidAngle(tri, query, singular_epsilon);
    }
    return solid_angle / kFourPi;
}

std::vector<double> FastWindingNumber::windingNumbers(
    const std::vector<Vec3>& queries,
    const QueryOptions& options
) const {
    std::vector<double> values;
    values.reserve(queries.size());
    for (const Vec3& q : queries) {
        values.push_back(windingNumber(q, options));
    }
    return values;
}

std::size_t FastWindingNumber::buildNode(std::size_t begin, std::size_t end) {
    const std::size_t node_index = nodes_.size();
    nodes_.push_back(Node());
    Node& node = nodes_[node_index];
    node.begin = begin;
    node.end = end;

    Vec3 weighted_center;
    for (std::size_t i = begin; i < end; ++i) {
        const TriangleData& tri = triangles_[order_[i]];
        node.bounds.expand(tri.bounds);
        node.area_vector += tri.area_vector;
        node.total_area += tri.area;
        weighted_center += tri.centroid * tri.area;
    }

    node.center = node.total_area > 0.0
        ? weighted_center / node.total_area
        : node.bounds.center();

    for (std::size_t i = begin; i < end; ++i) {
        const TriangleData& tri = triangles_[order_[i]];
        node.radius = std::max(node.radius, norm(tri.v0 - node.center));
        node.radius = std::max(node.radius, norm(tri.v1 - node.center));
        node.radius = std::max(node.radius, norm(tri.v2 - node.center));

        const Vec3 offset = tri.centroid - node.center;
        for (int j = 0; j < 3; ++j) {
            const double area_component = tri.area_vector[j];
            for (int k = 0; k < 3; ++k) {
                node.first_moment[moment1Index(j, k)] += area_component * offset[k];
                for (int l = 0; l < 3; ++l) {
                    node.second_moment[moment2Index(j, k, l)] += area_component * offset[k] * offset[l];
                }
            }
        }
    }

    const std::size_t count = end - begin;
    if (count <= leaf_size_) {
        return node_index;
    }

    Aabb centroid_bounds;
    for (std::size_t i = begin; i < end; ++i) {
        centroid_bounds.expand(triangles_[order_[i]].centroid);
    }
    const Vec3 extent = centroid_bounds.extent();
    int axis = 0;
    if (extent.y > extent.x && extent.y >= extent.z) {
        axis = 1;
    } else if (extent.z > extent.x && extent.z > extent.y) {
        axis = 2;
    }

    const std::size_t mid = begin + count / 2;
    std::nth_element(
        order_.begin() + static_cast<std::ptrdiff_t>(begin),
        order_.begin() + static_cast<std::ptrdiff_t>(mid),
        order_.begin() + static_cast<std::ptrdiff_t>(end),
        [&](std::size_t a, std::size_t b) {
            return triangles_[a].centroid[axis] < triangles_[b].centroid[axis];
        }
    );

    const int left = static_cast<int>(buildNode(begin, mid));
    const int right = static_cast<int>(buildNode(mid, end));
    nodes_[node_index].left = left;
    nodes_[node_index].right = right;
    return node_index;
}

double FastWindingNumber::solidAngleRecursive(
    std::size_t node_index,
    const Vec3& query,
    const QueryOptions& options
) const {
    const Node& node = nodes_[node_index];
    const double dist = norm(node.center - query);
    const double beta = std::max(0.0, options.beta);

    if (!node.isLeaf() && node.radius > 0.0 && dist > beta * node.radius) {
        return multipoleSolidAngle(node, query, options.expansion_order, options.singular_epsilon);
    }

    if (node.isLeaf()) {
        return exactSolidAngleRange(node.begin, node.end, query, options.singular_epsilon);
    }

    double solid_angle = 0.0;
    if (node.left >= 0) {
        solid_angle += solidAngleRecursive(static_cast<std::size_t>(node.left), query, options);
    }
    if (node.right >= 0) {
        solid_angle += solidAngleRecursive(static_cast<std::size_t>(node.right), query, options);
    }
    return solid_angle;
}

double FastWindingNumber::exactSolidAngleRange(
    std::size_t begin,
    std::size_t end,
    const Vec3& query,
    double eps
) const {
    double solid_angle = 0.0;
    for (std::size_t i = begin; i < end; ++i) {
        solid_angle += triangleSolidAngle(triangles_[order_[i]], query, eps);
    }
    return solid_angle;
}

double FastWindingNumber::triangleSolidAngle(const TriangleData& tri, const Vec3& query, double eps) {
    const Vec3 a = tri.v0 - query;
    const Vec3 b = tri.v1 - query;
    const Vec3 c = tri.v2 - query;

    const double la = norm(a);
    const double lb = norm(b);
    const double lc = norm(c);
    if (la <= eps || lb <= eps || lc <= eps) {
        return 0.0;
    }

    const double numerator = dot(a, cross(b, c));
    const double denominator =
        la * lb * lc +
        dot(a, b) * lc +
        dot(b, c) * la +
        dot(c, a) * lb;

    if (std::abs(numerator) <= eps && std::abs(denominator) <= eps) {
        return 0.0;
    }
    return 2.0 * std::atan2(numerator, denominator);
}

double FastWindingNumber::multipoleSolidAngle(const Node& node, const Vec3& query, int order, double eps) {
    const Vec3 r = node.center - query;
    const double r2 = squaredNorm(r);
    if (r2 <= eps * eps) {
        return 0.0;
    }
    const double r1 = std::sqrt(r2);
    const double r3 = r2 * r1;
    const double r5 = r3 * r2;
    const double r7 = r5 * r2;
    const int expansion_order = std::max(0, std::min(2, order));

    double solid_angle = dot(node.area_vector, r) / r3;

    if (expansion_order >= 1) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                const double df = delta(j, k) / r3 - 3.0 * r[j] * r[k] / r5;
                solid_angle += node.first_moment[moment1Index(j, k)] * df;
            }
        }
    }

    if (expansion_order >= 2) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    const double ddf =
                        15.0 * r[j] * r[k] * r[l] / r7 -
                        3.0 * (
                            delta(j, k) * r[l] +
                            delta(j, l) * r[k] +
                            r[j] * delta(k, l)
                        ) / r5;
                    solid_angle += 0.5 * node.second_moment[moment2Index(j, k, l)] * ddf;
                }
            }
        }
    }

    return solid_angle;
}

}  // namespace fwn
