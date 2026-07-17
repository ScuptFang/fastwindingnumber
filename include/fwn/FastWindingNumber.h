#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "fwn/Mesh.h"

namespace fwn {

struct QueryOptions {
    bool fast = true;
    double beta = 3.5;
    int expansion_order = 2;
    double singular_epsilon = 1e-14;
};

class FastWindingNumber {
public:
    explicit FastWindingNumber(const Mesh& mesh, std::size_t leaf_size = 16);

    double windingNumber(const Vec3& query, const QueryOptions& options = QueryOptions()) const;
    double exactWindingNumber(const Vec3& query, double singular_epsilon = 1e-14) const;
    std::vector<double> windingNumbers(
        const std::vector<Vec3>& queries,
        const QueryOptions& options = QueryOptions()
    ) const;

    std::size_t triangleCount() const { return triangles_.size(); }
    std::size_t nodeCount() const { return nodes_.size(); }
    std::size_t leafSize() const { return leaf_size_; }

private:
    struct Aabb {
        Vec3 min;
        Vec3 max;
        bool valid = false;

        void expand(const Vec3& p);
        void expand(const Aabb& box);
        Vec3 center() const;
        Vec3 extent() const;
    };

    struct TriangleData {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;
        Vec3 centroid;
        Vec3 area_vector;
        double area = 0.0;
        Aabb bounds;
    };

    struct Node {
        Aabb bounds;
        Vec3 center;
        Vec3 area_vector;
        std::array<double, 9> first_moment{};
        std::array<double, 27> second_moment{};
        double total_area = 0.0;
        double radius = 0.0;
        std::size_t begin = 0;
        std::size_t end = 0;
        int left = -1;
        int right = -1;

        bool isLeaf() const {
            return left < 0 && right < 0;
        }
    };

    std::size_t buildNode(std::size_t begin, std::size_t end);
    double solidAngleRecursive(std::size_t node_index, const Vec3& query, const QueryOptions& options) const;
    double exactSolidAngleRange(std::size_t begin, std::size_t end, const Vec3& query, double eps) const;

    static double triangleSolidAngle(const TriangleData& tri, const Vec3& query, double eps);
    static double multipoleSolidAngle(const Node& node, const Vec3& query, int order, double eps);

    std::vector<TriangleData> triangles_;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;
    std::size_t root_ = 0;
    std::size_t leaf_size_ = 16;
};

}  // namespace fwn
