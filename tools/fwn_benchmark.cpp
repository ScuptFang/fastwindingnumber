#include "fwn/FastWindingNumber.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct Options {
    int u_segments = 64;
    int v_segments = 24;
    int query_count = 128;
    int repeats = 5;
    int seed = 1337;
    std::size_t leaf_size = 8;
    double major_radius = 1.0;
    double minor_radius = 0.35;
    double beta = 5.0;
    int expansion_order = 2;
};

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

std::vector<fwn::Vec3> makeQueries(const Options& options) {
    std::mt19937 rng(static_cast<std::mt19937::result_type>(options.seed));
    std::uniform_real_distribution<double> xy_dist(
        -(options.major_radius + options.minor_radius),
        options.major_radius + options.minor_radius
    );
    std::uniform_real_distribution<double> z_dist(-options.minor_radius, options.minor_radius);

    std::vector<fwn::Vec3> queries;
    queries.reserve(static_cast<std::size_t>(options.query_count));
    for (int i = 0; i < options.query_count; ++i) {
        queries.emplace_back(xy_dist(rng), xy_dist(rng), z_dist(rng));
    }
    return queries;
}

double parseDouble(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0') {
        throw std::runtime_error("invalid " + name + ": " + value);
    }
    return parsed;
}

int parseInt(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0') {
        throw std::runtime_error("invalid " + name + ": " + value);
    }
    return static_cast<int>(parsed);
}

std::size_t parseSize(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0' || parsed == 0) {
        throw std::runtime_error("invalid " + name + ": " + value);
    }
    return static_cast<std::size_t>(parsed);
}

void printUsage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " [options]\n"
        << "  --queries N           Query count, default 128\n"
        << "  --repeats N           Timing repeats, default 5\n"
        << "  --seed N              Random seed, default 1337\n"
        << "  --u-segments N        Torus major-ring segments, default 64\n"
        << "  --v-segments N        Torus tube segments, default 24\n"
        << "  --leaf-size N         BVH leaf triangle count, default 8\n"
        << "  --beta VALUE          Fast far-field factor, default 5.0\n"
        << "  --expansion-order N   Fast expansion order 0, 1, or 2; default 2\n";
}

Options parseArgs(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--queries" && i + 1 < argc) {
            options.query_count = parseInt(argv[++i], "queries");
        } else if (arg == "--repeats" && i + 1 < argc) {
            options.repeats = parseInt(argv[++i], "repeats");
        } else if (arg == "--seed" && i + 1 < argc) {
            options.seed = parseInt(argv[++i], "seed");
        } else if (arg == "--u-segments" && i + 1 < argc) {
            options.u_segments = parseInt(argv[++i], "u segments");
        } else if (arg == "--v-segments" && i + 1 < argc) {
            options.v_segments = parseInt(argv[++i], "v segments");
        } else if (arg == "--leaf-size" && i + 1 < argc) {
            options.leaf_size = parseSize(argv[++i], "leaf size");
        } else if (arg == "--beta" && i + 1 < argc) {
            options.beta = parseDouble(argv[++i], "beta");
        } else if (arg == "--expansion-order" && i + 1 < argc) {
            options.expansion_order = parseInt(argv[++i], "expansion order");
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown or incomplete option: " + arg);
        }
    }

    if (options.u_segments < 3 || options.v_segments < 3) {
        throw std::runtime_error("torus segment counts must be at least 3");
    }
    if (options.query_count <= 0 || options.repeats <= 0) {
        throw std::runtime_error("queries and repeats must be positive");
    }
    if (options.expansion_order < 0 || options.expansion_order > 2) {
        throw std::runtime_error("expansion order must be 0, 1, or 2");
    }
    return options;
}

template <typename Func>
double timeMilliseconds(Func&& fn) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseArgs(argc, argv);
        const fwn::Mesh mesh = makeTorus(
            options.major_radius,
            options.minor_radius,
            options.u_segments,
            options.v_segments
        );
        const std::vector<fwn::Vec3> queries = makeQueries(options);

        double build_ms = 0.0;
        fwn::FastWindingNumber* raw = nullptr;
        build_ms = timeMilliseconds([&]() {
            raw = new fwn::FastWindingNumber(mesh, options.leaf_size);
        });
        const std::unique_ptr<fwn::FastWindingNumber> fwn(raw);

        fwn::QueryOptions fast_options;
        fast_options.fast = true;
        fast_options.beta = options.beta;
        fast_options.expansion_order = options.expansion_order;

        std::vector<double> exact_values(queries.size());
        std::vector<double> fast_values(queries.size());
        volatile double exact_checksum = 0.0;
        volatile double fast_checksum = 0.0;

        const double exact_ms = timeMilliseconds([&]() {
            for (int repeat = 0; repeat < options.repeats; ++repeat) {
                for (std::size_t i = 0; i < queries.size(); ++i) {
                    exact_values[i] = fwn->exactWindingNumber(queries[i]);
                    exact_checksum += exact_values[i];
                }
            }
        });

        const double fast_ms = timeMilliseconds([&]() {
            for (int repeat = 0; repeat < options.repeats; ++repeat) {
                for (std::size_t i = 0; i < queries.size(); ++i) {
                    fast_values[i] = fwn->windingNumber(queries[i], fast_options);
                    fast_checksum += fast_values[i];
                }
            }
        });

        double worst_error = 0.0;
        double total_error = 0.0;
        for (std::size_t i = 0; i < queries.size(); ++i) {
            const double error = std::abs(fast_values[i] - exact_values[i]);
            worst_error = std::max(worst_error, error);
            total_error += error;
        }

        const double total_queries = static_cast<double>(options.query_count * options.repeats);
        const double mean_error = total_error / static_cast<double>(queries.size());
        const double exact_us = exact_ms * 1000.0 / total_queries;
        const double fast_us = fast_ms * 1000.0 / total_queries;

        std::cout << std::setprecision(6) << std::fixed;
        std::cout << "mesh_vertices " << mesh.vertices.size() << '\n';
        std::cout << "mesh_triangles " << mesh.triangles.size() << '\n';
        std::cout << "bvh_nodes " << fwn->nodeCount() << '\n';
        std::cout << "queries " << options.query_count << '\n';
        std::cout << "repeats " << options.repeats << '\n';
        std::cout << "leaf_size " << options.leaf_size << '\n';
        std::cout << "beta " << options.beta << '\n';
        std::cout << "expansion_order " << options.expansion_order << '\n';
        std::cout << "build_ms " << build_ms << '\n';
        std::cout << "exact_total_ms " << exact_ms << '\n';
        std::cout << "fast_total_ms " << fast_ms << '\n';
        std::cout << "exact_us_per_query " << exact_us << '\n';
        std::cout << "fast_us_per_query " << fast_us << '\n';
        std::cout << "speedup " << (fast_ms > 0.0 ? exact_ms / fast_ms : 0.0) << '\n';
        std::cout << "worst_abs_error " << worst_error << '\n';
        std::cout << "mean_abs_error " << mean_error << '\n';
        std::cout << "exact_checksum " << exact_checksum << '\n';
        std::cout << "fast_checksum " << fast_checksum << '\n';
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
