#include "fwn/FastWindingNumber.h"
#include "fwn/IO.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct CliOptions {
    std::string mesh_path;
    std::string points_path;
    std::string output_path;
    bool fast = true;
    bool classify = false;
    double beta = 3.5;
    double threshold = 0.5;
    int expansion_order = 2;
    std::size_t leaf_size = 16;
};

void printUsage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " mesh.obj points.xyz [options]\n"
        << "\n"
        << "Options:\n"
        << "  --fast                 Use BVH dipole approximation (default)\n"
        << "  --exact                Sum exact triangle solid angles\n"
        << "  --beta VALUE           Far-field acceptance factor, default 3.5\n"
        << "  --expansion-order N    Multipole expansion order 0, 1, or 2; default 2\n"
        << "  --leaf-size VALUE      BVH leaf triangle count, default 16\n"
        << "  --classify             Add inside/outside column using abs(w) >= threshold\n"
        << "  --threshold VALUE      Classification threshold, default 0.5\n"
        << "  --output PATH          Write output to PATH instead of stdout\n";
}

double parseDouble(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0') {
        throw std::runtime_error("invalid " + name + ": " + value);
    }
    return parsed;
}

std::size_t parseSize(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0' || parsed == 0) {
        throw std::runtime_error("invalid " + name + ": " + value);
    }
    return static_cast<std::size_t>(parsed);
}

int parseExpansionOrder(const std::string& value) {
    char* end = nullptr;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0' || parsed < 0 || parsed > 2) {
        throw std::runtime_error("expansion order must be 0, 1, or 2: " + value);
    }
    return static_cast<int>(parsed);
}

CliOptions parseArgs(int argc, char** argv) {
    if (argc < 3) {
        printUsage(argv[0]);
        throw std::runtime_error("missing mesh or point file");
    }

    CliOptions opts;
    opts.mesh_path = argv[1];
    opts.points_path = argv[2];

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--fast") {
            opts.fast = true;
        } else if (arg == "--exact") {
            opts.fast = false;
        } else if (arg == "--classify") {
            opts.classify = true;
        } else if (arg == "--beta" && i + 1 < argc) {
            opts.beta = parseDouble(argv[++i], "beta");
        } else if (arg == "--expansion-order" && i + 1 < argc) {
            opts.expansion_order = parseExpansionOrder(argv[++i]);
        } else if (arg == "--threshold" && i + 1 < argc) {
            opts.threshold = parseDouble(argv[++i], "threshold");
        } else if (arg == "--leaf-size" && i + 1 < argc) {
            opts.leaf_size = parseSize(argv[++i], "leaf size");
        } else if (arg == "--output" && i + 1 < argc) {
            opts.output_path = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown or incomplete option: " + arg);
        }
    }
    return opts;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions cli = parseArgs(argc, argv);
        const fwn::Mesh mesh = fwn::loadObj(cli.mesh_path);
        const std::vector<fwn::Vec3> points = fwn::loadPointFile(cli.points_path);

        fwn::FastWindingNumber fwn(mesh, cli.leaf_size);
        fwn::QueryOptions query;
        query.fast = cli.fast;
        query.beta = cli.beta;
        query.expansion_order = cli.expansion_order;

        std::ofstream file_out;
        std::ostream* out = &std::cout;
        if (!cli.output_path.empty()) {
            file_out.open(cli.output_path);
            if (!file_out) {
                throw std::runtime_error("failed to open output file: " + cli.output_path);
            }
            out = &file_out;
        }

        *out << std::setprecision(17);
        for (const fwn::Vec3& point : points) {
            const double w = fwn.windingNumber(point, query);
            *out << point << ' ' << w;
            if (cli.classify) {
                *out << ' ' << (std::abs(w) >= cli.threshold ? 1 : 0);
            }
            *out << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
