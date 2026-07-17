#include "fwn/IO.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fwn {
namespace {

std::string trimLeading(const std::string& s) {
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) {
        ++i;
    }
    return s.substr(i);
}

int parseObjIndex(const std::string& token, std::size_t vertex_count) {
    const std::size_t slash = token.find('/');
    const std::string head = slash == std::string::npos ? token : token.substr(0, slash);
    if (head.empty()) {
        throw std::runtime_error("OBJ face token has no vertex index: " + token);
    }

    const int raw = std::stoi(head);
    if (raw > 0) {
        const int idx = raw - 1;
        if (static_cast<std::size_t>(idx) >= vertex_count) {
            throw std::runtime_error("OBJ face index is outside the vertex array");
        }
        return idx;
    }

    if (raw < 0) {
        const long idx = static_cast<long>(vertex_count) + raw;
        if (idx < 0 || static_cast<std::size_t>(idx) >= vertex_count) {
            throw std::runtime_error("negative OBJ face index is outside the vertex array");
        }
        return static_cast<int>(idx);
    }

    throw std::runtime_error("OBJ indices are 1-based; zero is invalid");
}

}  // namespace

Mesh loadObj(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open OBJ file: " + path);
    }

    Mesh mesh;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        line = trimLeading(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string tag;
        iss >> tag;

        if (tag == "v") {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            if (!(iss >> x >> y >> z)) {
                throw std::runtime_error("invalid vertex at line " + std::to_string(line_number));
            }
            mesh.vertices.emplace_back(x, y, z);
        } else if (tag == "f") {
            std::vector<int> face;
            std::string token;
            while (iss >> token) {
                face.push_back(parseObjIndex(token, mesh.vertices.size()));
            }
            if (face.size() < 3) {
                throw std::runtime_error("face has fewer than three vertices at line " + std::to_string(line_number));
            }
            for (std::size_t i = 1; i + 1 < face.size(); ++i) {
                mesh.triangles.emplace_back(face[0], face[i], face[i + 1]);
            }
        }
    }

    if (mesh.empty()) {
        throw std::runtime_error("OBJ file contains no usable triangle mesh: " + path);
    }
    return mesh;
}

std::vector<Vec3> loadPointFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open point file: " + path);
    }

    std::vector<Vec3> points;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        line = trimLeading(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(iss >> x >> y >> z)) {
            throw std::runtime_error("invalid point at line " + std::to_string(line_number));
        }
        points.emplace_back(x, y, z);
    }
    return points;
}

}  // namespace fwn
