#pragma once

#include <string>
#include <vector>

#include "fwn/Mesh.h"

namespace fwn {

Mesh loadObj(const std::string& path);
std::vector<Vec3> loadPointFile(const std::string& path);

}  // namespace fwn
