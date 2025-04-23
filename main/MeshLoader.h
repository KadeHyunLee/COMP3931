#ifndef MESHLOADER_H
#define MESHLOADER_H

#include "MeshUtils.h"
#include <string>

class MeshLoader {
public:
    // Loads a mesh from file into the provided mesh object
    static bool load(const std::string& filename, MyMesh& mesh);
};

#endif // MESHLOADER_H