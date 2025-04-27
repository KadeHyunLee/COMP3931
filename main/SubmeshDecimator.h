#ifndef SUBMESH_DECIMATOR_H
#define SUBMESH_DECIMATOR_H

#include "MeshUtils.h"
#include <unordered_map>

class SubmeshDecimator {
public:
    // decimationRatio: 0.0 ~ 1.0
    static void decimate(std::unordered_map<GridIndex, MyMesh>& subMeshes, float decimationRatio);
    void improveBoundary(std::unordered_map<GridIndex, MyMesh>& subMeshes);
};

#endif // SUBMESH_DECIMATOR_H