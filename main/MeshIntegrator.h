#ifndef MESH_INTEGRATOR_H
#define MESH_INTEGRATOR_H

#include "MeshUtils.h"
#include <unordered_map>

class MeshIntegrator {
public:
    static void integrate(
        const std::unordered_map<GridIndex, MyMesh>& subMeshes,
        MyMesh& finalMesh,
        const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes,
        const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes,
        float posEpsilon);
};

#endif // MESH_INTEGRATOR_H