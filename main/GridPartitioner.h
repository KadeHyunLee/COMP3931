#ifndef GRID_PARTITIONER_H
#define GRID_PARTITIONER_H

#include "MeshUtils.h"
#include <unordered_map>
#include <mutex>

class GridPartitioner {
public:
    GridPartitioner(float gridSize);

    GridIndex computeIndex(const OpenMesh::Vec3f& point) const;

    void mapVertices(const MyMesh& mesh,
        std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap) const;

private:
    float gridSize_;
    static constexpr float posEpsilon = 1e-6f;
};

#endif // GRID_PARTITIONER_H