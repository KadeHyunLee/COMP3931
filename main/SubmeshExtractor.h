#ifndef SUBMESH_EXTRACTOR_H
#define SUBMESH_EXTRACTOR_H

#include "MeshUtils.h"
#include <unordered_map>

class SubmeshExtractor {
public:
    static void extract(const MyMesh& original,
        const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap,
        std::unordered_map<GridIndex, MyMesh>& subMeshes,
        std::unordered_map<GridIndex, MyMesh>& emptySubMeshes,
        float& outposEpsilon);
};

#endif // SUBMESH_EXTRACTOR_H