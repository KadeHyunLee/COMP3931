#include "GridPartitioner.h"
#include <iostream>
#include <cfloat>
#ifdef _OPENMP
#include <omp.h>
#endif

GridPartitioner::GridPartitioner(float gridSize) : gridSize_(gridSize) {}

GridIndex GridPartitioner::computeIndex(const OpenMesh::Vec3f& point) const {
    return {
        static_cast<int>(point[0] / gridSize_),
        static_cast<int>(point[1] / gridSize_),
        static_cast<int>(point[2] / gridSize_)
    };
}

void GridPartitioner::mapVertices(const MyMesh& mesh,
    std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap) const {

    std::cout << "[Debug] Mapping vertices to grid...\n";
    std::cout << "[Debug] Total vertices in mesh: " << mesh.n_vertices() << "\n";
    
    // From OpenMP documentation
    std::mutex gridMapMutex;

    #pragma omp parallel
    {
        std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> privateMap;

        #pragma omp for
        for (int i = 0; i < mesh.n_vertices(); ++i) {
            auto vh = MyMesh::VertexHandle(i);
            auto point = mesh.point(vh);
            auto idx = computeIndex(point);
            privateMap[idx].push_back(vh);
        }

        #pragma omp critical
        {
            for (const auto& [idx, verts] : privateMap) {
                std::lock_guard<std::mutex> lock(gridMapMutex);
                gridMap[idx].insert(gridMap[idx].end(), verts.begin(), verts.end());
            }
        }
    }
}
