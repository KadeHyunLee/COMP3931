#include "MeshIntegrator.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <atomic>
#include <mutex>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

void MeshIntegrator::integrate(
    const std::unordered_map<GridIndex, MyMesh>& subMeshes,
    MyMesh& finalMesh,
    const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes,
    const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes)
{
    std::cout << "[Debug] Integrating Submeshes...\n";

    std::atomic<int> mergedFaces{0}, mergedVertices{0}, skippedSubmeshes{0};
    std::mutex finalMeshMutex;

    static constexpr float posEpsilon = 1e-6f;
    
    auto hashVec = [](const OpenMesh::Vec3f& v) {
        return std::make_tuple(
            std::round(v[0] / posEpsilon),
            std::round(v[1] / posEpsilon),
            std::round(v[2] / posEpsilon)
        );
    };

    std::vector<GridIndex> keys;
    for (const auto& [gridIdx, _] : subMeshes) {
        keys.push_back(gridIdx);
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(keys.size()); ++i) {
        const auto& gridIdx = keys[i];
        const auto& submesh = subMeshes.at(gridIdx);

        if (fixedSubMeshes.count(gridIdx)) {
            std::cout << "[Skip] Skipping fixed submesh from Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            continue;
        }

        if (emptySubMeshes.count(gridIdx)) {
            std::cerr << "[Warning] Skipping empty submesh from Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            skippedSubmeshes++;
            continue;
        }

        std::unordered_map<std::tuple<int,int,int>, MyMesh::VertexHandle> positionMap;
        std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> localMap;

        for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
            OpenMesh::Vec3f point = submesh.point(*v_it);
            auto key = hashVec(point);

            if (positionMap.count(key) == 0) {
                MyMesh::VertexHandle newVH;
                {
                    std::lock_guard<std::mutex> lock(finalMeshMutex);
                    newVH = finalMesh.add_vertex(point);
                }
                positionMap[key] = newVH;
                mergedVertices++;
            }

            localMap[*v_it] = positionMap[key];
        }

        for (auto f_it = submesh.faces_begin(); f_it != submesh.faces_end(); ++f_it) {
            std::vector<MyMesh::VertexHandle> vhs;
            for (auto fv_it = submesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
                vhs.push_back(localMap[*fv_it]);
            }

            if (vhs.size() == 3) {
                std::lock_guard<std::mutex> lock(finalMeshMutex);
                finalMesh.add_face(vhs);
                mergedFaces++;
            } else {
                std::cerr << "[Warning] Skipping invalid face from Grid (" 
                          << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z 
                          << ") with " << vhs.size() << " vertices.\n";
            }
        }
    }

    std::cout << "[Debug] Integration Completed!\n";
    std::cout << "[Summary] Skipped Empty Submeshes: " << skippedSubmeshes.load() << "\n";

    OpenMesh::IO::write_mesh(finalMesh, "final_integrated_mesh.ply");
    std::cout << "[Saved] Final Integrated Mesh saved as: final_integrated_mesh.ply\n";
}