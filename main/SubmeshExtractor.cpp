#include "SubmeshExtractor.h"
#include "GridPartitioner.h"
#include <iostream>
#include <unordered_set>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

void SubmeshExtractor::extract(const MyMesh& original,
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap,
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes)
{
    std::cout << "[Debug] Extracting submeshes (face-centric with OpenMP)...\n";

    float gridSize = calculateOptimalGridSize(original);
    GridPartitioner partitioner(gridSize);

    std::vector<MyMesh::FaceHandle> allFaces;
    for (auto f_it = original.faces_begin(); f_it != original.faces_end(); ++f_it)
        allFaces.push_back(*f_it);

    std::unordered_map<GridIndex, MyMesh> localSubMeshes;
    std::unordered_map<GridIndex, std::unordered_map<std::tuple<int, int, int>, MyMesh::VertexHandle>> globalVertexMaps;

    std::vector<std::tuple<GridIndex, MyMesh, std::unordered_map<std::tuple<int,int,int>, MyMesh::VertexHandle>>> mergeQueue;
    static constexpr float posEpsilon = 1e-6f;

    #pragma omp parallel
    {
        std::unordered_map<GridIndex, MyMesh> threadSubMeshes;
        std::unordered_map<GridIndex, std::unordered_map<std::tuple<int, int, int>, MyMesh::VertexHandle>> threadVertexMaps;

        auto hashVec = [](const OpenMesh::Vec3f& v) {
            return std::make_tuple(
                static_cast<int>(v[0] / posEpsilon),
                static_cast<int>(v[1] / posEpsilon),
                static_cast<int>(v[2] / posEpsilon)
            );
        };

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < static_cast<int>(allFaces.size()); ++i) {
            auto face = allFaces[i];
            std::vector<MyMesh::VertexHandle> originalVHs;
            OpenMesh::Vec3f faceCenter(0, 0, 0);

            for (auto fv_it = original.cfv_iter(face); fv_it.is_valid(); ++fv_it) {
                originalVHs.push_back(*fv_it);
                faceCenter += original.point(*fv_it);
            }

            if (originalVHs.size() != 3) continue;
            faceCenter /= 3.0f;

            GridIndex gridIdx = partitioner.computeIndex(faceCenter);
            auto& submesh = threadSubMeshes[gridIdx];
            auto& vhandleMap = threadVertexMaps[gridIdx];

            std::vector<MyMesh::VertexHandle> submeshVHs;
            std::unordered_set<int> uniqueIndices;

            for (const auto& vh : originalVHs) {
                auto point = original.point(vh);
                auto key = hashVec(point);
                if (vhandleMap.find(key) == vhandleMap.end()) {
                    vhandleMap[key] = submesh.add_vertex(point);
                }
                auto newVH = vhandleMap[key];
                submeshVHs.push_back(newVH);
                uniqueIndices.insert(newVH.idx());
            }

            if (submeshVHs.size() != 3 || uniqueIndices.size() != 3) continue;

            auto p0 = submesh.point(submeshVHs[0]);
            auto p1 = submesh.point(submeshVHs[1]);
            auto p2 = submesh.point(submeshVHs[2]);
            float area = OpenMesh::cross(p1 - p0, p2 - p0).norm() * 0.5f;

            if (area > 1e-8f) {
                submesh.add_face(submeshVHs);
            }
        }

        #pragma omp critical
        {
            for (auto& [gridIdx, mesh] : threadSubMeshes) {
                mergeQueue.emplace_back(gridIdx, std::move(mesh), std::move(threadVertexMaps[gridIdx]));
            }
        }
    }

    auto hashVec = [](const OpenMesh::Vec3f& v) {
        return std::make_tuple(
            static_cast<int>(v[0] / posEpsilon),
            static_cast<int>(v[1] / posEpsilon),
            static_cast<int>(v[2] / posEpsilon)
        );
    };

    for (auto& [gridIdx, mesh, vertexMap] : mergeQueue) {
        if (localSubMeshes.find(gridIdx) == localSubMeshes.end()) {
            localSubMeshes[gridIdx] = std::move(mesh);
            globalVertexMaps[gridIdx] = std::move(vertexMap);
        } else {
            auto& dst = localSubMeshes[gridIdx];
            auto& globalMap = globalVertexMaps[gridIdx];
            std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> handleMap;

            for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
                auto point = mesh.point(*v_it);
                auto key = hashVec(point);
                if (globalMap.find(key) == globalMap.end()) {
                    globalMap[key] = dst.add_vertex(point);
                }
                handleMap[*v_it] = globalMap[key];
            }

            for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
                std::vector<MyMesh::VertexHandle> vhs;
                for (auto fv_it = mesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
                    vhs.push_back(handleMap[*fv_it]);
                }
                dst.add_face(vhs);
            }
        }
    }

    for (auto& [gridIdx, mesh] : localSubMeshes) {
        if (mesh.n_faces() > 0) {
            subMeshes[gridIdx] = mesh;
        } else {
            emptySubMeshes[gridIdx] = mesh;
        }
    }

    std::cout << "[Debug] Extracted " << subMeshes.size() << " submeshes.\n";
    std::cerr << "[Warning] Found " << emptySubMeshes.size() << " empty submeshes.\n";
}