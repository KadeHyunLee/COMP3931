#include "SubmeshDecimator.h"
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModAspectRatioT.hh>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

void SubmeshDecimator::decimate(std::unordered_map<GridIndex, MyMesh>& subMeshes, float decimationRatio) {
    int totalFacesBefore = 0, totalVerticesBefore = 0;
    int totalFacesAfter = 0, totalVerticesAfter = 0;

    std::vector<GridIndex> gridIndices;
    for (const auto& [gridIdx, _] : subMeshes) {
        gridIndices.push_back(gridIdx);
    }

    #pragma omp parallel for reduction(+:totalFacesBefore, totalVerticesBefore, totalFacesAfter, totalVerticesAfter)
    for (int i = 0; i < static_cast<int>(gridIndices.size()); ++i) {
        GridIndex gridIdx = gridIndices[i];
        auto& submesh = subMeshes[gridIdx];

        int facesBefore = submesh.n_faces();
        if (facesBefore <= 3) {
            std::cout << "[Skip] Skipping Submesh at Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            continue;
        }

        int verticesBefore = submesh.n_vertices();

        totalFacesBefore += facesBefore;
        totalVerticesBefore += verticesBefore;

        submesh.request_vertex_status();
        submesh.request_halfedge_status();

        for (auto he_it = submesh.halfedges_begin(); he_it != submesh.halfedges_end(); ++he_it) {
            if (submesh.is_boundary(*he_it)) {
                submesh.status(submesh.from_vertex_handle(*he_it)).set_locked(true);
                submesh.status(submesh.to_vertex_handle(*he_it)).set_locked(true);
            }
        }

        OpenMesh::Decimater::DecimaterT<MyMesh> decimater(submesh);
        OpenMesh::Decimater::ModQuadricT<MyMesh>::Handle modQuadric;
        OpenMesh::Decimater::ModAspectRatioT<MyMesh>::Handle modAspect;

        if (!decimater.add(modQuadric) || !decimater.add(modAspect)) {
            std::cerr << "[Error] Failed to add decimation module!\n";
            continue;
        }

        float dynamicErr = (1.0f - decimationRatio) * 0.15f;
        int target_faces = std::max(static_cast<int>(facesBefore * decimationRatio), 3);

        decimater.module(modQuadric).set_max_err(dynamicErr);
        decimater.initialize();
        decimater.decimate_to_faces(target_faces);
        submesh.garbage_collection();

        totalFacesAfter += submesh.n_faces();
        totalVerticesAfter += submesh.n_vertices();
    }

    std::cout << "[Summary] Decimation Completed!\n";
    std::cout << " - Faces: " << totalFacesBefore << " -> " << totalFacesAfter << "\n";
    std::cout << " - Vertices: " << totalVerticesBefore << " -> " << totalVerticesAfter << "\n";

    float faceRate = 100.0f * (1.0f - float(totalFacesAfter) / totalFacesBefore);
    float vertexRate = 100.0f * (1.0f - float(totalVerticesAfter) / totalVerticesBefore);

    std::cout << " - Face Reduction: " << faceRate << "%\n";
    std::cout << " - Vertex Reduction: " << vertexRate << "%\n";
    std::cout << " - Target Decimation Ratio: " << decimationRatio << "\n";
}