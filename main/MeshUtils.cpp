#include "MeshUtils.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModAspectRatioT.hh>
#include <iostream>
#include <unordered_set>
#include <omp.h> 

using std::cout;
using std::cerr;
using std::endl;


// Function to automatically calculate an optimal grid cell size based on the spatial bounds of the mesh
// it computes the average dimension of the mesh bounding box
float calculateOptimalGridSize(const MyMesh& mesh) {
    // Separate variables for reduction (OpenMP doesn't support Vec3f directly)
    float minX = FLT_MAX, minY = FLT_MAX, minZ = FLT_MAX;
    float maxX = -FLT_MAX, maxY = -FLT_MAX, maxZ = -FLT_MAX;

    #pragma omp parallel for reduction(min:minX, minY, minZ) reduction(max:maxX, maxY, maxZ)
    for (int i = 0; i < mesh.n_vertices(); ++i) {
        OpenMesh::Vec3f point = mesh.point(MyMesh::VertexHandle(i));
        minX = std::min(minX, point[0]);
        minY = std::min(minY, point[1]);
        minZ = std::min(minZ, point[2]);

        maxX = std::max(maxX, point[0]);
        maxY = std::max(maxY, point[1]);
        maxZ = std::max(maxZ, point[2]);
    }

    // Create global min/max vectors
    OpenMesh::Vec3f globalMin(minX, minY, minZ);
    OpenMesh::Vec3f globalMax(maxX, maxY, maxZ);

    // Compute average bounding box size
    float avgSize = (globalMax[0] - globalMin[0] +
                     globalMax[1] - globalMin[1] +
                     globalMax[2] - globalMin[2]) / 3.0f;

    float gridSize = avgSize * 0.2f;

    std::cout << "[Debug] Auto-Calculated Grid Size: " << gridSize << "\n";

    return gridSize;
}


float computeDynamicEpsilon(const MyMesh& mesh) {
    if (mesh.n_vertices() == 0) return 1e-6f;

    OpenMesh::Vec3f minPt = mesh.point(*mesh.vertices_begin());
    OpenMesh::Vec3f maxPt = minPt;

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f p = mesh.point(*v_it);
        for (int i = 0; i < 3; ++i) {
            minPt[i] = std::min(minPt[i], p[i]);
            maxPt[i] = std::max(maxPt[i], p[i]);
        }
    }

    float diagonal = (maxPt - minPt).length();
    return diagonal * 1e-4f;
}
