#include "MeshUtils.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

// Function to load a mesh from the given file into the provided mesh object.
// Returns true if the mesh is successfully loaded, false otherwise.
bool loadMesh(const std::string& filename, MyMesh& mesh) {
    // Attempt to read the mesh from the given filename and load it 
    if (!OpenMesh::IO::read_mesh(mesh, filename)) {
        std::cerr << "[Error] Failed to load mesh: " << filename << "\n";
        std::cerr << "[Debug] Please check the file path and format. Supported formats include: .obj, .ply, .stl\n";
        // Return false if file was not load
        return false;
    }
    std::cout << "[Debug] Loaded mesh: " << filename << " with "
              << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";
    // Return true to indicate successful loading.
    return true;
}

// The point is divided by the grid cell size to determine which cell it belongs to.
GridIndex computeGridIndex(const OpenMesh::Vec3f& point, float gridSize) {
    GridIndex idx {
        static_cast<int>(point[0] / gridSize),// Compute grid index along X-axis
        static_cast<int>(point[1] / gridSize),// Compute grid index along Y-axis
        static_cast<int>(point[2] / gridSize)// Compute grid index along Z-axis
    };
    /*std::cout << "[Debug] Computed GridIndex(" << idx.x << ", " << idx.y << ", " << idx.z
              << ") from point (" << point[0] << ", " << point[1] << ", " << point[2]
              << ") with gridSize " << gridSize << "\n";
              */
    return idx;
}

// Function to map each vertex of the mesh into its corresponding grid cell.
void mapVerticesToGrid(const MyMesh& mesh, std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, float gridSize) {
    std::cout << "[Debug] Mapping vertices to grid...\n";
    std::cout << "[Debug] Total vertices in mesh: " << mesh.n_vertices() << "\n";

    // Initialize min and max indices to track the range of grid cells used
    int minIdxX = INT_MAX, minIdxY = INT_MAX, minIdxZ = INT_MAX;
    int maxIdxX = INT_MIN, maxIdxY = INT_MIN, maxIdxZ = INT_MIN;

    // Iterate over each vertex in the mesh
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        // Get 3D Position of the current vertex
        OpenMesh::Vec3f point = mesh.point(*v_it);
        // Compute gridcell to belong points
        GridIndex gridIdx = computeGridIndex(point, gridSize);

        // Add the current vertex handle to the corresponding grid cell
        gridMap[gridIdx].push_back(*v_it);
        //std::cout << "[Debug] Added vertex handle " << (*v_it).idx()
          //        << " to GridIndex(" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";

        // Track the minimum and maximum indices in X, Y, Z
        minIdxX = std::min(minIdxX, gridIdx.x);
        minIdxY = std::min(minIdxY, gridIdx.y);
        minIdxZ = std::min(minIdxZ, gridIdx.z);

        maxIdxX = std::max(maxIdxX, gridIdx.x);
        maxIdxY = std::max(maxIdxY, gridIdx.y);
        maxIdxZ = std::max(maxIdxZ, gridIdx.z);
    }

    std::cout << "[Debug] Current min/max index ranges: "
              << "X(" << minIdxX << " ~ " << maxIdxX << "), "
              << "Y(" << minIdxY << " ~ " << maxIdxY << "), "
              << "Z(" << minIdxZ << " ~ " << maxIdxZ << ")\n";
    std::cout << "[Debug] Finished mapping vertices to grid.\n";
    std::cout << "[Debug] Grid Index Range: X(" << minIdxX << " ~ " << maxIdxX 
              << "), Y(" << minIdxY << " ~ " << maxIdxY 
              << "), Z(" << minIdxZ << " ~ " << maxIdxZ << ")\n";
    
    std::cout << "[Debug] Grid Distribution: \n";
    for (const auto& cell : gridMap) {
        std::cout << "  Grid (" << cell.first.x << ", " 
                << cell.first.y << ", " << cell.first.z 
                << ") -> " << cell.second.size() << " vertices\n";
    }
}

// Function to automatically calculate an optimal grid cell size based on the spatial bounds of the mesh
// it computes the average dimension of the mesh bounding box
float calculateOptimalGridSize(const MyMesh& mesh) {
    // Initialize min and max bounds
    OpenMesh::Vec3f minBounds(FLT_MAX, FLT_MAX, FLT_MAX);
    OpenMesh::Vec3f maxBounds(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    // Iterate over each vertex to find the min and max bounds
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        // Get the position of current vertex
        OpenMesh::Vec3f point = mesh.point(*v_it);
        // For each axis, update min and max bounds
        for (int i = 0; i < 3; i++) {
            minBounds[i] = std::min(minBounds[i], point[i]);
            maxBounds[i] = std::max(maxBounds[i], point[i]);
        }
    }

    // Compute the size of bounding box
    float avgSize = (maxBounds[0] - minBounds[0] +
                     maxBounds[1] - minBounds[1] +
                     maxBounds[2] - minBounds[2]) / 3.0f;
    
    // Calculate the average isze across all three axes
    float gridSize = avgSize * 0.05f;  
    std::cout << "[Debug] Auto-Calculated Grid Size: " << gridSize << "\n";

    return gridSize;
}

/*
void removeEmptyGrids(std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
                      const MyMesh& mesh, int minThreshold) {
    std::cout << "[Debug] Before Removing: " << gridMap.size() << " grids\n";
    
    auto it = gridMap.begin();
    while (it != gridMap.end()) {
        const auto& vertices = it->second;

        if (vertices.size() < minThreshold) {
            bool hasConnectedFaces = false;

            for (const auto& v : vertices) {
                if (mesh.valence(v) > 0) { 
                    hasConnectedFaces = true;
                    break;
                }
            }

            if (!hasConnectedFaces) {
                std::cout << "[Debug] Removing Grid (" << it->first.x << ", " 
                          << it->first.y << ", " << it->first.z 
                          << ") -> " << vertices.size() << " vertices (No Faces Attached)\n";

                it = gridMap.erase(it);
                continue;
            }
        }
        ++it;
    }

    std::cout << "[Debug] After Removing: " << gridMap.size() << " grids remain\n";
}
*/

void extractSubMeshes(const MyMesh& original, 
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes) { 
    std::cout << "[Debug] Extracting submeshes from grids...\n";

    // Counter for vaild submesh
    int subMeshCount = 0;
    // Counter for submeshes without faces
    int emptySubMeshCount = 0;
    
    // Iterate through each grid cell in the grid map
    for (const auto& cell : gridMap) {
        const GridIndex& gridIdx = cell.first;
        const auto& vertices = cell.second;

        MyMesh submesh; // New submesh for the grid
        // Map to track old vertex handles to new one in submesh
        std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> vhandleMap;

        // STEP 1 : Add vertices to submesh
        for (const auto& v : vertices) {
            // If vertex has not been added yet
            if (vhandleMap.find(v) == vhandleMap.end()) {
                // Get postion from the orginal mesh
                OpenMesh::Vec3f point = original.point(v);
                // Add vertex to submesh
                MyMesh::VertexHandle new_vhandle = submesh.add_vertex(point);
                // Map old one to new one
                vhandleMap[v] = new_vhandle;
            }
        }

        // Track how many faces are added for this submesh
        int faceAddCount = 0;

        for (const auto& v : vertices) {
            // Iterate over face connected to vertex v in original mesh
            for (MyMesh::ConstVertexFaceIter vf_it = original.cvf_iter(v); vf_it.is_valid(); ++vf_it) {
                // Current face
                MyMesh::FaceHandle face = *vf_it;
                // Collect new vertex for face
                std::vector<MyMesh::VertexHandle> face_vhandles;
                // Flag to check if all vertices are inside vhandleMap
                bool validFace = true;
                
                // Iterate over vertices of face
                for (MyMesh::ConstFaceVertexIter fv_it = original.cfv_iter(face); fv_it.is_valid(); ++fv_it) {
                    // If the vertex is in the submesh, add the new vertex
                    if (vhandleMap.find(*fv_it) != vhandleMap.end()) {
                        face_vhandles.push_back(vhandleMap[*fv_it]);
                    } else {
                        // At least one vertex is missing in this submesh, skip the face
                        validFace = false;
                        break;
                    }
                }
                // If all vertices are vaild, and it can be triangle
                if (validFace && face_vhandles.size() == 3) {
                    // Add face to submesh
                    MyMesh::FaceHandle newFace = submesh.add_face(face_vhandles);
                    if (newFace.is_valid()) {
                        faceAddCount++;
                    }
                }
            }
            }

        std::cout << "[Debug] Submesh for Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z
                  << ") -> " << submesh.n_vertices() << " vertices, " 
                  << submesh.n_faces() << " faces\n";
        
        // STEP 3: Map submesh into submeshes map
        if (submesh.n_vertices() > 0) {
            subMeshes[gridIdx] = submesh;
            subMeshCount++;

            if (submesh.n_faces() == 0) {
                std::cerr << "[Warning] Submesh for Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z
                          << ") has 0 faces but " << submesh.n_vertices() << " vertices!\n";
                emptySubMeshes[gridIdx] = submesh; 
                emptySubMeshCount++;
            }
        }
    }

    std::cout << "[Debug] Extracted " << subMeshCount << " valid submeshes.\n";
    std::cerr << "[Warning] Found " << emptySubMeshCount << " submeshes with 0 faces!\n";
}

void decimateSubMeshes(std::unordered_map<GridIndex, MyMesh>& subMeshes) {
    std::cout << "[Debug] Performing Decimation on Submeshes...\n";

    int totalFacesBefore = 0;
    int totalVerticesBefore = 0;
    int totalFacesAfter = 0;
    int totalVerticesAfter = 0;

    for (auto& [gridIdx, submesh] : subMeshes) {
        int facesBefore = submesh.n_faces();
        int verticesBefore = submesh.n_vertices();
        
        totalFacesBefore += facesBefore;
        totalVerticesBefore += verticesBefore;

        if (facesBefore == 0) {
            std::cout << "[Skip] Skipping Decimation for empty Submesh at Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            continue;
        }

        OpenMesh::Decimater::DecimaterT<MyMesh> decimater(submesh);
        OpenMesh::Decimater::ModQuadricT<MyMesh>::Handle modQuadric;

        if (!decimater.add(modQuadric)) {
            std::cerr << "[Error] Failed to add Decimation module!\n";
            return;
        }

        decimater.initialize();
        auto& modQuadricRef = decimater.module(modQuadric);
        double maxError = facesBefore > 100 ? 1e-2 : 1e-4;
        modQuadricRef.set_max_err(maxError);

        int target_faces = std::max(static_cast<int>(facesBefore * 0.3), 3);
        decimater.decimate_to_faces(target_faces);
        submesh.garbage_collection();

        int facesAfter = submesh.n_faces();
        int verticesAfter = submesh.n_vertices();

        totalFacesAfter += facesAfter;
        totalVerticesAfter += verticesAfter;

        std::cout << "[Debug] Decimated Submesh at Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z
                  << ") -> " << facesAfter << " faces remaining (Before: " << facesBefore << " faces, "
                  << "Vertices: " << verticesBefore << " → " << verticesAfter << ")\n";

        if (facesBefore == facesAfter) {
            std::cerr << "[Warning] No faces were removed during Decimation! Check Decimation settings.\n";
        }
    }

    std::cout << "=============================\n";
    std::cout << "[Summary] Decimation Completed!\n";
    std::cout << " - Total Faces Before: " << totalFacesBefore << "\n";
    std::cout << " - Total Faces After: " << totalFacesAfter << "\n";
    std::cout << " - Total Vertices Before: " << totalVerticesBefore << "\n";
    std::cout << " - Total Vertices After: " << totalVerticesAfter << "\n";
    std::cout << "=============================\n";
}

void integrateSubMeshes(const std::unordered_map<GridIndex, MyMesh>& subMeshes, 
    MyMesh& finalMesh, 
    const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes, 
    const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes) { 
    std::cout << "[Debug] Integrating Submeshes...\n";

    std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> globalVertexMap;
    int mergedFaces = 0, mergedVertices = 0;
    int skippedSubmeshes = 0; // 스킵된 서브메쉬 개수

    for (const auto& [gridIdx, submesh] : subMeshes) {
        // ✅ 보정된 서브메쉬(`fixedSubMeshes`)에 포함된 그리드는 스킵
        if (fixedSubMeshes.find(gridIdx) != fixedSubMeshes.end()) {
            std::cout << "[Skip] Skipping fixed submesh from Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            continue;
        }

        // ✅ 빈 서브메쉬는 메인에서 처리하므로 여기서는 스킵만 함
        if (emptySubMeshes.find(gridIdx) != emptySubMeshes.end()) {
            std::cerr << "[Warning] Skipping empty submesh from Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            skippedSubmeshes++;
            continue;
        }

        for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
            OpenMesh::Vec3f point = submesh.point(*v_it);
            MyMesh::VertexHandle new_vhandle = finalMesh.add_vertex(point);
            globalVertexMap[*v_it] = new_vhandle;
            mergedVertices++;
        }

        for (auto f_it = submesh.faces_begin(); f_it != submesh.faces_end(); ++f_it) {
            std::vector<MyMesh::VertexHandle> face_vhandles;
            for (auto fv_it = submesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
                face_vhandles.push_back(globalVertexMap[*fv_it]);
            }

            if (face_vhandles.size() == 3) {
                finalMesh.add_face(face_vhandles);
                mergedFaces++;
            } else {
                std::cerr << "[Warning] Skipping invalid face from Grid (" 
                          << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z 
                          << ") with " << face_vhandles.size() << " vertices.\n";
            }
        }
    }

    std::cout << "[Debug] Integration Completed!\n";
    std::cout << "[Summary] Skipped Empty Submeshes: " << skippedSubmeshes << "\n";

    OpenMesh::IO::write_mesh(finalMesh, "final_integrated_mesh.ply");
    std::cout << "[Saved] Final Integrated Mesh saved as: final_integrated_mesh.ply\n";
}

