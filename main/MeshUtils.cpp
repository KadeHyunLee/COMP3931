#include "MeshUtils.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <iostream>
#include <unordered_set>

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
    float gridSize = avgSize * 0.5f;  
    std::cout << "[Debug] Auto-Calculated Grid Size: " << gridSize << "\n";

    return gridSize;
}

void extractSubMeshes(const MyMesh& original,
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap,
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes) {

    std::cout << "[Debug] Extracting submeshes (face-centric)...\n";
    float gridSize = calculateOptimalGridSize(original);

    std::unordered_map<GridIndex, MyMesh> localSubMeshes;
    std::unordered_map<GridIndex, std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle>> localVertexMaps;

    for (auto f_it = original.faces_begin(); f_it != original.faces_end(); ++f_it) {
        MyMesh::FaceHandle face = *f_it;

        std::vector<MyMesh::VertexHandle> originalVHs;
        OpenMesh::Vec3f faceCenter(0, 0, 0);

        for (auto fv_it = original.cfv_iter(face); fv_it.is_valid(); ++fv_it) {
            originalVHs.push_back(*fv_it);
            faceCenter += original.point(*fv_it);
        }

        if (originalVHs.size() != 3) continue;
        faceCenter /= 3.0f;

        GridIndex gridIdx = computeGridIndex(faceCenter, gridSize);
        auto& submesh = localSubMeshes[gridIdx];
        auto& vhandleMap = localVertexMaps[gridIdx];

        std::vector<MyMesh::VertexHandle> submeshVHs;
        std::unordered_set<int> unique;

        for (const auto& vh : originalVHs) {
            if (vhandleMap.find(vh) == vhandleMap.end()) {
                vhandleMap[vh] = submesh.add_vertex(original.point(vh));
            }
            auto newVH = vhandleMap[vh];
            submeshVHs.push_back(newVH);
            unique.insert(newVH.idx());
        }

        if (submeshVHs.size() == 3 && unique.size() == 3) {
            auto p0 = submesh.point(submeshVHs[0]);
            auto p1 = submesh.point(submeshVHs[1]);
            auto p2 = submesh.point(submeshVHs[2]);
            float area = OpenMesh::cross(p1 - p0, p2 - p0).norm() * 0.5f;

            if (area > 1e-8f) {
                submesh.add_face(submeshVHs);
            }
        }
    }

    // Copy results to output
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

void decimateSubMeshes(std::unordered_map<GridIndex, MyMesh>& subMeshes) {
    std::cout << "[Debug] Number of Submeshes: " << subMeshes.size() << "\n";
    std::cout << "[Debug] Performing Decimation on Submeshes...\n";

    int totalFacesBefore = 0;
    int totalVerticesBefore = 0;
    int totalFacesAfter = 0;
    int totalVerticesAfter = 0;

    for (auto& [gridIdx, submesh] : subMeshes) {
        std::cout << "[Debug] Processing Submesh at Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
        std::cout << " - Vertices: " << submesh.n_vertices() << ", Faces: " << submesh.n_faces() << "\n";
        int facesBefore = submesh.n_faces();
        if (facesBefore <= 3) {
            std::cout << "[Skip] Skipping Decimation for small Submesh at Grid ("
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z
                      << ") with " << facesBefore << " faces\n";
            continue;
        }
        int verticesBefore = submesh.n_vertices();
        
        totalFacesBefore += facesBefore;
        totalVerticesBefore += verticesBefore;
 
        if (facesBefore == 0) {
            std::cout << "[Skip] Skipping Decimation for empty Submesh at Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            continue;
        }

        submesh.request_vertex_status();
        submesh.request_halfedge_status();
        for (auto he_it = submesh.halfedges_begin(); he_it != submesh.halfedges_end(); ++he_it) {
            if (submesh.is_boundary(*he_it)) {
                auto from = submesh.from_vertex_handle(*he_it);
                auto to = submesh.to_vertex_handle(*he_it);
                submesh.status(from).set_locked(true);
                submesh.status(to).set_locked(true);
            }
        }

        OpenMesh::Decimater::DecimaterT<MyMesh> decimater(submesh);
        OpenMesh::Decimater::ModQuadricT<MyMesh>::Handle modQuadric;

        if (!decimater.add(modQuadric)) {
            std::cerr << "[Error] Failed to add Decimation module!\n";
            return;
        }

        decimater.initialize();
        auto& modQuadricRef = decimater.module(modQuadric);
        double maxError = facesBefore > 100 ? 1e-3 : 1e-5;
        modQuadricRef.set_max_err(maxError);

        int target_faces = std::max(static_cast<int>(facesBefore * 0.3), 3);
        decimater.decimate_to_faces(target_faces);
        std::cout << "[Debug] Decimation complete for Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
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

std::vector<MyMesh::VertexHandle> findBoundaryVertices(const MyMesh& mesh) {
    std::vector<MyMesh::VertexHandle> boundaryVertices;

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        if (mesh.is_boundary(*v_it)) {
            boundaryVertices.push_back(*v_it);
        }
    }

    std::cout << "[Debug] Found " << boundaryVertices.size() << " boundary vertices.\n";
    return boundaryVertices;
}

void seamFixBoundaryVertices(MyMesh& mesh, float mergeThreshold) {
    std::cout << "[Debug] Seam Fixing (Boundary Vertices) Started...\n";

    auto boundaryVertices = findBoundaryVertices(mesh);

    int mergeCount = 0;

    for (size_t i = 0; i < boundaryVertices.size(); ++i) {
        for (size_t j = i + 1; j < boundaryVertices.size(); ++j) {
            OpenMesh::Vec3f p1 = mesh.point(boundaryVertices[i]);
            OpenMesh::Vec3f p2 = mesh.point(boundaryVertices[j]);
            float distance = (p1 - p2).length();

            if (distance < mergeThreshold) {
                mesh.set_point(boundaryVertices[j], p1);
                mergeCount++;
            }
        }
    }

    std::cout << "[Debug] Seam Fixing Completed! Boundary vertices merged: " << mergeCount << "\n";

    mesh.garbage_collection();
}

void integrateSubMeshes(const std::unordered_map<GridIndex, MyMesh>& subMeshes, 
    MyMesh& finalMesh, 
    const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes, 
    const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes) { 
    std::cout << "[Debug] Integrating Submeshes...\n";

    std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> globalVertexMap;
    int mergedFaces = 0, mergedVertices = 0;
    int skippedSubmeshes = 0; 

    for (const auto& [gridIdx, submesh] : subMeshes) {
        
        if (fixedSubMeshes.find(gridIdx) != fixedSubMeshes.end()) {
            std::cout << "[Skip] Skipping fixed submesh from Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            continue;
        }

        if (emptySubMeshes.find(gridIdx) != emptySubMeshes.end()) {
            std::cerr << "[Warning] Skipping empty submesh from Grid (" 
                      << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
            skippedSubmeshes++;
            continue;
        }

        static constexpr float posEpsilon = 1e-6f;
        auto hashVec = [](const OpenMesh::Vec3f& v) {
            return std::make_tuple(
                std::round(v[0] / posEpsilon),
                std::round(v[1] / posEpsilon),
                std::round(v[2] / posEpsilon)
            );
        };

        std::unordered_map<std::tuple<int, int, int>, MyMesh::VertexHandle> positionMap;

        for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
            OpenMesh::Vec3f point = submesh.point(*v_it);
            auto key = hashVec(point);

            MyMesh::VertexHandle new_vhandle;

            auto it = positionMap.find(key);
            if (it != positionMap.end()) {
                new_vhandle = it->second;
            } else {
                new_vhandle = finalMesh.add_vertex(point);
                positionMap[key] = new_vhandle;
                mergedVertices++;
            }

            globalVertexMap[*v_it] = new_vhandle;
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

    // Perform seam-fixing after integrating submeshes
    // seamFixBoundaryVertices(finalMesh, 1e-5f);

    OpenMesh::IO::write_mesh(finalMesh, "final_integrated_mesh.ply");
    std::cout << "[Saved] Final Integrated Mesh saved as: final_integrated_mesh.ply\n";
}





