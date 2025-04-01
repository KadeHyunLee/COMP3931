#include "MeshUtils.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModAspectRatioT.hh>
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

    return idx;
}

// Function to map each vertex of the mesh into its corresponding grid cell.
void mapVerticesToGrid(const MyMesh& mesh, std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, float gridSize) {
    // Debug: Start mapping vertices to grid
    std::cout << "[Debug] Mapping vertices to grid...\n";
    // Debug: Output total number of vertices in the mesh
    std::cout << "[Debug] Total vertices in mesh: " << mesh.n_vertices() << "\n";

    // Mutex to protect shared gridMap during merging from parallel threads
    std::mutex gridMapMutex;
    
    // Begin OpenMP parallel region to process vertices concurrently
    #pragma omp parallel
    {
        // Each thread maintains its own private map to collect vertex handles per grid cell
        std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> privateGridMap;

        // Distribute the vertex loop among threads using OpenMP for directive
        #pragma omp for
        for (int i = 0; i < mesh.n_vertices(); ++i) {
            // Retrieve vertex handle for current index
            auto vertex = MyMesh::VertexHandle(i);
            // Get the 3D position of the vertex
            OpenMesh::Vec3f point = mesh.point(vertex);
            // Compute the corresponding grid cell index for the vertex
            GridIndex gridIdx = computeGridIndex(point, gridSize);
            // Add the vertex handle to this thread's local map for the computed grid cell
            privateGridMap[gridIdx].push_back(vertex);
        }

        // Enter critical section to merge the private map into the shared gridMap safely
        #pragma omp critical
        {
            // Iterate over each grid cell in the thread's private map
            for (const auto& [gridIdx, verts] : privateGridMap) {
                // Lock the shared gridMap using a lock_guard for thread safety
                std::lock_guard<std::mutex> lock(gridMapMutex);
                // Merge the vertex handles from the local map into the global gridMap for this grid cell
                gridMap[gridIdx].insert(gridMap[gridIdx].end(), verts.begin(), verts.end());
            }
        }
    }
}

// Function to automatically calculate an optimal grid cell size based on the spatial bounds of the mesh
// it computes the average dimension of the mesh bounding box
float calculateOptimalGridSize(const MyMesh& mesh) {
    // Initialize thread-local min and max bounds
    OpenMesh::Vec3f globalMin(FLT_MAX, FLT_MAX, FLT_MAX);
    OpenMesh::Vec3f globalMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    // Parallel reduction using OpenMP
    #pragma omp parallel
    {
        OpenMesh::Vec3f localMin(FLT_MAX, FLT_MAX, FLT_MAX);
        OpenMesh::Vec3f localMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);

        #pragma omp for nowait
        for (int i = 0; i < mesh.n_vertices(); ++i) {
            OpenMesh::Vec3f point = mesh.point(MyMesh::VertexHandle(i));
            for (int j = 0; j < 3; ++j) {
                localMin[j] = std::min(localMin[j], point[j]);
                localMax[j] = std::max(localMax[j], point[j]);
            }
        }

        // Merge local results into global bounds
        #pragma omp critical
        {
            for (int j = 0; j < 3; ++j) {
                globalMin[j] = std::min(globalMin[j], localMin[j]);
                globalMax[j] = std::max(globalMax[j], localMax[j]);
            }
        }
    }

    // Compute average bounding box size
    float avgSize = (globalMax[0] - globalMin[0] +
                     globalMax[1] - globalMin[1] +
                     globalMax[2] - globalMin[2]) / 3.0f;

    float gridSize = avgSize * 0.2f;
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

        // 정점 수집 및 중심 계산
        for (auto fv_it = original.cfv_iter(face); fv_it.is_valid(); ++fv_it) {
            originalVHs.push_back(*fv_it);
            faceCenter += original.point(*fv_it);
        }

        if (originalVHs.size() != 3) continue; // 삼각형이 아닌 경우 무시
        faceCenter /= 3.0f;

        // 격자 위치 계산
        GridIndex gridIdx = computeGridIndex(faceCenter, gridSize);
        auto& submesh = localSubMeshes[gridIdx];
        auto& vhandleMap = localVertexMaps[gridIdx];

        std::vector<MyMesh::VertexHandle> submeshVHs;
        std::unordered_set<int> uniqueIndices;

        // 정점을 서브메쉬에 매핑
        for (const auto& vh : originalVHs) {
            if (vhandleMap.find(vh) == vhandleMap.end()) {
                vhandleMap[vh] = submesh.add_vertex(original.point(vh));
            }
            auto newVH = vhandleMap[vh];
            submeshVHs.push_back(newVH);
            uniqueIndices.insert(newVH.idx());
        }

        // 중복 정점이 있는 삼각형은 제외
        if (submeshVHs.size() != 3 || uniqueIndices.size() != 3) continue;

        // 면적이 0에 가까운 삼각형도 제외
        auto p0 = submesh.point(submeshVHs[0]);
        auto p1 = submesh.point(submeshVHs[1]);
        auto p2 = submesh.point(submeshVHs[2]);
        float area = OpenMesh::cross(p1 - p0, p2 - p0).norm() * 0.5f;

        if (area > 1e-8f) {
            submesh.add_face(submeshVHs);
        }
    }

    // 결과 정리
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

        OpenMesh::Decimater::ModAspectRatioT<MyMesh>::Handle modAspect;
        if (!decimater.add(modAspect)) {
            std::cerr << "[Error] Failed to add Aspect Ratio module!\n";
            return;
        }
        decimater.module(modAspect).set_aspect_ratio(4.0); // Rejects triangles worse than 1:5 ratio

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





