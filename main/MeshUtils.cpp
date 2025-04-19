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

    // Mutex map for per-grid locking
    std::mutex gridMapMutex;
    
    // Pre-initialize mutexes to avoid concurrent insertion
    for (int i = 0; i < mesh.n_vertices(); ++i) {
        auto vertex = MyMesh::VertexHandle(i);
        OpenMesh::Vec3f point = mesh.point(vertex);
        GridIndex gridIdx = computeGridIndex(point, gridSize);
    }
    
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
        #pragma omp critical
        {
            // Merge the private map into the shared gridMap with locking
            for (const auto& [gridIdx, verts] : privateGridMap) {
                std::lock_guard<std::mutex> lock(gridMapMutex);
                gridMap[gridIdx].insert(gridMap[gridIdx].end(), verts.begin(), verts.end());
            }
        }
    }
}

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

void extractSubMeshes(const MyMesh& original,
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap,
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes) {
    
    std::unordered_map<GridIndex, std::mutex> submeshLocks;
    for (const auto& [gridIdx, _] : gridMap) {
        submeshLocks[gridIdx];  // Ensure mutex exists for each grid index
    }

    std::cout << "[Debug] Extracting submeshes (face-centric with OpenMP)...\n";

    float gridSize = calculateOptimalGridSize(original);

    std::vector<MyMesh::FaceHandle> allFaces;
    for (auto f_it = original.faces_begin(); f_it != original.faces_end(); ++f_it) {
        allFaces.push_back(*f_it);
    }

    std::unordered_map<GridIndex, MyMesh> localSubMeshes;
    std::unordered_map<GridIndex, std::unordered_map<std::tuple<int, int, int>, MyMesh::VertexHandle>> globalVertexMaps;

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
            MyMesh::FaceHandle face = allFaces[i];
            std::vector<MyMesh::VertexHandle> originalVHs;
            OpenMesh::Vec3f faceCenter(0, 0, 0);

            for (auto fv_it = original.cfv_iter(face); fv_it.is_valid(); ++fv_it) {
                originalVHs.push_back(*fv_it);
                faceCenter += original.point(*fv_it);
            }

            if (originalVHs.size() != 3) continue;
            faceCenter /= 3.0f;

            GridIndex gridIdx = computeGridIndex(faceCenter, gridSize);
            auto& submesh = threadSubMeshes[gridIdx];
            auto& vhandleMap = threadVertexMaps[gridIdx];

            std::vector<MyMesh::VertexHandle> submeshVHs;
            std::unordered_set<int> uniqueIndices;

            for (const auto& vh : originalVHs) {
                OpenMesh::Vec3f point = original.point(vh);
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

        for (auto& [gridIdx, mesh] : threadSubMeshes) {
            std::lock_guard<std::mutex> lock(submeshLocks[gridIdx]);
            if (localSubMeshes.find(gridIdx) == localSubMeshes.end()) {
                localSubMeshes[gridIdx] = mesh;
                globalVertexMaps[gridIdx] = threadVertexMaps[gridIdx];
            } else {
                auto& dst = localSubMeshes[gridIdx];
                auto& globalMap = globalVertexMaps[gridIdx];

                std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> handleMap;
                for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
                    OpenMesh::Vec3f point = mesh.point(*v_it);
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

void decimateSubMeshes(std::unordered_map<GridIndex, MyMesh>& subMeshes, float decimationRatio) {
    //std::cout << "[Debug] Number of Submeshes: " << subMeshes.size() << "\n";
    //std::cout << "[Debug] Performing Decimation on Submeshes...\n";

    int totalFacesBefore = 0;
    int totalVerticesBefore = 0;
    int totalFacesAfter = 0;
    int totalVerticesAfter = 0;

    // Convert keys to indexable vector for parallel iteration
    std::vector<GridIndex> gridIndices;
    for (const auto& [gridIdx, _] : subMeshes) {
        gridIndices.push_back(gridIdx);
    }

    #pragma omp parallel for reduction(+:totalFacesBefore, totalVerticesBefore, totalFacesAfter, totalVerticesAfter)
    for (int i = 0; i < static_cast<int>(gridIndices.size()); ++i) {
        GridIndex gridIdx = gridIndices[i];
        auto& submesh = subMeshes[gridIdx];

        //std::cout << "[Debug] Processing Submesh at Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z << ")\n";
        //std::cout << " - Vertices: " << submesh.n_vertices() << ", Faces: " << submesh.n_faces() << "\n";

        int facesBefore = submesh.n_faces();
        if (facesBefore <= 3 || facesBefore == 0) {
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
        // Optional: Set aspect ratio constraint if needed
        // decimater.module(modAspect).set_aspect_ratio(4.0);
 
        // Dynamically adjust maximum error tolerance based on decimation ratio
        float dynamicErr = (1.0f - decimationRatio) * 0.15f;
        int target_faces = std::max(static_cast<int>(facesBefore * decimationRatio), 3);
        decimater.module(modQuadric).set_max_err(dynamicErr);

        decimater.initialize();

        decimater.decimate_to_faces(target_faces);
        submesh.garbage_collection();

        int facesAfter = submesh.n_faces();
        int verticesAfter = submesh.n_vertices();

        totalFacesAfter += facesAfter;
        totalVerticesAfter += verticesAfter;
    }

    std::cout << "[Summary] Decimation Completed!\n";
    std::cout << " - Total Faces Before: " << totalFacesBefore << "\n";
    std::cout << " - Total Faces After: " << totalFacesAfter << "\n";
    std::cout << " - Total Vertices Before: " << totalVerticesBefore << "\n";
    std::cout << " - Total Vertices After: " << totalVerticesAfter << "\n";

    float faceReductionRate = 100.0f * (1.0f - static_cast<float>(totalFacesAfter) / totalFacesBefore);
    float vertexReductionRate = 100.0f * (1.0f - static_cast<float>(totalVerticesAfter) / totalVerticesBefore);
    
    std::cout << " - Face Reduction Rate: " << faceReductionRate << "%\n";
    std::cout << " - Vertex Reduction Rate: " << vertexReductionRate << "%\n";
    std::cout << " - Target Decimation Ratio: " << decimationRatio << "\n";
}


void integrateSubMeshes(const std::unordered_map<GridIndex, MyMesh>& subMeshes, 
    MyMesh& finalMesh, 
    const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes, 
    const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes) {

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

        std::vector<std::pair<OpenMesh::Vec3f, MyMesh::VertexHandle>> newVertices;
        std::vector<std::vector<MyMesh::VertexHandle>> newFaces;

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
                newVertices.emplace_back(point, newVH);
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
                {
                    std::lock_guard<std::mutex> lock(finalMeshMutex);
                    finalMesh.add_face(vhs);
                }
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




