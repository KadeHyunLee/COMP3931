#include "MeshUtils.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

// 🟢 메시 파일을 로드하는 함수 구현
bool loadMesh(const std::string& filename, MyMesh& mesh) {
    if (!OpenMesh::IO::read_mesh(mesh, filename)) {
        std::cerr << "[Error] Failed to load mesh: " << filename << "\n";
        return false;
    }
    std::cout << "[Debug] Loaded mesh: " << filename << " with "
              << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";
    return true;
}

// 🟢 그리드 인덱스 계산
GridIndex computeGridIndex(const OpenMesh::Vec3f& point, float gridSize) {
    return { 
        static_cast<int>(point[0] / gridSize), 
        static_cast<int>(point[1] / gridSize), 
        static_cast<int>(point[2] / gridSize) 
    };
}

// 🟢 정점을 그리드에 매핑하는 함수
void mapVerticesToGrid(const MyMesh& mesh, std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, float gridSize) {
    std::cout << "[Debug] Mapping vertices to grid...\n";

    int minIdxX = INT_MAX, minIdxY = INT_MAX, minIdxZ = INT_MAX;
    int maxIdxX = INT_MIN, maxIdxY = INT_MIN, maxIdxZ = INT_MIN;

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f point = mesh.point(*v_it);
        GridIndex gridIdx = computeGridIndex(point, gridSize);

        // 🔹 빈 그리드가 생기는 것을 방지하기 위해 체크
        if (gridMap[gridIdx].empty()) {
            gridMap[gridIdx] = std::vector<MyMesh::VertexHandle>();
        }

        gridMap[gridIdx].push_back(*v_it);

        // 그리드 인덱스 범위 추적
        minIdxX = std::min(minIdxX, gridIdx.x);
        minIdxY = std::min(minIdxY, gridIdx.y);
        minIdxZ = std::min(minIdxZ, gridIdx.z);

        maxIdxX = std::max(maxIdxX, gridIdx.x);
        maxIdxY = std::max(maxIdxY, gridIdx.y);
        maxIdxZ = std::max(maxIdxZ, gridIdx.z);
    }

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

float calculateOptimalGridSize(const MyMesh& mesh) {
    OpenMesh::Vec3f minBounds(FLT_MAX, FLT_MAX, FLT_MAX);
    OpenMesh::Vec3f maxBounds(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f point = mesh.point(*v_it);
        for (int i = 0; i < 3; i++) {
            minBounds[i] = std::min(minBounds[i], point[i]);
            maxBounds[i] = std::max(maxBounds[i], point[i]);
        }
    }

    float avgSize = (maxBounds[0] - minBounds[0] +
                     maxBounds[1] - minBounds[1] +
                     maxBounds[2] - minBounds[2]) / 3.0f;

    float gridSize = avgSize * 0.05f;  // 메시 크기의 5%를 셀 크기로 설정
    std::cout << "[Debug] Auto-Calculated Grid Size: " << gridSize << "\n";

    return gridSize;
}

// 정점 - 지금 경계에 정점 5 미만인 그리드 존재, 나중에 따로 처리하는 로직 구현해야함

void removeEmptyGrids(std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
                      const MyMesh& mesh, int minThreshold) {
    std::cout << "[Debug] Before Removing: " << gridMap.size() << " grids\n";
    
    auto it = gridMap.begin();
    while (it != gridMap.end()) {
        const auto& vertices = it->second;

        if (vertices.size() < minThreshold) {
            bool hasConnectedFaces = false;

            // 🔹 현재 그리드 정점이 연결된 삼각형 검사
            for (const auto& v : vertices) {
                if (mesh.valence(v) > 0) { 
                    hasConnectedFaces = true;
                    break;
                }
            }

            // 🔹 연결된 삼각형이 없으면 정점과 함께 삭제
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

void extractSubMeshes(const MyMesh& original, 
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes) { // ✅ 추가된 emptySubMeshes
    std::cout << "[Debug] Extracting submeshes from grids...\n";

    int subMeshCount = 0; 
    int emptySubMeshCount = 0;

    for (const auto& cell : gridMap) {
        const GridIndex& gridIdx = cell.first;
        const auto& vertices = cell.second;

        MyMesh submesh;
        std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> vhandleMap;

        // 정점 추가
        for (const auto& v : vertices) {
            if (vhandleMap.find(v) == vhandleMap.end()) {
                OpenMesh::Vec3f point = original.point(v);
                MyMesh::VertexHandle new_vhandle = submesh.add_vertex(point);
                vhandleMap[v] = new_vhandle;
            }
        }

        int faceAddCount = 0;
        for (const auto& v : vertices) {
            for (MyMesh::ConstVertexFaceIter vf_it = original.cvf_iter(v); vf_it.is_valid(); ++vf_it) {
                MyMesh::FaceHandle face = *vf_it;
                std::vector<MyMesh::VertexHandle> face_vhandles;

                bool validFace = true;
                for (MyMesh::ConstFaceVertexIter fv_it = original.cfv_iter(face); fv_it.is_valid(); ++fv_it) {
                    if (vhandleMap.find(*fv_it) != vhandleMap.end()) {
                        face_vhandles.push_back(vhandleMap[*fv_it]);
                    } else {
                        validFace = false;
                        break;
                    }
                }

                if (validFace && face_vhandles.size() == 3) {
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

        if (submesh.n_vertices() > 0) {
            subMeshes[gridIdx] = submesh;
            subMeshCount++;

            if (submesh.n_faces() == 0) {
                std::cerr << "[Warning] Submesh for Grid (" << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z
                          << ") has 0 faces but " << submesh.n_vertices() << " vertices!\n";
                emptySubMeshes[gridIdx] = submesh; // ✅ 빈 서브메쉬 저장
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

    // 🔹 최종 병합된 메시 저장 (PLY 파일)
    OpenMesh::IO::write_mesh(finalMesh, "final_integrated_mesh.ply");
    std::cout << "[Saved] Final Integrated Mesh saved as: final_integrated_mesh.ply\n";
}

void integrateFixedSubMeshes(const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes, 
    MyMesh& finalMesh_fixed) {
    std::cout << "[Debug] Integrating Fixed Submeshes...\n";

    std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> globalVertexMap;
    int mergedFaces = 0, mergedVertices = 0;

    for (const auto& [gridIdx, submesh] : fixedSubMeshes) {
        for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
            OpenMesh::Vec3f point = submesh.point(*v_it);
            MyMesh::VertexHandle new_vhandle = finalMesh_fixed.add_vertex(point);
            globalVertexMap[*v_it] = new_vhandle;
            mergedVertices++;
        }

        for (auto f_it = submesh.faces_begin(); f_it != submesh.faces_end(); ++f_it) {
            std::vector<MyMesh::VertexHandle> face_vhandles;
            for (auto fv_it = submesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
                face_vhandles.push_back(globalVertexMap[*fv_it]);
            }

            if (face_vhandles.size() == 3) {
                finalMesh_fixed.add_face(face_vhandles);
                mergedFaces++;
            } else {
                std::cerr << "[Warning] Skipping invalid face from Grid (" 
                          << gridIdx.x << ", " << gridIdx.y << ", " << gridIdx.z 
                          << ") with " << face_vhandles.size() << " vertices.\n";
            }
        }
    }

    std::cout << "[Debug] Integration of Fixed Submeshes Completed!\n";
    std::cout << "[Summary] Fixed Submesh Integration Result: " << finalMesh_fixed.n_vertices() 
              << " vertices, " << finalMesh_fixed.n_faces() << " faces.\n";

    // 🔹 최종 보정된 메시 저장 (PLY 파일)
    OpenMesh::IO::write_mesh(finalMesh_fixed, "fixed_integrated_mesh.ply");
    std::cout << "[Saved] Fixed Integrated Mesh saved as: fixed_integrated_mesh.ply\n";
}