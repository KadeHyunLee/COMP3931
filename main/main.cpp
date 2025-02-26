#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Core/Geometry/QuadricT.hh>
#include <vector>
#include "MeshUtils.h"
#include <omp.h>
#include <unordered_map>

typedef OpenMesh::Decimater::ModQuadricT<MyMesh> QuadricErrorModule;
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
typedef OpenMesh::Decimater::DecimaterT<MyMesh> MyDecimater;
typedef OpenMesh::Decimater::ModQuadricT<MyMesh>::Handle HModQuadric;

int determineGridCell(const OpenMesh::Vec3f& position, const Grid& grid) {
    int x = std::clamp(static_cast<int>((position[0] - grid.minX) / grid.cellSizeX), 0, grid.dimX - 1);
    int y = std::clamp(static_cast<int>((position[1] - grid.minY) / grid.cellSizeY), 0, grid.dimY - 1);
    int z = std::clamp(static_cast<int>((position[2] - grid.minZ) / grid.cellSizeZ), 0, grid.dimZ - 1);
    
    int gridIndex = x + y * grid.dimX + z * grid.dimX * grid.dimY;

    return gridIndex;
}

void decimateMesh(MyMesh& mesh) {
    MyDecimater decimater(mesh); 
    HModQuadric hModQuadric;
    decimater.add(hModQuadric);
    decimater.module(hModQuadric).set_max_err(0.1);

    decimater.initialize();
    decimater.decimate();
    mesh.garbage_collection();
}

void quantizeVertex(MyMesh::Point& point, float quantizationLevel) {
    
    point[0] = std::floor(point[0] / quantizationLevel + 0.5) * quantizationLevel;
    point[1] = std::floor(point[1] / quantizationLevel + 0.5) * quantizationLevel;
    point[2] = std::floor(point[2] / quantizationLevel + 0.5) * quantizationLevel;
}

void quantizeSubMesh(MyMesh& submesh, float quantizationLevel) {
    
    for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
        auto& point = submesh.point(*v_it);
        quantizeVertex(point, quantizationLevel);
    }
}

const float QUADRIC_ERROR_THRESHOLD = 0.0001f;  // 적절한 임계값 설정 필요

// 쿼드릭 오차 행렬을 저장하는 맵
std::unordered_map<MyMesh::VertexHandle, OpenMesh::Geometry::QuadricT<double>> vertexQuadrics;

// 정점이 기존 정점과 유사한지 판단하는 함수 (쿼드릭 에러 기반)
MyMesh::VertexHandle findNearestVertexByQEM(const OpenMesh::Vec3f& point, 
    MyMesh& mesh,  
    std::unordered_map<MyMesh::VertexHandle, OpenMesh::Geometry::QuadricT<double>>& vertexQuadrics, 
    float threshold) 
{
    MyMesh::VertexHandle bestMatch;
    float minError = FLT_MAX;

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f existingPoint = mesh.point(*v_it);
        float error = (point - existingPoint).sqrnorm();

        if (vertexQuadrics.find(*v_it) != vertexQuadrics.end()) {
            float quadricError = vertexQuadrics.at(*v_it)(point);
            error += quadricError;
        } else {
            std::cerr << "[Warning] Missing quadric error for vertex!\n";
        }

        if (error < minError && error < threshold) {
            bestMatch = *v_it;
            minError = error;
        }
    }

    if (!bestMatch.is_valid()) {
        std::cerr << "[Warning] No valid nearest vertex found for point: " << point << ".\n";
        bestMatch = mesh.add_vertex(point);  // ✅ `const` 제거 후 사용 가능
        vertexQuadrics[bestMatch] = OpenMesh::Geometry::QuadricT<double>(); 
    }

    return bestMatch;
}

// 서브 메시를 병합하는 함수
void integrateSubMesh(MyMesh& original, const MyMesh& submesh, 
                      std::unordered_map<MyMesh::VertexHandle, OpenMesh::Geometry::QuadricT<double>>& vertexQuadrics)
{
    std::unordered_map<MyMesh::VertexHandle, MyMesh::VertexHandle> vhandle_map;

    for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f point = submesh.point(*v_it);
        MyMesh::VertexHandle nearest = findNearestVertexByQEM(point, original, vertexQuadrics, QUADRIC_ERROR_THRESHOLD);

        MyMesh::VertexHandle new_vhandle;
        
        // ✅ 기존에 있는 정점을 사용하거나, 새로 추가하는 부분
        if (nearest.is_valid()) {
            new_vhandle = nearest;
        } else {
            if (vhandle_map.find(*v_it) == vhandle_map.end()) {  // ✅ 추가한 부분
                new_vhandle = original.add_vertex(point);
                vhandle_map[*v_it] = new_vhandle;
                vertexQuadrics[new_vhandle] = vertexQuadrics[*v_it];  // ✅ 쿼드릭 데이터 복사
            } else {
                new_vhandle = vhandle_map[*v_it];
            }
        }

        vhandle_map[*v_it] = new_vhandle;

        if (!new_vhandle.is_valid()) 
        {
            std::cerr << "[Error] Invalid vertex handle detected in integrateSubMesh()!\n";
        }
    }

    // ✅ 면 추가하는 부분
    for (auto f_it = submesh.faces_begin(); f_it != submesh.faces_end(); ++f_it) {
        std::set<MyMesh::VertexHandle> unique_vhandles;
        std::vector<MyMesh::VertexHandle> face_vhandles;

        for (auto fv_it = submesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
            if (vhandle_map.find(*fv_it) == vhandle_map.end()) {
                std::cerr << "[Error] Missing vertex mapping in integrateSubMesh()! Vertex index: " 
                        << *fv_it << std::endl;
                continue;
            }
            unique_vhandles.insert(vhandle_map[*fv_it]); // 중복 방지
        }

        if (unique_vhandles.size() == 3) {
            face_vhandles.assign(unique_vhandles.begin(), unique_vhandles.end());
            MyMesh::FaceHandle newFace = original.add_face(face_vhandles);
            if (!newFace.is_valid()) {
                std::cerr << "[Error] Failed to add face in integrateSubMesh().\n";
            } else {
                std::cout << "[Debug] Successfully added face in integrateSubMesh().\n";
            }
        } else {
            std::cerr << "[Error] Face skipped due to incorrect vertex count: " 
                    << unique_vhandles.size() << std::endl;
        }
    }
}

int main() {
    MyMesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, "data/bun000_clean.ply")) {
        std::cerr << "Error: Unable to read mesh from 'path_to_mesh_file'. Please check the file path and permissions." << std::endl;
        return 1; // Return an error code
    }

    OpenMesh::Vec3f minBound(FLT_MAX, FLT_MAX, FLT_MAX);
    OpenMesh::Vec3f maxBound(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f position = mesh.point(*v_it);
        minBound.minimize(position);
        maxBound.maximize(position);
    }

    Grid grid;
    grid.minX = minBound[0];
    grid.minY = minBound[1];
    grid.minZ = minBound[2];
    grid.cellSizeX = (maxBound[0] - minBound[0]) / 6;  
    grid.cellSizeY = (maxBound[1] - minBound[1]) / 6;
    grid.cellSizeZ = (maxBound[2] - minBound[2]) / 6;
    grid.dimX = 6;
    grid.dimY = 6;
    grid.dimZ = 6; 

    float quantizationLevel = 0.1; 
    std::vector<int> vertexCellIndices(mesh.n_vertices());

    int idx = 0;
    for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) 
    {
        OpenMesh::Vec3f position = mesh.point(*v_it);
        int gridIndex = determineGridCell(position, grid);
        vertexCellIndices[idx++] = gridIndex;
    }

    MyMesh finalMesh;
    std::unordered_map<MyMesh::VertexHandle, OpenMesh::Geometry::QuadricT<double>> vertexQuadrics;

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < grid.dimX * grid.dimY * grid.dimZ; ++i) 
        {
            MyMesh submesh = extractSubMesh(mesh, grid, i);
            std::cout << "[Debug] Before decimation: " << submesh.n_vertices() << " vertices\n";
            decimateMesh(submesh);
            std::cout << "[Debug] After decimation: " << submesh.n_vertices() << " vertices\n";
        }

        #pragma omp single
        {
            for (int i = 0; i < grid.dimX * grid.dimY * grid.dimZ; ++i) 
            {
                integrateSubMesh(finalMesh, extractSubMesh(mesh, grid, i), vertexQuadrics);
            }
        }
    }
    std::cout << "[Debug] Final mesh has " << finalMesh.n_faces() << " faces." << std::endl;

    if (!OpenMesh::IO::write_mesh(finalMesh, "output/final_mesh.ply")) {
        std::cerr << "Error: Unable to write final mesh to 'output/final_mesh.ply'. Check file path and permissions." << std::endl;
        return 1;
    }

    std::cout << "Final mesh saved as 'output/final_mesh.ply' with " << finalMesh.n_vertices() << " vertices.\n";
    
    return 0;
}