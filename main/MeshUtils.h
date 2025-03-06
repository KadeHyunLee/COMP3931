#ifndef MESHUTILS_H
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <unordered_map>
#include <string>

// 메시 타입 정의
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

// 3D 그리드 인덱스 구조체
struct GridIndex {
    int x, y, z;
    bool operator==(const GridIndex& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// 해시 함수 정의 (unordered_map에서 사용하기 위함)
namespace std {
    template<> struct hash<GridIndex> {
        size_t operator()(const GridIndex& g) const {
            return hash<int>()(g.x) ^ hash<int>()(g.y) ^ hash<int>()(g.z);
        }
    };

    template <>
    struct hash<OpenMesh::Vec3f> {
        size_t operator()(const OpenMesh::Vec3f& v) const noexcept {
            size_t h1 = hash<float>()(v[0]);
            size_t h2 = hash<float>()(v[1]);
            size_t h3 = hash<float>()(v[2]);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

// 🟢 메시 불러오기 함수 선언
bool loadMesh(const std::string& filename, MyMesh& mesh);

// 🟢 정점 좌표를 그리드 인덱스로 변환하는 함수 선언
GridIndex computeGridIndex(const OpenMesh::Vec3f& point, float gridSize);

// 🟢 메시의 정점을 그리드에 매핑하는 함수 선언
void mapVerticesToGrid(const MyMesh& mesh, std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, float gridSize);

float calculateOptimalGridSize(const MyMesh& mesh);

void removeEmptyGrids(std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
                      const MyMesh& mesh, int minThreshold);

// 서브메쉬 추출
void extractSubMeshes(const MyMesh& original, 
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes); 

// 🟢 서브메쉬 Decimation 함수 선언
void decimateSubMeshes(std::unordered_map<GridIndex, MyMesh>& subMeshes);

void integrateSubMeshes(const std::unordered_map<GridIndex, MyMesh>& subMeshes, 
    MyMesh& finalMesh, 
    const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes, 
    const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes);

void integrateFixedSubMeshes(const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes, 
    MyMesh& finalMesh_fixed);

#endif // MESHUTILS_H