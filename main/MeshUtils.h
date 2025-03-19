#ifndef MESHUTILS_H
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <unordered_map>
#include <string>
#include <unordered_set>

// Mesh type definition
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

// 3D grid index struct
struct GridIndex {
    int x, y, z;
    bool operator==(const GridIndex& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// Hash function definitions (for use with unordered_map)
namespace std {
    template<> struct hash<GridIndex> {
        size_t operator()(const GridIndex& g) const {
            return hash<int>()(g.x) ^ hash<int>()(g.y) ^ hash<int>()(g.z);
        }
    };

    template <>
    // Hash function for OpenMesh::Vec3f
    struct hash<OpenMesh::Vec3f> {
        size_t operator()(const OpenMesh::Vec3f& v) const noexcept {
            size_t h1 = hash<float>()(v[0]);
            size_t h2 = hash<float>()(v[1]);
            size_t h3 = hash<float>()(v[2]);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

// Function to load mesh from file
bool loadMesh(const std::string& filename, MyMesh& mesh);

// Convert vertex position to grid index
GridIndex computeGridIndex(const OpenMesh::Vec3f& point, float gridSize);

// Map mesh vertices to grid cells
void mapVerticesToGrid(const MyMesh& mesh, std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, float gridSize);

// Automatically calculate an optimal grid size based on mesh dimensions
float calculateOptimalGridSize(const MyMesh& mesh);

// Remove grids that have fewer vertices than a given threshold or no connected faces
//void removeEmptyGrids(std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
//                     const MyMesh& mesh, int minThreshold);

// Extract submeshes from original mesh
void extractSubMeshes(const MyMesh& original, 
    const std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>>& gridMap, 
    std::unordered_map<GridIndex, MyMesh>& subMeshes,
    std::unordered_map<GridIndex, MyMesh>& emptySubMeshes); 

// Perform decimation on submeshes to reduce face count and optimize geometry
void decimateSubMeshes(std::unordered_map<GridIndex, MyMesh>& subMeshes);

// Integrate all submeshes into a single final mesh (optionally excluding fixed submeshes)
void integrateSubMeshes(const std::unordered_map<GridIndex, MyMesh>& subMeshes, 
    MyMesh& finalMesh, 
    const std::unordered_map<GridIndex, MyMesh>& emptySubMeshes, 
    const std::unordered_map<GridIndex, MyMesh>& fixedSubMeshes);

// New function declaration for seam-fixing submeshes
//void seamFixSubMeshes(MyMesh& finalMesh, float mergeThreshold);
void processFailedClusters(const std::unordered_map<GridIndex, std::vector<std::vector<MyMesh::VertexHandle>>>& failedFaceClusters);
void mergeCloseVerticesInGrid(MyMesh& mesh, float epsilon = 1e-5f);
void analyzeAndListClusterVertices(const std::unordered_map<GridIndex, std::vector<std::vector<MyMesh::VertexHandle>>>& failedFaceClusters);

#endif // MESHUTILS_H