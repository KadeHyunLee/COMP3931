#ifndef MESHUTILS_H
#define MESHUTILS_H

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

// Hash function definitions
namespace std {
    template<> struct hash<GridIndex> {
        size_t operator()(const GridIndex& g) const {
            return hash<int>()(g.x) ^ hash<int>()(g.y) ^ hash<int>()(g.z);
        }
    };

    template <>
    // Hash function 
    struct hash<std::tuple<int, int, int>> {
        size_t operator()(const std::tuple<int, int, int>& key) const noexcept {
            auto [x, y, z] = key;
            size_t h1 = std::hash<int>()(x);
            size_t h2 = std::hash<int>()(y);
            size_t h3 = std::hash<int>()(z);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

// Automatically calculate an optimal grid size based on mesh dimensions
float calculateOptimalGridSize(const MyMesh& mesh);
float computeDynamicEpsilon(const MyMesh& mesh);

#endif // MESHUTILS_H