#include "MeshUtils.h"
#include <iostream>

int main() {
    MyMesh mesh;
    std::string filename = "data/bun000_clean.ply";

    //  Load mesh
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    if (!loadMesh(filename, mesh)) {
        return 1;
    }

    // Automatically calculate grid size
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    float gridSize = calculateOptimalGridSize(mesh);

    // Map vertices to grid
    std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> gridMap;
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    mapVerticesToGrid(mesh, gridMap, gridSize);

    // Remove empty grids
    //std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    //std::cin.get();
    //removeEmptyGrids(gridMap, mesh, 5);

    // Extract submeshes
    std::unordered_map<GridIndex, MyMesh> subMeshes;
    std::unordered_map<GridIndex, MyMesh> emptySubMeshes;

    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    extractSubMeshes(mesh, gridMap, subMeshes, emptySubMeshes);

    // Perform decimation
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    decimateSubMeshes(subMeshes);

    // Integrate submeshes (excluding fixed submeshes)
    MyMesh finalMesh;
    std::unordered_map<GridIndex, MyMesh> fixedSubMeshes; // Store fixed submeshes
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    integrateSubMeshes(subMeshes, finalMesh, emptySubMeshes, fixedSubMeshes); // Include fixedSubMeshes

    std::cout << "[Complete] Program finished successfully.\n";
    return 0;
}
