#include "MeshUtils.h"
#include <iostream>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char ** argv) {
    std::string filename;
    float decimationRatio = 0.5f;
    int threadCount = omp_get_max_threads(); // Default thread count

    if (argc >= 2) {
        filename = std::string("data/") + argv[1];
    } else {
        std::cout << "[Input] Enter mesh filename (inside data/): ";
        std::string userInput;
        std::cin >> userInput;
        filename = std::string("data/") + userInput;
    }

    if (argc >= 3) {
        decimationRatio = std::stof(argv[2]);
    } else {
        std::cout << "[Input] Enter decimation ratio (e.g., 0.3): ";
        std::cin >> decimationRatio;
    }
    
    if (argc >= 4) {
        threadCount = std::stoi(argv[3]);
    } else {
        std::cout << "[Input] Enter number of threads to use: ";
        std::cin >> threadCount;
    }
    
    // Try loading the mesh, allow retry if failed
    MyMesh mesh;
    while (!loadMesh(filename, mesh)) {
        std::cout << "[Retry] Enter mesh filename (inside data/): ";
        std::string userInput;
        std::cin >> userInput;
        filename = std::string("data/") + userInput;
    }

    //  Load mesh
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    if (!loadMesh(filename, mesh)) {
        return 1;
    }

    // Automatically calculate grid size
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    auto start = std::chrono::high_resolution_clock::now();
    #ifdef _OPENMP
    omp_set_num_threads(threadCount);
    std::cout << "[Debug] OpenMP Thread Count Set To: " << threadCount << "\n";
    std::cout << "[Debug] OpenMP Max Threads: " << omp_get_max_threads() << "\n";
    #else
    std::cout << "[Debug] OpenMP not enabled.\n";
    #endif
    float gridSize = calculateOptimalGridSize(mesh);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "[Timing] Grid Size Calculation took " << duration << " ms\n";

    // Map vertices to grid
    std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> gridMap;
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();

    auto mapStart = std::chrono::high_resolution_clock::now();
    mapVerticesToGrid(mesh, gridMap, gridSize);
    auto mapEnd = std::chrono::high_resolution_clock::now();
    auto mapDuration = std::chrono::duration_cast<std::chrono::milliseconds>(mapEnd - mapStart).count();
    std::cout << "[Timing] Vertex-to-Grid Mapping took " << mapDuration << " ms\n";

    // Extract submeshes
    std::unordered_map<GridIndex, MyMesh> subMeshes;
    std::unordered_map<GridIndex, MyMesh> emptySubMeshes;

    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    auto extractStart = std::chrono::high_resolution_clock::now();
    extractSubMeshes(mesh, gridMap, subMeshes, emptySubMeshes);
    auto extractEnd = std::chrono::high_resolution_clock::now();
    auto extractDuration = std::chrono::duration_cast<std::chrono::milliseconds>(extractEnd - extractStart).count();
    std::cout << "[Timing] Submesh Extraction took " << extractDuration << " ms\n";
    std::cout << "[Debug] Total Extracted Submeshes: " << subMeshes.size() << "\n";

    // Perform decimation
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    auto decimateStart = std::chrono::high_resolution_clock::now();
    decimateSubMeshes(subMeshes, decimationRatio);
    auto decimateEnd = std::chrono::high_resolution_clock::now();
    auto decimateDuration = std::chrono::duration_cast<std::chrono::milliseconds>(decimateEnd - decimateStart).count();
    std::cout << "[Timing] Decimation took " << decimateDuration << " ms\n";

    // Integrate submeshes 
    MyMesh finalMesh;
    std::unordered_map<GridIndex, MyMesh> fixedSubMeshes; // Store fixed submeshes
    std::cout << "\n[Press SPACE and ENTER to continue...]\n";
    std::cin.get();
    integrateSubMeshes(subMeshes, finalMesh, emptySubMeshes, fixedSubMeshes); // Include fixedSubMeshes

    std::cout << "[Complete] Program finished successfully.\n";
    return 0;
}
