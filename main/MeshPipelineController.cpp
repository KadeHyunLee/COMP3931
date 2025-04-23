#include "MeshPipelineController.h"
#include "MeshLoader.h"
#include "MeshUtils.h"
#include "GridPartitioner.h"
#include "SubmeshExtractor.h"
#include "SubmeshDecimator.h"
#include "MeshIntegrator.h"

#include <iostream>
#include <chrono>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

void MeshPipelineController::runPipeline() {
    
    std::string filename;
    float decimationRatio;
    
    std::cout << "[Input] Enter mesh filename (inside data/): ";
    std::string userInput;
    std::cin >> userInput;
    filename = std::string("data/") + userInput;

    std::cout << "[Input] Enter decimation ratio (e.g., 0.3): ";
    std::cin >> decimationRatio;

    #ifdef _OPENMP
    int threadCount = omp_get_max_threads();
    omp_set_num_threads(threadCount);
    std::cout << "[Debug] Auto thread count set to: " << threadCount << "\n";
    #endif
    
    MyMesh mesh;
    if (!MeshLoader::load(filename, mesh)) {
        std::cerr << "[Error] Failed to load mesh from: " << filename << "\n";
        return;
    }

    #ifdef _OPENMP
    omp_set_num_threads(threadCount);
    std::cout << "[Debug] OpenMP Thread Count Set To: " << threadCount << "\n";
    #endif

    float gridSize = calculateOptimalGridSize(mesh);
    
    GridPartitioner partitioner(gridSize);

    std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> gridMap;
    partitioner.mapVertices(mesh, gridMap);

    std::unordered_map<GridIndex, MyMesh> subMeshes, emptySubMeshes;
    SubmeshExtractor::extract(mesh, gridMap, subMeshes, emptySubMeshes);

    SubmeshDecimator::decimate(subMeshes, decimationRatio);

    MyMesh finalMesh;
    std::unordered_map<GridIndex, MyMesh> fixedSubMeshes;
    MeshIntegrator::integrate(subMeshes, finalMesh, emptySubMeshes, fixedSubMeshes);

    std::cout << "[Complete] Pipeline finished successfully.\n";
}

void MeshPipelineController::runPipelineInTesting() {

    std::string filename;
    float decimationRatio;
    int threadCount;
    
    std::cout << "[Input] Enter mesh filename (inside data/): ";
    std::string userInput;
    std::cin >> userInput;
    filename = std::string("data/") + userInput;
    
    std::cout << "[Input] Enter decimation ratio (e.g., 0.3): ";
    std::cin >> decimationRatio;
    
    std::cout << "[Input] Enter number of threads to use: ";
    std::cin >> threadCount;
    
    #ifdef _OPENMP
        omp_set_num_threads(threadCount);
        std::cout << "[Debug] OpenMP Thread Count Set To: " << threadCount << "\n";
        std::cout << "[Debug] Max Threads Available: " << omp_get_max_threads() << "\n";
    #endif

    MyMesh mesh;
    if (!MeshLoader::load(filename, mesh)) {
        std::cerr << "[Error] Failed to load mesh.\n";
        return;
    }

    #ifdef _OPENMP
    omp_set_num_threads(threadCount);
    std::cout << "[Debug] OpenMP Thread Count: " << threadCount << "\n";
    std::cout << "[Debug] Max Threads Available: " << omp_get_max_threads() << "\n";
    #endif

    auto start = std::chrono::high_resolution_clock::now();
    float gridSize = calculateOptimalGridSize(mesh);
    GridPartitioner partitioner(gridSize);

    std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> gridMap;
    start = std::chrono::high_resolution_clock::now();
    partitioner.mapVertices(mesh, gridMap);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "[Timing] Vertex Mapping took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    std::cout << "[Timing] Grid Size: " << gridSize << ", Time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    std::unordered_map<GridIndex, MyMesh> submeshes, emptySubmeshes;
    start = std::chrono::high_resolution_clock::now();
    SubmeshExtractor::extract(mesh, gridMap, submeshes, emptySubmeshes);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[Timing] Submesh Extraction took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    SubmeshDecimator::decimate(submeshes, decimationRatio);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[Timing] Decimation took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    MyMesh finalMesh;
    std::unordered_map<GridIndex, MyMesh> fixedSubMeshes;
    start = std::chrono::high_resolution_clock::now();
    MeshIntegrator::integrate(submeshes, finalMesh, emptySubmeshes, fixedSubMeshes);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[Timing] Integration took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    std::cout << "[Complete] Testing pipeline finished.\n";
}