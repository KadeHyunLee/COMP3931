#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <vector>

#include <omp.h>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
typedef OpenMesh::Decimater::DecimaterT<MyMesh> MyDecimater;
typedef OpenMesh::Decimater::ModQuadricT<MyMesh>::Handle HModQuadric;

struct Grid {
    float minX, minY, minZ;
    float cellSizeX, cellSizeY, cellSizeZ;
    int dimX, dimY, dimZ;
};

int determineGridCell(const OpenMesh::Vec3f& position, const Grid& grid) {
    int x = static_cast<int>((position[0] - grid.minX) / grid.cellSizeX);
    int y = static_cast<int>((position[1] - grid.minY) / grid.cellSizeY);
    int z = static_cast<int>((position[2] - grid.minZ) / grid.cellSizeZ);
    return x + y * grid.dimX + z * grid.dimX * grid.dimY;
}

void decimateMesh(MyMesh& mesh) {
    MyDecimater decimater(mesh); 
    HModQuadric hModQuadric;
    decimater.add(hModQuadric);

    decimater.module(hModQuadric).set_max_err(0.005); 

    decimater.initialize();
    decimater.decimate();
    mesh.garbage_collection();
}

void quantizeVertex(MyMesh::Point& point, float quantizationLevel) {
    // 각 좌표를 퀀타이제이션 레벨에 따라 반올림합니다.
    point[0] = std::floor(point[0] / quantizationLevel + 0.5) * quantizationLevel;
    point[1] = std::floor(point[1] / quantizationLevel + 0.5) * quantizationLevel;
    point[2] = std::floor(point[2] / quantizationLevel + 0.5) * quantizationLevel;
}

void quantizeSubMesh(MyMesh& submesh, float quantizationLevel) {
    // 모든 정점을 순회하면서 각 정점의 좌표를 퀀타이제이션합니다.
    for (auto v_it = submesh.vertices_begin(); v_it != submesh.vertices_end(); ++v_it) {
        auto& point = submesh.point(*v_it);
        quantizeVertex(point, quantizationLevel);
    }
}

int main() {
    MyMesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, "path_to_mesh_file")) {
        std::cerr << "Error: Unable to read mesh from 'path_to_mesh_file'. Please check the file path and permissions." << std::endl;
        return 1; // Return an error code
    }

    Grid grid = {0, 0, 0, 1, 1, 1, 10, 10, 10};
    float quantizationLevel = 0.1; 
    std::vector<int> vertexCellIndices(mesh.n_vertices());

    int idx = 0;
    for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f position = mesh.point(*v_it);
        vertexCellIndices[idx++] = determineGridCell(position, grid);
    }

    #pragma omp parallel for
    for (int i = 0; i < grid.dimX * grid.dimY * grid.dimZ; ++i) {
        MyMesh submesh = extractSubMesh(mesh, grid, i); 
        decimateMesh(submesh);
        quantizeSubMesh(submesh, quantizationLevel);

        // integrate
        
        
    }

    for (int cellIndex : vertexCellIndices) {
        std::cout << "Vertex is in cell: " << cellIndex << std::endl;
    }

    return 0;
}