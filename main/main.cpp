#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

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

int main() {
    MyMesh mesh;
    
    OpenMesh::IO::read_mesh(mesh, "path_to_your_mesh_file.off");

    Grid grid = {0, 0, 0, 1, 1, 1, 10, 10, 10}; 

    std::vector<int> vertexCellIndices(mesh.n_vertices());

    int idx = 0;
    for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f position = mesh.point(*v_it);
        vertexCellIndices[idx++] = determineGridCell(position, grid);
    }

    for (int cellIndex : vertexCellIndices) {
        std::cout << "Vertex is in cell: " << cellIndex << std::endl;
    }

    return 0;
}