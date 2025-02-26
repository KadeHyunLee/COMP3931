#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <vector>
#include <map>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

struct Grid {
    float minX, minY, minZ;
    float maxX, maxY, maxZ;
    float cellSizeX, cellSizeY, cellSizeZ;
    int dimX, dimY, dimZ;
};

MyMesh extractSubMesh(const MyMesh& mesh, const Grid& grid, int gridIndex);

MyMesh integrateMeshes(const std::vector<MyMesh>& submeshes, const Grid& grid);