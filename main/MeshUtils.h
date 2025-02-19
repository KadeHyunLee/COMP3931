#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

MyMesh extractSubMesh(const MyMesh& mesh, const Grid& grid, int gridIndex);