#include "MeshLoader.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>

bool MeshLoader::load(const std::string& filename, MyMesh& mesh) {
    if (!OpenMesh::IO::read_mesh(mesh, filename)) {
        std::cerr << "[Error] Failed to load mesh: " << filename << "\n";
        std::cerr << "[Debug] Please check the file path and format. Supported formats include: .obj, .ply, .stl\n";
        return false;
    }

    std::cout << "[Debug] Loaded mesh: " << filename << " with "
              << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";
    return true;
}