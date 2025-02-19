#include "MeshUtils.h"

MyMesh extractSubMesh(const MyMesh& mesh, const Grid& grid, int gridIndex) {
    MyMesh submesh;
    int startX = (gridIndex % (grid.dimX * grid.dimY)) % grid.dimX * grid.cellSizeX + grid.minX;
    int startY = (gridIndex % (grid.dimX * grid.dimY)) / grid.dimX * grid.cellSizeY + grid.minY;
    int startZ = gridIndex / (grid.dimX * grid.dimY) * grid.cellSizeZ + grid.minZ;
    int endX = startX + grid.cellSizeX;
    int endY = startY + grid.cellSizeY;
    int endZ = startZ + grid.cellSizeZ;

    // map to store index, vertex
    std::map<MyMesh::VertexHandle, MyMesh::VertexHandle> vhandle_map;

    for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        bool all_vertices_in = true;
        std::vector<MyMesh::VertexHandle> face_vhandles;

        for (auto fv_it = mesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
            auto point = mesh.point(*fv_it);
            if (point[0] >= startX && point[0] < endX &&
                point[1] >= startY && point[1] < endY &&
                point[2] >= startZ && point[2] < endZ) {
                if (vhandle_map.find(*fv_it) == vhandle_map.end()) {
                    // add vertex to submesh
                    MyMesh::VertexHandle new_vh = submesh.add_vertex(mesh.point(*fv_it));
                    vhandle_map[*fv_it] = new_vh;
                }
                face_vhandles.push_back(vhandle_map[*fv_it]);
            } else {
                all_vertices_in = false;
                break;
            }
        }

        if (all_vertices_in && face_vhandles.size() == 3) {
            // if all vertices in grid, add face
            submesh.add_face(face_vhandles);
        }
    }

    return submesh;
}

