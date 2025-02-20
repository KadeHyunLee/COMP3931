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

    for (auto face_iter = mesh.faces_begin(); face_iter != mesh.faces_end(); ++face_iter) {
        bool allVerticesIn = true;
        std::vector<MyMesh::VertexHandle> face_vhandles;

        for (auto face_vertex_iter = mesh.cfv_iter(*face_iter); fv_it.is_valid(); ++face_vertex_iter) {
            auto point = mesh.point(*face_vertex_iter);
            if (point[0] >= startX && point[0] < endX &&
                point[1] >= startY && point[1] < endY &&
                point[2] >= startZ && point[2] < endZ) {
                if (vhandle_map.find(*face_vertex_iter) == vhandle_map.end()) {
                    // add vertex to submesh
                    MyMesh::VertexHandle new_vhandles = submesh.add_vertex(mesh.point(*face_vertex_iter));
                    vhandle_map[*face_vertex_iter] = new_vhandles;
                }
                face_vhandles.push_back(vhandle_map[*face_vertex_iter]);
            } else {
                allVerticesIn = false;
                break;
            }
        }

        if (allVerticesIn && face_vhandles.size() == 3) {
            // if all vertices in grid, add face
            submesh.add_face(face_vhandles);
        }
    }

    return submesh;
}

