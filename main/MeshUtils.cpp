#include "MeshUtils.h"

MyMesh extractSubMesh(const MyMesh& mesh, const Grid& grid, int gridIndex) {
    MyMesh submesh;

    int gridX = gridIndex % grid.dimX;
    int gridY = (gridIndex / grid.dimX) % grid.dimY;
    int gridZ = (gridIndex / (grid.dimX * grid.dimY)) % grid.dimZ;  // ✅ 수정

    float startX = grid.minX + gridX * grid.cellSizeX;
    float startY = grid.minY + gridY * grid.cellSizeY;
    float startZ = grid.minZ + gridZ * grid.cellSizeZ;

    float endX = grid.minX + (gridX + 1) * grid.cellSizeX;  // ✅ 수정
    float endY = grid.minY + (gridY + 1) * grid.cellSizeY;  // ✅ 수정
    float endZ = grid.minZ + (gridZ + 1) * grid.cellSizeZ;  // ✅ 수정

    // map to store index, vertex
    std::map<MyMesh::VertexHandle, MyMesh::VertexHandle> vhandle_map;

    int vertexCount = 0;
    int faceCount = 0;

    for (auto face_iter = mesh.faces_begin(); face_iter != mesh.faces_end(); ++face_iter) {
        bool allVerticesIn = true;
        std::vector<MyMesh::VertexHandle> face_vhandles;

        int vertexInGrid = 0; 

        for (auto face_vertex_iter = mesh.cfv_iter(*face_iter); face_vertex_iter.is_valid(); ++face_vertex_iter) {
            auto point = mesh.point(*face_vertex_iter);
            bool inside = (point[0] >= startX && point[0] <= endX &&
                        point[1] >= startY && point[1] <= endY &&
                        point[2] >= startZ && point[2] <= endZ);

            if (inside) {
                vertexInGrid++;
                if (vhandle_map.find(*face_vertex_iter) == vhandle_map.end()) 
                {
                    // add vertex to submesh
                    MyMesh::VertexHandle new_vhandles = submesh.add_vertex(mesh.point(*face_vertex_iter));
                    vhandle_map[*face_vertex_iter] = new_vhandles;
                    vertexCount++; // 디버깅용
                }
                face_vhandles.push_back(vhandle_map[*face_vertex_iter]);
            } else {
                allVerticesIn = false;
                break;
            }
        }

        if (vertexInGrid > 0 && face_vhandles.size() == 3) 
        {
            MyMesh::FaceHandle newFace = submesh.add_face(face_vhandles);
            if (!newFace.is_valid()) {
                std::cerr << "[Error] Failed to add face in extractSubMesh()! Vertex count: " 
                        << face_vhandles.size() << std::endl;
            } else {
                std::cout << "[Debug] Added face successfully in extractSubMesh().\n";
            }
        }
    }
    return submesh;
}

