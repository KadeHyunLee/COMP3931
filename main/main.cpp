#include "MeshUtils.h"
#include <iostream>

int main() {
    MyMesh mesh;
    std::string filename = "data/bun000_clean.ply";

    // 🟢 메시 로드
    if (!loadMesh(filename, mesh)) {
        return 1;
    }

    // 🟢 자동으로 그리드 크기 계산
    float gridSize = calculateOptimalGridSize(mesh);

    // 🟢 그리드에 정점 매핑
    std::unordered_map<GridIndex, std::vector<MyMesh::VertexHandle>> gridMap;
    mapVerticesToGrid(mesh, gridMap, gridSize);

    // 🟢 빈 그리드 제거
    removeEmptyGrids(gridMap, mesh, 5);

    // 🟢 Submesh 분할 수행
    std::unordered_map<GridIndex, MyMesh> subMeshes;
    std::unordered_map<GridIndex, MyMesh> emptySubMeshes;

    extractSubMeshes(mesh, gridMap, subMeshes, emptySubMeshes);


    // 🟢 Decimation 수행
    decimateSubMeshes(subMeshes);

    // 🟢 서브메쉬 통합 수행 (보정된 서브메쉬 제외)
    MyMesh finalMesh;
    std::unordered_map<GridIndex, MyMesh> fixedSubMeshes; // 🔹 보정된 서브메쉬 저장
    integrateSubMeshes(subMeshes, finalMesh, emptySubMeshes, fixedSubMeshes); // ✅ fixedSubMeshes 추가

    // 🟢 보정된 서브메쉬 통합 (별도 파일로 저장)
    MyMesh finalMesh_fixed;
    integrateFixedSubMeshes(fixedSubMeshes, finalMesh_fixed);

    return 0;
}
