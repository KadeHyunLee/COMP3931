#ifndef MESHPIPELINECONTROLLER_H
#define MESHPIPELINECONTROLLER_H

#include "MeshUtils.h"
#include <string>
#include <unordered_map>

class MeshPipelineController {
public:
    void runPipeline();
    void runPipelineInTesting();

    //void setImproveQuality(int flag) { improveQuality_ = flag; }
    /*
    private:
        void improveBoundaryIfNeeded(std::unordered_map<GridIndex, MyMesh>& subMeshes);
        int improveQuality_ = 0;
    */

};

#endif // MESHPIPELINECONTROLLER_H
