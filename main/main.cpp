#include "MeshUtils.h"
#include "MeshPipelineController.h"
#include <iostream>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char ** argv) {
    int mode = 0;
    std::cout << "[Input] Choose mode: (0 = Normal, 1 = Testing for execution time): ";
    std::cin >> mode;

    MeshPipelineController controller;
    
    /*
    int improveQuality = 0;
    std::cout << "[Input] Improve quality? (Y/N): ";
    char qualityAnswer;
    std::cin >> qualityAnswer;
    if (qualityAnswer == 'Y' || qualityAnswer == 'y') {
        improveQuality = 1;
    }
    controller.setImproveQuality(improveQuality);
    */

    if (mode == 1) {
        controller.runPipelineInTesting();
    } else {
        controller.runPipeline();
    }

    std::cout << "[Complete] Program finished successfully.\n";
    return 0;
}
