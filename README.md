# COMP3931 - Parallel Mesh Simplification on macOS

This project implements a parallel 3D mesh simplification pipeline using uniform grid-based spatial partitioning and OpenMP. It is designed and tested exclusively on macOS (Apple Silicon, M1/M2) and uses a `Makefile` for building the project.

## Environment
- OS: macOS 
- Chip: Apple M2 Pro 
- Build System: Makefile 

## Dependencies
- [OpenMesh]
- [OpenMP] 

## Build Instructions
1. Install OpenMesh and OpenMP:

   brew install open-mesh 
   brew install llvm

2. Build the project:
   make

## Execution 

1.Download test model
https://leeds365-my.sharepoint.com/:u:/g/personal/sc22cl_leeds_ac_uk/EcFzkak9L_pMjvd49F1Cb9wBN-aoY-2mRl80i0oEh3Nx_A?e=hIyn0g

2. Please move the file to data foldera and unzip the compressed file in the directory before running

./main -> 0(mode) -> model data name (eg. bunny.ply, dragon.ply etc) -> deciamtion ratio (eg. 0.1 or 0.2 (recommend))
./main -> 1(mode) -> model data name (eg. bunny.ply, dragon.ply etc) -> decimation ratio eg. 0.1 or 0.2 (recommend)-> number of threds you want to use(eg. 10) 

## Notes
- Tested only on Apple Silicon (macOS). Performance may vary on Intel or other architectures on other OS


