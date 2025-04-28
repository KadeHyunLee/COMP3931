
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
    
    CXXFLAGS = -std=c++20 -Xpreprocessor -fopenmp -I/opt/homebrew/include -I/opt/homebrew/Cellar/libomp/19.1.7/include -I$(shell brew --prefix open-mesh)/include -I./nanoflann/include
    LDFLAGS = -L/opt/homebrew/lib -L/opt/homebrew/Cellar/libomp/19.1.7/lib -L$(shell brew --prefix open-mesh)/lib -lOpenMeshCore -lOpenMeshTools -lomp -Wl,-rpath,$(shell brew --prefix open-mesh)/lib
else

    CXXFLAGS = -std=c++20 -fopenmp -I/usr/local/include -I./nanoflann/include
    LDFLAGS = -L/usr/local/lib -lOpenMeshCore -lOpenMeshTools -lomp
endif

CXX = g++

SRC = main.cpp MeshUtils.cpp MeshPipelineController.cpp MeshLoader.cpp GridPartitioner.cpp SubmeshExtractor.cpp SubmeshDecimator.cpp MeshIntegrator.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = main

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(OBJ) -o $(EXEC) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	./$(EXEC)

clean:
	rm -f $(OBJ) $(EXEC)