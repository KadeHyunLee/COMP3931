# 컴파일러 설정
CXX = g++
CXXFLAGS = -std=c++20 -Xpreprocessor -fopenmp -I/opt/homebrew/include -I/opt/homebrew/Cellar/libomp/19.1.7/include -I$(shell brew --prefix open-mesh)/include
LDFLAGS = -L/opt/homebrew/lib -L/opt/homebrew/Cellar/libomp/19.1.7/lib -L$(shell brew --prefix open-mesh)/lib -lOpenMeshCore -lOpenMeshTools -lomp -Wl,-rpath,$(shell brew --prefix open-mesh)/lib

# 소스 파일 및 출력 파일
SRC = main.cpp MeshUtils.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = main

# 기본 빌드 명령
all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(OBJ) -o $(EXEC) $(LDFLAGS)

# 개별 파일 컴파일
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 실행
run: all
	./$(EXEC)

# 정리
clean:
	rm -f $(OBJ) $(EXEC)