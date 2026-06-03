# Makefile for ReaxTools without RDKIT
# 使用方法:
#   make          - 编译 Release 版本
#   make debug    - 编译 Debug 版本
#   make clean    - 清理编译文件
#   make install  - 安装到系统目录 (可选)

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -I. -DFMT_HEADER_ONLY -DGRAPHVIZ 
LDFLAGS = -pthread -lgvc -lcgraph -lcdt 
TARGET = reax_tools.exe

# 源文件
SOURCES = $(wildcard *.cpp)


# 生成的目标文件
OBJECTS = $(SOURCES:.cpp=.o)

# 头文件依赖
HEADERS = $(SOURCES:.cpp=.h)

# 默认目标: Release
all: $(TARGET)

# 编译可执行文件
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)
	@echo "=========================================="
	@echo "编译成功！生成: $(TARGET)"
	@echo "使用方法: ./$(TARGET) --traj input.xyz --output output_dir"
	@echo "=========================================="

# 编译对象文件
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Debug 版本
debug: CXXFLAGS = -std=c++17 -g -DWASM_MODE -Wall -Wextra
debug: clean all

# 清理
clean:
	rm -f $(OBJECTS) $(TARGET)
	@echo "清理完成"

# 安装到系统路径 (可选，需要管理员权限)
install: $(TARGET)
	@echo "安装 $(TARGET) 到 /usr/local/bin/ 需要管理员权限"
	sudo cp $(TARGET) /usr/local/bin/ || echo "安装失败，请以管理员身份运行"

# 帮助信息
help:
	@echo "ReaxTools Makefile"
	@echo ""
	@echo "可用目标:"
	@echo "  make        - 编译 Release 版本"
	@echo "  make debug  - 编译 Debug 版本 (包含调试信息)"
	@echo "  make clean  - 清理编译生成的文件"
	@echo "  make help   - 显示此帮助信息"
	@echo ""
	@echo "依赖要求:"
	@echo "  - g++ (MinGW-w64 或 GCC)"
	@echo "  - pthread 库"
	@echo "  - fmt 库 (如果需要完整功能，去掉 -DWASM_MODE)"
	@echo ""
	@echo "示例用法:"
	@echo "  ./$(TARGET) --traj trajectory.xyz --output results --threads 8"

# 声明伪目标
.PHONY: all debug clean install help
