CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

SRC_DIR := src
INC_DIR := include
GRPH_DIR := graphics

TARGET := bjorken_solver
TABLE_TARGET := generate_Ki_table

SOURCES := \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/rhs.cpp \
	$(SRC_DIR)/rk4.cpp \
	$(SRC_DIR)/KiTable.cpp \
	$(SRC_DIR)/bickleynaylor.cpp \
	$(SRC_DIR)/utils.cpp

TABLE_SOURCES := \
	$(SRC_DIR)/generate_Ki_table.cpp \
	$(SRC_DIR)/bickleynaylor.cpp

table_file ?= data/Ki_table.dat

.PHONY: all run bickleynaylor plot clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(SOURCES) -o $(TARGET)

$(TABLE_TARGET): $(TABLE_SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(TABLE_SOURCES) -o $(TABLE_TARGET)

#run: $(TARGET)
#	./$(TARGET) $(l) $(t0) $(tf) $(n)

run: $(TARGET)
	./$(TARGET) $(l) $(w0) $(wf) $(n)

table: $(TABLE_TARGET)
	mkdir -p data
	./$(TABLE_TARGET) $(nmin) $(nmax) $(xmin) $(xmax) $(npoints) $(table_file)

plot:
	python3 $(GRPH_DIR)/plot_moments.py $(l)

clean:
	rm -f $(TARGET)
