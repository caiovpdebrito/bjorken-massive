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
	$(SRC_DIR)/thermodynamics.cpp \
	$(SRC_DIR)/utils.cpp

TABLE_SOURCES := \
	$(SRC_DIR)/generate_Ki_table.cpp \
	$(SRC_DIR)/bickleynaylor.cpp

.PHONY: help all run table plot clean

.PHONY: help

help:
	@echo "Available targets:"
	@echo "  build   Build the project"
	@echo "  run     Execute ODEs solver"
	@echo "  table   Generate table of Bickley-Naylor functions"
	@echo "  plot    Create plots with matplotlib"
	@echo "  clean   Remove generated files"

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(SOURCES) -o $(TARGET)

$(TABLE_TARGET): $(TABLE_SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(TABLE_SOURCES) -o $(TABLE_TARGET)

#==========================================#
#=========== EXECUTE ODE SOLVER ===========#
#==========================================#

run: $(TARGET)
	./$(TARGET) $(l) $(w0) $(wf) $(n)

#==========================================#
#======== GENERATE TABLE OF ĸ_n(z) ========#
#==========================================#

table_file ?= data/Ki_table.dat
npoints ?= 1000

table: $(TABLE_TARGET)
	mkdir -p data
	./$(TABLE_TARGET) $(nmin) $(nmax) $(xmin) $(xmax) $(npoints) $(table_file)

#==========================================#
#============== PLOT MOMENTS ==============#
#==========================================#

plot:
	python3 $(GRPH_DIR)/plot_moments.py $(l)

#==========================================#
#========== REMOVE CREATED FILES ==========#
#==========================================#

clean:
	rm -f $(TARGET) $(TABLE_TARGET)

