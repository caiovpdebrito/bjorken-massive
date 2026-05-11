CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

SRC_DIR := src
INC_DIR := include
GRPH_DIR := graphics

TARGET := bjorken_solver
TABLE_TARGET := generate_Ki_table
SCAN_TARGET := mass_scan 

SOURCES := \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/rhs.cpp \
	$(SRC_DIR)/rk4.cpp \
	$(SRC_DIR)/KiTable.cpp \
	$(SRC_DIR)/bickleynaylor.cpp \
	$(SRC_DIR)/thermodynamics.cpp \
	$(SRC_DIR)/utils.cpp

SCAN_SOURCES := \
	$(SRC_DIR)/mass_scan.cpp \
	$(SRC_DIR)/rhs.cpp \
	$(SRC_DIR)/rk4.cpp \
	$(SRC_DIR)/KiTable.cpp \
	$(SRC_DIR)/bickleynaylor.cpp \
	$(SRC_DIR)/thermodynamics.cpp \
	$(SRC_DIR)/utils.cpp

TABLE_SOURCES := \
	$(SRC_DIR)/generate_Ki_table.cpp \
	$(SRC_DIR)/bickleynaylor.cpp

.PHONY: help all run table plot mscan plot-mass clean

help:
	@echo "Available targets:"
	@echo "  build			Build the project"
	@echo "  table   		Generate table of Bickley-Naylor functions"
	@echo "  run     		Execute ODEs solver"
	@echo "  plot    		Create plots with matplotlib"
	@echo "  mscan   		Solve ODEs for a range of masses"
	@echo "  plot-mass  	Create plots with matplotlib"
	@echo "  clean   		Remove generated files"

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(SOURCES) -o $(TARGET)

$(SCAN_TARGET): $(SCAN_SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(SCAN_SOURCES) -o $(SCAN_TARGET)

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
#======== MASS SCAN: RUN AND PLOTS ========#
#==========================================#

# Initial conditions
tau0 ?= 0.5
tau ?= 10
mmin ?= 0.0
mmax ?= 1.0

# Truncation parameters
mode ?= n

ifeq ($(mode),n)
nfirst ?= 0
nlast ?= 10
lfirst ?= 1
llast ?= $(lfirst)
else ifeq ($(mode),l)
nfirst ?= 1
nlast ?= $(nfirst)
lfirst ?= 1
llast ?= 10
endif

# Numerical parameters
steps ?= 1000
nmass ?= 10
mass_file ?= data/mass_scan_$(mode).dat

mscan: $(SCAN_TARGET)
	mkdir -p data
	./$(SCAN_TARGET) $(tau0) $(tau) $(steps) $(mmin) $(mmax) $(nmass) $(mode) $(nfirst) $(nlast) $(lfirst) $(llast) $(mass_file)

plot-mass:
	python3 $(GRPH_DIR)/plot_mass_scan.py --input $(mass_file) --tau-label $(tau)

#==========================================#
#========== REMOVE CREATED FILES ==========#
#==========================================#

clean:
	rm -f $(TARGET) $(TABLE_TARGET) $(SCAN_TARGET)
