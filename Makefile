CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

SRC_DIR := src
INC_DIR := include
TARGET := bjorken_solver

SOURCES := \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/rhs.cpp \
	$(SRC_DIR)/rk4.cpp \
	$(SRC_DIR)/KiTable.cpp \
	$(SRC_DIR)/bickleynaylor.cpp \
	$(SRC_DIR)/utils.cpp

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $(SOURCES) -o $(TARGET)

run: $(TARGET)
	./$(TARGET) $(l) $(t0) $(tf) $(n)

clean:
	rm -f $(TARGET)
