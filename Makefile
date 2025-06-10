CXX = g++
CXXFLAGS = -Wall -std=c++17 -g  # Added -g for debugging symbols

SRC_DIRS = src/mimiq src/mimiq/qcircuit src/Experiments
SRC_MAIN = src/main.cpp

SOURCES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp)) $(SRC_MAIN)

OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))


TARGET = build/app

all: $(TARGET)

build:
	mkdir -p build/mimiq build/Experiments

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

clean:
	rm -rf build/app build/main.o build/mimiq/*.o build/Experiments/*.o

.PHONY: all clean
