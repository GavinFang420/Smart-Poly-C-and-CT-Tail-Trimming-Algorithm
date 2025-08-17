CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
DEBUG_FLAGS = -g -DDEBUG

# Source files
SOURCES = src/smarttrim.cpp src/mergeread.cpp
HEADERS = src/smarttrim.h src/mergeread.h
TEST_SOURCES = tests/unit_test.cpp

# Object files (put in build directory)
BUILD_DIR = build
OBJECTS = $(addprefix $(BUILD_DIR)/, $(notdir $(SOURCES:.cpp=.o)))
TEST_OBJECTS = $(addprefix $(BUILD_DIR)/, $(notdir $(TEST_SOURCES:.cpp=.o)))

# Executable names
MAIN_EXEC = smarttrim
TEST_EXEC = unit_test

# Default target
all: $(MAIN_EXEC) $(TEST_EXEC)

# Main executable (if you want to create a standalone version later)
$(MAIN_EXEC): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Test executable
$(TEST_EXEC): $(OBJECTS) $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

# Create build directory and compile
$(BUILD_DIR)/%.o: src/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD_DIR)/%.o: tests/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(TEST_EXEC)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(MAIN_EXEC) $(TEST_EXEC)

# Run tests
test: $(TEST_EXEC)
	./$(TEST_EXEC)

# Run tests with valgrind (memory check)
memcheck: $(TEST_EXEC)
	valgrind --leak-check=full --show-leak-kinds=all ./$(TEST_EXEC)

# Install (copy to /usr/local/bin, requires sudo)
install: $(MAIN_EXEC)
	sudo cp $(MAIN_EXEC) /usr/local/bin/

# Uninstall
uninstall:
	sudo rm -f /usr/local/bin/$(MAIN_EXEC)

# Generate parameter matrix test (will create 100+ configurations)
benchmark: $(TEST_EXEC)
	./$(TEST_EXEC) --benchmark

# Format code (requires clang-format)
format:
	clang-format -i *.cpp *.h

# Check code style
lint:
	cppcheck --enable=all --std=c++11 *.cpp *.h

.PHONY: all debug clean test memcheck install uninstall benchmark format lint
