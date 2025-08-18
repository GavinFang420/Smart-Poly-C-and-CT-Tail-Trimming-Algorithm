# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
DEBUG_FLAGS = -std=c++17 -Wall -Wextra -g -DDEBUG
RELEASE_FLAGS = -std=c++17 -Wall -Wextra -O3 -DNDEBUG

# Directories
SRC_DIR = src
BUILD_DIR = build
TEST_DIR = tests

# Source files (based on your actual source files in src/)
SOURCES = mergeread.cpp smarttrim.cpp
HEADERS = mergeread.h smarttrim.h

# Add src/ prefix to sources and headers
SRC_FILES = $(SOURCES:%=$(SRC_DIR)/%)
HEADER_FILES = $(HEADERS:%=$(SRC_DIR)/%)

# Object files
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)

# Test files
TEST_SOURCES = try_Ctail.cpp
TEST_FILES = $(TEST_SOURCES:%=$(TEST_DIR)/%)

# Executables
MAIN_TARGET = ctail_trimmer
TEST_TARGET = try_ctail_test

# Default target
all: $(TEST_TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Test executable (links with your source files)
$(TEST_TARGET): $(BUILD_DIR) $(OBJECTS) $(BUILD_DIR)/try_Ctail.o
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(BUILD_DIR)/try_Ctail.o -o $@

# Main executable (if you have a main.cpp later)
$(MAIN_TARGET): $(BUILD_DIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

# Compile source files from src/ directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) -c $< -o $@

# Compile test files from tests/ directory
$(BUILD_DIR)/try_Ctail.o: $(TEST_DIR)/try_Ctail.cpp $(HEADER_FILES) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) -c $< -o $@

# Debug build
debug: CXXFLAGS = $(DEBUG_FLAGS)
debug: $(TEST_TARGET)

# Release build
release: CXXFLAGS = $(RELEASE_FLAGS)
release: $(TEST_TARGET)

# Run tests
test: $(TEST_TARGET)
	./$(TEST_TARGET)

# Run tests with valgrind (if available)
test-memory: $(TEST_TARGET)
	valgrind --leak-check=full --track-origins=yes ./$(TEST_TARGET)

# Extract test data from FASTQ files (adjusted for Windows Git Bash)
prepare-test-data:
	@echo "Extracting test data from FASTQ files..."
	@if [ -f "/c/Users/fangy/Desktop/Ctail/normal_R1.fastq" ]; then \
		head -n 200 "/c/Users/fangy/Desktop/Ctail/normal_R1.fastq" > normal_R1_200.txt; \
		echo "Created normal_R1_200.txt"; \
	else \
		echo "Warning: /c/Users/fangy/Desktop/Ctail/normal_R1.fastq not found"; \
	fi
	@if [ -f "/c/Users/fangy/Desktop/Ctail/normal_R2.fastq" ]; then \
		head -n 200 "/c/Users/fangy/Desktop/Ctail/normal_R2.fastq" > normal_R2_200.txt; \
		echo "Created normal_R2_200.txt"; \
	else \
		echo "Warning: /c/Users/fangy/Desktop/Ctail/normal_R2.fastq not found"; \
	fi
	@if [ -f "/c/Users/fangy/Desktop/Ctail/tumor_R1.fastq" ]; then \
		head -n 200 "/c/Users/fangy/Desktop/Ctail/tumor_R1.fastq" > tumor_R1_200.txt; \
		echo "Created tumor_R1_200.txt"; \
	else \
		echo "Warning: /c/Users/fangy/Desktop/Ctail/tumor_R1.fastq not found"; \
	fi
	@if [ -f "/c/Users/fangy/Desktop/Ctail/tumor_R2.fastq" ]; then \
		head -n 200 "/c/Users/fangy/Desktop/Ctail/tumor_R2.fastq" > tumor_R2_200.txt; \
		echo "Created tumor_R2_200.txt"; \
	else \
		echo "Warning: /c/Users/fangy/Desktop/Ctail/tumor_R2.fastq not found"; \
	fi

# Alternative method using PowerShell (if Git Bash head doesn't work)
prepare-test-data-ps:
	@echo "Extracting test data using PowerShell..."
	powershell -Command "if (Test-Path 'C:\\Users\\fangy\\Desktop\\Ctail\\normal_R1.fastq') { Get-Content 'C:\\Users\\fangy\\Desktop\\Ctail\\normal_R1.fastq' -Head 200 | Out-File 'normal_R1_200.txt' -Encoding utf8; Write-Host 'Created normal_R1_200.txt' }"
	powershell -Command "if (Test-Path 'C:\\Users\\fangy\\Desktop\\Ctail\\normal_R2.fastq') { Get-Content 'C:\\Users\\fangy\\Desktop\\Ctail\\normal_R2.fastq' -Head 200 | Out-File 'normal_R2_200.txt' -Encoding utf8; Write-Host 'Created normal_R2_200.txt' }"
	powershell -Command "if (Test-Path 'C:\\Users\\fangy\\Desktop\\Ctail\\tumor_R1.fastq') { Get-Content 'C:\\Users\\fangy\\Desktop\\Ctail\\tumor_R1.fastq' -Head 200 | Out-File 'tumor_R1_200.txt' -Encoding utf8; Write-Host 'Created tumor_R1_200.txt' }"
	powershell -Command "if (Test-Path 'C:\\Users\\fangy\\Desktop\\Ctail\\tumor_R2.fastq') { Get-Content 'C:\\Users\\fangy\\Desktop\\Ctail\\tumor_R2.fastq' -Head 200 | Out-File 'tumor_R2_200.txt' -Encoding utf8; Write-Host 'Created tumor_R2_200.txt' }"

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(MAIN_TARGET) $(TEST_TARGET)
	rm -f *.o *.exe

# Clean test data
clean-test-data:
	rm -f *_200.txt

# Full clean
clean-all: clean clean-test-data

# Install (optional)
install: $(TEST_TARGET)
	cp $(TEST_TARGET) /usr/local/bin/

# Help
help:
	@echo "Available targets:"
	@echo "  all                  - Build test program (default)"
	@echo "  debug                - Build with debug flags"
	@echo "  release              - Build with optimization flags"
	@echo "  test                 - Run test suite"
	@echo "  test-memory          - Run tests with memory checking"
	@echo "  prepare-test-data    - Extract first 200 lines from FASTQ files (Git Bash)"
	@echo "  prepare-test-data-ps - Extract first 200 lines using PowerShell"
	@echo "  clean                - Remove build files"
	@echo "  clean-test-data      - Remove extracted test data"
	@echo "  clean-all            - Remove all generated files"
	@echo "  install              - Install to /usr/local/bin"
	@echo "  help                 - Show this help"

# Phony targets
.PHONY: all debug release test test-memory prepare-test-data prepare-test-data-ps clean clean-test-data clean-all install help
