CXX = g++
CXXFLAGS = -std=c++11 -O3 -pthread
SOURCES = src/main.cpp src/fastqreader.cpp src/writer.cpp src/read.cpp src/sequence.cpp
TARGET = smarttrim

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) -lz

clean:
	rm -f $(TARGET)
