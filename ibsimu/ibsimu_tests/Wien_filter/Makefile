CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

wien: wien.o
	$(CC) -o wien wien.o $(LDFLAGS)
wien.o: wien.cpp

clean:
	$(RM) *.o wien 

