AFLAGS =  -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo
LFLAGS = -std=c++11 `pkg-config --libs lib/glfw3.pc` 
N=55


build/Collision.o: src/Collision.cpp src/Collision.h src/disk.h
	g++ -g -std=c++11 -O2 -c -o build/Collision.o src/Collision.cpp

build/circle.o: src/circle.cpp src/circle.h
	g++ -g -std=c++11 -O2 -c -o build/circle.o src/circle.cpp

build/cross.o: src/cross.cpp src/cross.h
	g++ -g -std=c++11 -O2 -c -o build/cross.o src/cross.cpp

build/Disks.o: src/Disks.cpp src/animation.h src/stats.h src/Collision.h src/animation.cpp src/stats.cpp src/Collision.cpp src/disk.h
	g++ -g -std=c++11 -O2 -c -o build/Disks.o src/Disks.cpp

build/animation.o: src/animation.cpp src/animation.h src/circle.h src/cross.h
	g++ -g -std=c++11 -O2 -c -o build/animation.o src/animation.cpp

build/RotationSim.o: src/RotationSim.cpp src/animation.h src/Disks.h src/Collision.h src/animation.cpp src/Collision.cpp src/Disks.cpp src/stats.h src/stats.cpp
	g++ -g -std=c++11 -O2 -c -o build/RotationSim.o src/RotationSim.cpp

build/stats.o: src/stats.cpp src/stats.h src/disk.h
	g++ -g -std=c++11 -O2 -c -o build/stats.o src/stats.cpp



compile: build/cross.o build/Collision.o build/circle.o build/Disks.o build/animation.o build/RotationSim.o build/stats.o
	g++ -g -O2 $(LFLAGS) $(AFLAGS) build/stats.o build/cross.o build/Collision.o build/Disks.o build/RotationSim.o build/animation.o build/circle.o build/lib/glew.o -o RotationSim

clean:
	rm build/*.o
test:
	./RotationSim t
animate:
	./RotationSim a ${N}
run:
	./RotationSim r ${N}
