CC = g++

GL_INC_DIR = /usr/include
GL_LIB_DIR = /usr/lib

GL_LIBS = -L$(GL_LIB_DIR) -lglut -lGLU -lGL

CFLAGS = -Wall -pg -pedantic -O3
BOOST_FLAGS = -DNDEBUG -DBOOST_UBLAS_NDEBUG

OBJ = Grid.o main.o FluidSim.o Renderer.o Particles.o

.PHONY:clean

LIBS = $(GL_LIBS)
COMP_HEAEDERS = Vec.h MetaConfig.h

all: $(OBJ) $(TEST_CLASS_OBJ) sim
	
sim: $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $@

main.o: FluidSim.o Grid.o Renderer.o main.cpp Dimensions.h MetaConfig.h
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) main.cpp -o $@

FluidSim.o: Grid.o FluidSim.cpp FluidSim.h $(COMP_HEADERS)
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) FluidSim.cpp -o $@

Renderer.o: Grid.o Renderer.cpp Renderer.h $(COMP_HEADERS)
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) Renderer.cpp -o $@

Grid.o: Particles.h Particles.cpp Grid.cpp Grid.h GridTest.h $(COMP_HEADERS)
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) Grid.cpp -o $@
	
Particles.o: Particles.h Particles.cpp Vec.h
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) Particles.cpp -o $@

clean:
	rm -f  ./*~ ./sim $(OBJ) 


