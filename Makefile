CC = g++

GL_INC_DIR = /usr/include
GL_LIB_DIR = /usr/lib/include
GL_LIBS = -L$(GL_LIB_DIR) -lglut -lGLU -lGL

OPTIFLAGS = -pedantic -O3
CFLAGS = -Wall -pg 
BOOST_FLAGS = -DNDEBUG -DBOOST_UBLAS_NDEBUG

CLASSES_O = Grid.o FluidSim.o Renderer.o Particles.o sim.o
MODULES = $(CLASSES_O) main.o config.tab.o config.yy.o 	

.PHONY:clean

LIBS = $(GL_LIBS)
COMP_HEAEDERS = Vec.h MetaConfig.h sim.h

all: $(MODULES) $(TEST_CLASS_OBJ) firesim
	
firesim: $(MODULES)
	$(CC) $(MODULES) $(LIBS) -o $@

main.o: FluidSim.o Grid.o Renderer.o main.cpp MetaConfig.h sim.o
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) $(OPTIFLAGS) main.cpp -o $@

FluidSim.o: Grid.o FluidSim.cpp FluidSim.h $(COMP_HEADERS)
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) $(OPTIFLAGS) FluidSim.cpp -o $@

Renderer.o: Grid.o Renderer.cpp Renderer.h $(COMP_HEADERS)
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) $(OPTIFLAGS) Renderer.cpp -o $@

Grid.o: Particles.h Particles.cpp Grid.cpp Grid.h GridTest.h $(COMP_HEADERS)
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) $(OPTIFLAGS) Grid.cpp -o $@
	
Particles.o: Particles.h Particles.cpp Vec.h sim.h
	$(CC) -c $(CFLAGS) $(BOOST_FLAGS) $(OPTIFLAGS) Particles.cpp -o $@
	
sim.o: sim.h sim.cpp
	$(CC) -c $(CFLAGS) $(OPTIFLAGS) sim.cpp -o $@
	

config.yy.c : config.l config.h
	lex -o config.yy.c config.l

config.tab.c : config.y
	yacc -o config.tab.c --defines=config.tab.h config.y

config.tab.o: config.tab.c
	gcc -c -w $(CFLAGS) $< -o $@

config.yy.o: config.yy.c
	gcc -c -w $(CFLAGS) $< -o $@
	

clean:
	rm -f  ./*~ 
	rm -f ./firesim 
	rm -f $(MODULES) 
	rm -f config.*.* 


