# -*- Makefile -*-
PROGRAM = mcalc
CXXFILES = $(wildcard src/*.cpp)
CFILES = $(wildcard src/*.c)
CXXFLAGS = `pkg-config --cflags gtk+-3.0`
LIBS = `pkg-config --libs gtk+-3.0`
COBJECTS = $(patsubst src/%.c, %.o, $(CFILES))
CXXOBJECTS = $(patsubst src/%.cpp, %.o, $(CXXFILES))

debug:	CXXFLAGS+=-g -D_DEBUG
debug:	$(PROGRAM)

release:	CXXFLAGS+=-O
release:	$(PROGRAM)

$(COBJECTS): $(CFILES)
	$(CC) $(CXXFLAGS) -march=core2 -c $?

$(CXXOBJECTS): $(CXXFILES)
	$(CXX) $(CXXFLAGS) -march=core2 -c $?

$(PROGRAM): $(COBJECTS) $(CXXOBJECTS)
	$(CXX) -no-pie $(COBJECTS) $(CXXOBJECTS) -o $(PROGRAM) $(LIBS)

#$(PROGRAM): $(CFILES) $(CXXFILES)
#	$(CXX) -no-pie $(CXXFLAGS) -march=core2 $(CFILES) $(CXXFILES) -o $(PROGRAM) $(LIBS)

run:
	./$(PROGRAM)

clean:
	rm $(COBJECTS) $(CXXOBJECTS) $(PROGRAM)
