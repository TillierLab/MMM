NAME = mmmWithL
CC = g++
FLAGS = -O2 -Wall $(INCLUDE)
INCLUDE = -I include -I shiny -I ports
LDLIBS = -lm -lz -lboost_iostreams

#vpath %.o bin
vpath %.cpp src
vpath %.hpp include
vpath %.h include

OBJS = classCmdLineArgParser.o ProteinGroup.o Alignment.o DistanceMatrix.o mmm_algorithm.o MMML.o treeManip.o sighand.o hw_iface.o timer.o ports/tmports.o

all: MatrixMatchMaker.cpp $(OBJS) 
	$(CC) $(FLAGS) -g -pg $(OBJS) $< -o $(NAME)_debug $(LDLIBS) -D VERBOSE
	$(CC) $(FLAGS) $(OBJS) $< -o $(NAME)_verbose $(LDLIBS) -D VERBOSE
	$(CC) $(FLAGS) $(OBJS) $< -o $(NAME) $(LDLIBS)
mmm_algorithm.o: mmm_algorithm.cpp mmm_algorithm.h bonus.h
	$(CC) $(FLAGS) -c -o $@ $<

ProteinGroup.o: ProteinGroup.cpp ProteinGroup.hpp
	$(CC) $(FLAGS) -c -o $@ $<

Alignment.o: Alignment.cpp Alignment.hpp 
	$(CC) $(FLAGS) -c -o $@ $<

DistanceMatrix.o: DistanceMatrix.cpp DistanceMatrix.hpp
	$(CC) $(FLAGS) -c -o $@ $<

%.o: %.cpp %.h
	$(CC) $(FLAGS) -c -o $@ $<

clean:
	@rm -f $(OBJS) $(NAME) $(NAME)_verbose *~

pmbmaker: PmbMaker.cpp $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $< -o pmbmaker $(LDLIBS)
