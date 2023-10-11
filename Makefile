CPP=g++
CFLAGS=-O2
DIR_PATH = .
SOURCES = $(wildcard $(DIR_PATH)/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: 
	projetprelim projet	projetopt

projetprelim: $(eval DIR_PATH = ./projet_prelim) $(OBJECTS)
	$(CPP) $(CFLAGS) -o $@ $(eval DIR_PATH = ./projet_prelim) $(OBJECTS)

projet: $(eval DIR_PATH = ./projet_) $(OBJECTS)
	$(CPP) $(CFLAGS) -o $@ $(eval DIR_PATH = ./projet_) $(OBJECTS)
	
projetopt: $(eval DIR_PATH = ./projet_opt) $(OBJECTS)
	$(CPP) $(CFLAGS) -o $@ $(eval DIR_PATH = ./projet_opt) $(OBJECTS)

%.o: $(DIR_PATH)/%.cpp $(DIR_PATH)/%.hpp
	$(CPP) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm $(wildcard */*.o)
