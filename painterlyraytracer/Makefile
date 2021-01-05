# set up how compilation will happen
CC = g++
CFLAGS = -Wall -g
LDFLAGS =
LIBS = -lpng -lX11 -lpthread
DEPS = CImg.h lights.h objects.h painter.h particle.h ray.h style.h utilities.h vec.h

# define the set of files used
objects = pathtracer.o

# default 
all: ghibli

ghibli: $(objects) 
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) $< -o $@

# remove files created
clean:
	rm -f ghibli $(objects)

# targets that do not produce files
.PHONEY: all clean