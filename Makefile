# set up how compilation will happen
CC = g++
CFLAGS = -Wall -g
LDFLAGS =
LIBS = -lpng -lX11 -lpthread
DEPS = CImg.h brush.h light.h material.h object.h painter.h particle.h plane.h ray.h sphere.h stroke.h style.h tri.h utilities.h vec3.h vertex.h

# define the set of files used
objects = raytracer.o light.o material.o plane.o ray.o sphere.o tri.o vec3.o vertex.o

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