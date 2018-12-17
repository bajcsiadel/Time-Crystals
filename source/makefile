all: clean bin softsim
bin: 
	mkdir -p ../build
softsim: globaldata.h globaldata.c initializer.h initializer.c running.h running.c main.c
	gcc globaldata.c initializer.c running.c jsmn.c jsmn.h main.c -o ../build/softsim -lm -ljsmn -L../jsmn
clean:
	rm -f ../build/softsim
