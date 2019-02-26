EXE = softsim
SOURCES = globaldata.c initializer.c running.c jsmn.c main.c
HEADERS = globaldata.h initializer.h running.h jsmn.h

OUT_DIR = ../build
RESULT_DIRS = ../results/movies ../results/outs ../results/stats ../results/graphs

all: clean bin softsim

.PHONY: bin

bin: $(OUT_DIR) $(RESULT_DIRS)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

$(RESULT_DIRS):
	mkdir -p $(RESULT_DIRS)

$(EXE): $(SOURCES) $(HEADERS)
	gcc $^ -o $(OUT_DIR)/$@ -lm -ljsmn -L../jsmn

clean:
	rm -f $(OUT_DIR)/$(EXE)
