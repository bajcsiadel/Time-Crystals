EXE = softsim

LIB_SOURCES = jsmn.c
LIB_HEADERS = jsmn.h

MY_SOURCES = globaldata.c initializer.c running.c utils.c main.c
MY_HEADERS = globaldata.h initializer.h running.h utils.h

OUT_DIR = ../build
RESULT_DIRS = ../results/movies ../results/outs ../results/stats ../results/graphs
LIBS = -lm -ljsmn -L../jsmn

INDENT_GENERAL_OPTIONS = -bl -bli0 -nce -nbad -ss -cdw  -cli4 -cbi0 -cs -npsl -bls -blf -lp -npcs -ut -ts4 -i4
INDENT_HEADER_OPTIONS = -di8 -ppi4

ifeq ($(OS),Windows_NT) #WINDOWS
	ECHO_MESSAGE = "Windows"

	GREEN = "[42m"
	DEFAULT = "[0m"
else
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux) #LINUX
	ECHO_MESSAGE = "Linux"
endif

ifeq ($(UNAME_S), Darwin) #APPLE
	ECHO_MESSAGE = "Mac OS X"
endif

GREEN = "\033[0;32m"
DEFAULT = "\033[0m"
endif

all: bin $(EXE)

.PHONY: bin indent

bin: $(OUT_DIR) $(RESULT_DIRS)

$(OUT_DIR):
	@mkdir -p $(OUT_DIR)
	@echo "Binary (named build) directory created"

$(RESULT_DIRS):
	@mkdir -p $(RESULT_DIRS)
	@echo "Result directory created\n"

$(EXE): $(MY_SOURCES) $(MY_HEADERS) $(LIB_SOURCES) $(LIB_HEADERS)
	@$(CC) $^ -o $(OUT_DIR)/$@ $(LIBS)
	@echo "Create \t" $(GREEN) $@ $(DEFAULT) "\n"
	@echo $(GREEN)"Build complete for $(ECHO_MESSAGE)\n"$(DEFAULT)

clean:
	@rm -f $(OUT_DIR)/$(EXE)
	@echo $(GREEN)"Existing binary files removed\n"$(DEFAULT)

indent:
	@indent $(INDENT_GENERAL_OPTIONS) $(INDENT_HEADER_OPTIONS) $(MY_HEADERS)
	@indent $(INDENT_GENERAL_OPTIONS) $(MY_SOURCES)
	@echo "Headers and sources indented"

install_libraries:
	sudo apt-get update
	sudo apt-get install gcc
	sudo apt-get install indent
	@echo $(GREEN)"Libraries has been installed!\n"$(DEFAULT)