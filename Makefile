NO_RAW_PARSER ?=0
NO_GRAPHVIZ ?=0

FILE ?=

CC = gcc
INCDIR = ./include
LIBDIR = ./lib
LDLIBS = -lm
CFLAGS = -I$(INCDIR) -fPIC -O3 -Wall -Wno-unused-variable

ifeq ($(OS),Windows_NT)
  OS_DETECTED := Windows
else
  OS_DETECTED := $(shell uname -s)
endif
$(warning $(OS_DETECTED))
ifeq ($(OS_DETECTED),Darwin)
	LDFLAGS += -dynamiclib
else
	LDFLAGS += -shared
endif

# Raw parser
ifeq ($(NO_RAW_PARSER),1)
	CFLAGS += -DNO_RAW_PARSER
else
	LDLIBS += -lraw_parser
endif

# Graphviz
ifeq ($(NO_GRAPHVIZ),1)
	CFLAGS += -DNO_GRAPHVIZ
else
	LDLIBS += -lgvc -lcgraph
	CFLAGS += -I$(GRAPHVIZ)/include # -lgvc -lcgraph
endif

SOURCES_LIB = $(shell echo src/*/*.c src/*/*/*.c)
OBJECTS_LIB = $(SOURCES_LIB:.c=.o)
TARGET_LIB = $(LIBDIR)/libpfnet.so

SOURCES_TEST = $(shell echo tests/*.c)
OBJECTS_TEST = $(SOURCES_TEST:.c=.o)
TARGET_TEST = $(SOURCES_TEST:.c=.out)

.PHONY: all
all : lib test

.PHONY: lib
lib : $(TARGET_LIB)

$(TARGET_LIB) : $(OBJECTS_LIB)
	mkdir -p lib
	$(CC) $(CFLAGS) -o $@ $(OBJECTS_LIB) $(LDFLAGS)

.PHONY: test
test : $(TARGET_TEST)
tests/%.out: tests/%.c
	$(CC) $(CFLAGS) -L$(LIBDIR) -o $@ $< -lpfnet $(LDLIBS)
	./tests/run_tests.out ./data/ieee14.mat

.PHONY: single
single: 
ifndef FILE
	$(error error: 'FILE' must be specified)
else
	$(CC) -c $(CFLAGS) -Wfatal-errors $(FILE) -o $(basename $(FILE)).o
endif

.PHONY: docs
docs :
ifndef PFNET_DOCS
	$(error error: 'PFNET_DOCS' must be set to the location to put documentation files)
else
	doxygen ./docs/Doxyfile
endif

.PHONY: clean
clean :
	rm -f $(OBJECTS_LIB) $(TARGET_LIB)
	rm -f $(OBJECTS_TEST) $(TARGET_TEST)
