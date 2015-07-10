DEBUG ?= 0
NO_RAW_PARSER ?=0
NO_GRAPHVIZ ?=0

CC = gcc
INCDIR = ./include
LIBDIR = ./lib
LDFLAGS = -shared
LDLIBS = -lpfnet -lm
CFLAGS = -I$(INCDIR) -fPIC -O3 -Wall -Wno-unused-variable

# Debug
ifeq ($(DEBUG),1)
	CFLAGS += -DDEBUG
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
	LDLIBS += -lgvc
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
	$(CC) $(CFLAGS) -L$(LIBDIR) -o $@ $< $(LDLIBS)
	./tests/run_tests.out ./data/ieee14.mat

.PHONY: clean
clean : 
	rm -f $(OBJECTS_LIB) $(TARGET_LIB)
	rm -f $(OBJECTS_TEST) $(TARGET_TEST)

