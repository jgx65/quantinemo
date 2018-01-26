#CC=icpc -static -O3
#CC=g++ -static #-O3
#CFLAGS= -O3 #-static
#CFLAGS_DB="" #-static
CC=g++ # for mac os x

SOURCES=$(shell ls src/*.cpp)
OBJECTS=$(shell for file in $(SOURCES);\
do echo $$file | sed -e "s/\(.*\)\.cpp/\1\.o/"; echo " ";\
done)

GIT_VERSION := $(shell git describe --abbrev=7 --always --tags --dirty)

PRGNAME=quantiNemo2
PRGDIR=bin

INCS      =     -I.  -D_SHOW_MEMORY

.PHONY: all clean debug release profile

all:   touch release

#so that git version and compiled time is always correct
touch:
	touch src/tsim_manager.cpp  

thread: CFLAGS += -std=c++11
thread: INCS += -D_THREAD
thread: release

debug: INCS += -D_DEBUG
debug: CFLAGS  += -Wall
debug: bin

release: CFLAGS  += -O3 #-static
release: CFLAGS  += -DVERSIONGIT=\"$(GIT_VERSION)\"
release: bin

profile: CFLAGS += -pg
profile: bin

bin : CFLAGS  += -Iinclude
bin : objects
	echo $(OBJECTS)
	$(CC)$(CFLAGS) src/*.o -o $(PRGDIR)/$(PRGNAME) $(LIBS)
	
objects : $(OBJECTS)

%.o : %.cpp
	$(CC)$(CFLAGS) -c -w $< -o $@ $(INCS)

clean:
	rm -f src/*.o $(PRGNAME) 
	

depend:
	$(CC)$(CFLAGS) -M *.cpp > $@

-include depend



