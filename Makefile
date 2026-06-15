# Attempt to figure out what kind of machine we're on.
UNAME=$(shell uname -a | awk '{if(/CYGWIN/){V="CYGWIN"}else if(/Darwin/){if(/arm64/)V="arm64";else V="Darwin"}else if(/Linux/){V="Linux"}}END{if(V){print V;exit}else{print "unknown OS" > "/dev/stderr"; exit 1}}')

# Number of cores to use when invoking parallelism
ifndef CORES
    CORES := 2
endif

# Waywe needs gcc-6 on MacOS:
ifndef GCC_VER
    GCC_VER=$(shell echo $(UNAME) $(HOME) | awk '/Darwin/&&/Users.wayne/{V="-6"}END{if(V)print V;else{printf "using default gcc: " > "/dev/null"; exit 1}}')
endif
GCC=gcc$(GCC_VER) # gcc gcc-4.2 gcc-6 gcc-7 gcc-8 gcc-9 # Possibilities on Darwin
CXX=g++$(GXX_VER)
STACKSIZE=$(shell ($(GCC) -v 2>&1; uname -a) | awk '/CYGWIN/{print "-Wl,--stack,83886080"}/gcc-/{actualGCC=1}/Darwin/&&actualGCC{print "-Wl,-stack_size -Wl,0x5000000"}')
CC=$(GCC) $(OPT) $(GDB) $(DEBUG) -Wno-unused-variable -Wno-unused-but-set-variable -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow $(PG) $(STACKSIZE)
LIBWAYNE_HOME:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

# Safe to call whenever a previous full make call has succeeded, call this after a source file has changed.
incremental:
	cd src; $(MAKE) -f Makefile.incremental 'CC=$(GCC)'

make-successful:
	touch make-started # create this file so it is OLDER than .a files
	make all
	mv make-started make-successful

all:
	@if [ -f do-not-make ]; then echo "not making; assuming they already exist"; ls -l libwayne*.a; else /bin/rm -f *.a; $(MAKE) libwayne_all; fi
	$(MAKE) testlib


# Run 'make watch' in a separate terminal for auto-recompilation 
watch:
	ls src/*.c include/*.h | entr $(MAKE) incremental

libwayne_all: parallel
	/bin/rm -f *.a src/*.a
	$(CC) -o bin/intSizes intSizes.c && ./bin/intSizes > include/intSizes.h

	$(MAKE) libwayne_one VARIANT=libwayne       PG_FLAG=    DBG=  NDEBUG_FLAG=          OPT_FLAGS=-O2
	$(MAKE) libwayne_one VARIANT=libwayne-nd    PG_FLAG=    DBG=  NDEBUG_FLAG=-DNDEBUG  OPT_FLAGS=-O2
	$(MAKE) libwayne_one VARIANT=libwayne-g     PG_FLAG=    DBG=1 NDEBUG_FLAG=          OPT_FLAGS=
	$(MAKE) libwayne_one VARIANT=libwayne-pg    PG_FLAG=-pg DBG=  NDEBUG_FLAG=          OPT_FLAGS=-O2
	$(MAKE) libwayne_one VARIANT=libwayne-pg-nd PG_FLAG=-pg DBG=  NDEBUG_FLAG=-DNDEBUG  OPT_FLAGS=-O2
	$(MAKE) libwayne_one VARIANT=libwayne-pg-g  PG_FLAG=-pg DBG=1 NDEBUG_FLAG=          OPT_FLAGS=

empty :=
space := $(empty) $(empty)

VARIANT := $(subst $(space),-,$(strip \
	$(if $(PG_FLAG),pg) \
	$(if $(DBG),g) \
	$(if $(NDEBUG_FLAG),nd) \
	$(if $(OPT_FLAGS),$(subst -,,$(OPT_FLAGS))) \
))

VARIANT := libwayne$(if $(VARIANT),-$(VARIANT),)

libwayne_one:
	mkdir -p build/$(VARIANT)
	echo "$(strip $(OPT_FLAGS) $(PG_FLAG) $(if $(DBG),-ggdb -DDEBUG=1,$(NDEBUG_FLAG)))" > build/$(VARIANT)/.cflags
	$(MAKE) -j$(CORES) \
		'PG=$(PG_FLAG)' \
		'OPT=$(OPT_FLAGS)' \
		'GDB=$(if $(DBG),-ggdb)' \
		'DEBUG=$(if $(DBG),-DDEBUG=1,$(NDEBUG_FLAG))' \
		'OBJDIR=../build/$(VARIANT)' \
		'LIBOUT=$(VARIANT).a' \
		libwayne
	/bin/rm -f src/$(VARIANT).a
	mkdir -p build/.stamps && for o in build/$(VARIANT)/*.o; do touch build/.stamps/$$(basename $${o%.o}).stamp; done

parallel: parallel.c
	$(CC) -o bin/parallel parallel.c

testlib:
	export LIBWAYNE_HOME=$(LIBWAYNE_HOME); for x in ebm covar stats hash raw_hashmap htree-test avltree-test bintree-test CI graph-sanity tinygraph-sanity graph-weighted circ_buf; do rm -f bin/$$x tests/$$x.o; ( cd tests; $(MAKE) $$x; mv $$x ../bin; IN=/dev/null; [ -f $$x.in ] && IN=$$x.in; cat $$IN | ../bin/$$x $$x.in > /tmp/$$x.test$$$$ 2>&1 || exit 1; cat /tmp/$$x.test$$$$ | if [ -f $$x.out ]; then cmp - $$x.out; else wc; fi; /bin/rm -f /tmp/$$x.test$$$$); done #sim_anneal

opt:
	$(MAKE) 'OPT=-O2' 'DEBUG=-DNDEBUG' 'LIBOUT=libwayne.a' libwayne

debug:
	$(MAKE) 'GDB=-ggdb' 'DEBUG=-DDEBUG=1' 'LIBOUT=libwayne-g.a' libwayne

shared-opt:
	$(MAKE) 'OPT=-O2 -fPIC' 'LIBOUT=libwayne.a' libwayne

shared-debug:
	$(MAKE) 'OPT=-fPIC' 'GDB=-ggdb' 'DEBUG=-DDEBUG=1' 'LIBOUT=libwayne-g.a' libwayne

ndebug:
	$(MAKE) 'OPT=-O2' 'DEBUG=-DNDEBUG=1' 'LIBOUT=libwayne-nd.a' libwayne

libwayne:
	$(MAKE) $(LIBOUT)
	# add misc.o below since adding nothing fails on Mac M1, and misc.o is pretty much necessary for libwayne to work.
	[ "$(UNAME)" = Darwin ] || ar r $(LIBOUT) build/$(notdir $(OBJDIR))/misc.o

debug_clean:
	@$(MAKE) 'GDB=-ggdb' 'DEBUG=-DDEBUG=1' 'LIBOUT=libwayne-g.a' raw_clean

opt_clean:
	@$(MAKE) 'OPT=-O2' 'LIBOUT=libwayne.a' raw_clean

ndebug_clean:
	@$(MAKE) 'OPT=-O2' 'LIBOUT=libwayne-nd.a' raw_clean

raw_clean:
	rm -f make-successful
	@/bin/rm -f src/*.[oa] $(LIBOUT)
	@cd C++; $(MAKE) clean

clean:
	rm -f make-successful
	rm -rf build
	@/bin/rm -rf src/libwayne*/
	@# The following is meant to remove the non-Windows binary, ie stats but not stats.exe.
	@/bin/rm -f bin/stats bin/raw_hashmap bin/hash
	@/bin/rm -f *.a src/*.a
	@$(MAKE) debug_clean
	@$(MAKE) opt_clean

$(LIBOUT): src/$(LIBOUT)
	if ranlib src/$(LIBOUT); then :; else echo "ranlib failed but it's not crucial" >&2; fi
	cp -p src/$(LIBOUT) .

.PHONY: src/$(LIBOUT)
src/$(LIBOUT):
	cd src; $(MAKE) -f Makefile.all 'CC=$(CC)' 'LIBOUT=$(LIBOUT)' 'GDB=$(GDB)' 'OBJDIR=$(OBJDIR)' '$(LIBOUT)'

#.PHONY: libwayne-pg-g src/$(LIBOUT)