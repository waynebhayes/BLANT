CC=gcc $(OPT) $(DEBUG) -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow #-pg

default:
	# Make the same thing we made most recently
	if ls -rtlL *.a | tail -1 | grep .-g; then make debug; else make opt; fi

all:
	make opt_clean
	make -j4 opt
	make debug_clean
	make -j4 debug
	(cd tests; make stats; mv stats ../bin)

opt:
	make 'OPT=-O1' 'LIBOUT=libwayne.a' libwayne

debug:
	make 'DEBUG=-ggdb' 'LIBOUT=libwayne-g.a' libwayne

libwayne:
	make $(LIBOUT)
	mv src/$(LIBOUT) .
	ar r $(LIBOUT)

debug_clean:
	make 'DEBUG=-ggdb' 'LIBOUT=libwayne-g.a' raw_clean

opt_clean:
	make 'OPT=-O1' 'LIBOUT=libwayne.a' raw_clean

raw_clean:
	/bin/rm -f src/*.[oa] $(LIBOUT)

clean:
	make debug_clean
	make opt_clean

$(LIBOUT): src/$(LIBOUT)
	ranlib src/$(LIBOUT)

src/$(LIBOUT):
	cd src; make 'CC=$(CC)' 'LIBOUT=$(LIBOUT)' 'DEBUG=$(DEBUG)'
