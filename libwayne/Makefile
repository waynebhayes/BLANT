CC=gcc $(OPT) $(DEBUG) -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow $(PG)

default:
	# Make the same thing we made most recently
	if ls -rtlL *.a | tail -1 | grep .-g; then make debug; else make opt; fi

all:
	/bin/rm -f *.a
	# Make the pg versions (for profiling)
	make debug_clean
	make 'PG=-pg' -j4 debug
	make opt_clean
	make 'PG=-pg' -j4 opt
	for i in *.a; do b=`basename $$i .a`; mv $$i $$b-pg.a; done
	# Make the non-pg versions (for profiling)
	make debug_clean
	make -j4 debug
	make opt_clean
	make -j4 opt
	#(cd tests; make stats; mv stats ../bin)
	touch made

opt:
	make 'OPT=-O1' 'LIBOUT=libwayne.a' libwayne

debug:
	make 'DEBUG=-ggdb' 'LIBOUT=libwayne-g.a' libwayne

libwayne:
	make $(LIBOUT)
	mv src/$(LIBOUT) .
	[ `arch` = Darwin ] || ar r $(LIBOUT)

debug_clean:
	make 'DEBUG=-ggdb' 'LIBOUT=libwayne-g.a' raw_clean

opt_clean:
	make 'OPT=-O1' 'LIBOUT=libwayne.a' raw_clean

raw_clean:
	/bin/rm -f src/*.[oa] $(LIBOUT) made

clean:
	/bin/rm -f *.a
	make debug_clean
	make opt_clean
	/bin/rm -f made

$(LIBOUT): src/$(LIBOUT)
	ranlib src/$(LIBOUT)

src/$(LIBOUT):
	cd src; make 'CC=$(CC)' 'LIBOUT=$(LIBOUT)' 'DEBUG=$(DEBUG)'
