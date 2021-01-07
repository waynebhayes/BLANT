# Number of cores to use when invoking parallelism
ifndef CORES
    CORES := 2
endif

CC=gcc $(OPT) $(DEBUG) -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow $(PG)

default:
	# Make the same thing we made most recently
	if ls -rtlL *.a | tail -1 | grep .-g; then $(MAKE) debug; else $(MAKE) opt; fi

all:
	/bin/rm -f *.a
	# Make the pg versions (for profiling)
	$(MAKE) debug_clean
	$(MAKE) -j$(CORES) 'PG=-pg' debug
	$(MAKE) opt_clean
	$(MAKE) -j$(CORES) 'PG=-pg' opt
	for i in *.a; do b=`basename $$i .a`; mv $$i $$b-pg.a; done
	# Make the non-pg versions (for profiling)
	$(MAKE) debug_clean
	$(MAKE) -j$(CORES) debug
	$(MAKE) opt_clean
	$(MAKE) -j$(CORES) opt
	for x in stats hashtest htree-test; do if [ ! -x bin/$$x ]; then (cd tests; $(MAKE) $$x; mv $$x ../bin; [ -f $$x.in ] && cat $$x.in | ../bin/$$x); fi; done
	touch made

opt:
	$(MAKE) 'OPT=-O3' 'LIBOUT=libwayne.a' libwayne

debug:
	$(MAKE) 'DEBUG=-ggdb' 'LIBOUT=libwayne-g.a' libwayne

libwayne:
	$(MAKE) $(LIBOUT)
	mv src/$(LIBOUT) .
	[ `arch` = Darwin ] || ar r $(LIBOUT)

debug_clean:
	@$(MAKE) 'DEBUG=-ggdb' 'LIBOUT=libwayne-g.a' raw_clean

opt_clean:
	@$(MAKE) 'OPT=-O3' 'LIBOUT=libwayne.a' raw_clean

raw_clean:
	@/bin/rm -f src/*.[oa] $(LIBOUT) made
	@cd MT19937; $(MAKE) clean

clean:
	@# The following is meant to remove the non-Windows binary, ie stats but not stats.exe.
	@/bin/rm -f bin/stats bin/hashtest
	@/bin/rm -f *.a
	@$(MAKE) debug_clean
	@$(MAKE) opt_clean
	@/bin/rm -f made

$(LIBOUT): src/$(LIBOUT)
	ranlib src/$(LIBOUT)

src/$(LIBOUT):
	cd src; $(MAKE) 'CC=$(CC)' 'LIBOUT=$(LIBOUT)' 'DEBUG=$(DEBUG)'
