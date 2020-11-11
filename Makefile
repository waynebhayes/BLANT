# Number of cores to use when invoking parallelism
ifndef CORES
    CORES := 2
endif

# Some architectures, eg CYGWIN 32-bit and MacOS("Darwin") need an 80MB stack.
STACKSIZE=$(shell arch | awk '/CYGWIN/{print "-Wl,--stack,83886080"}/Darwin/{print "-Wl,-stack_size -Wl,0x5000000"}')
export LIBWAYNE_HOME=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/libwayne
LIBWAYNE_OPTS=-O3 -I $(LIBWAYNE_HOME)/include -L $(LIBWAYNE_HOME) -lwayne -lm $(STACKSIZE) # -static OPTIMIZED
#LIBWAYNE_OPTS=-O0 -I $(LIBWAYNE_HOME)/include -L $(LIBWAYNE_HOME) -lwayne-g  -lm -ggdb $(STACKSIZE) # for debugging
#LIBWAYNE_OPTS=-I $(LIBWAYNE_HOME)/include -L $(LIBWAYNE_HOME) -lwayne-pg -lm -pg  # for profiling

# Name of BLANT source directory
SRCDIR = src
# Put all c files in SRCDIR below.
BLANT_SRCS = blant.c \
			 blant-window.c \
			 blant-output.c \
			 blant-kovacs.c \
			 blant-utils.c \
			 blant-sampling.c \
			 blant-synth-graph.c

OBJDIR = _objs
OBJS = $(addprefix $(OBJDIR)/, $(BLANT_SRCS:.c=.o))
CC=gcc -O3 #-ggdb
CXX=g++

### Generated File Lists ###
EIGHT := 8# COMMENT OUT THIS LINE to save "make" time (and disable k=8 sized graphlets)
SEVEN := 7#
K := 3 4 5 6 $(SEVEN) $(EIGHT)
canon_map_bins := $(foreach k,$(K), canon_maps/canon_map$(k).bin)
perm_map_bins := $(foreach k,$(K), canon_maps/perm_map$(k).bin)
canon_map_txts := $(foreach k,$(K), canon_maps/canon_map$(k).txt)
canon_list_txts := $(foreach k,$(K), canon_maps/canon_list$(k).txt)
canon_ordinal_to_signature_txts := $(foreach k,$(K), canon_maps/canon-ordinal-to-signature$(k).txt)
orbit_map_txts := $(foreach k,$(K), canon_maps/orbit_map$(k).txt)
canon_map_files := $(canon_map_bins) $(perm_map_bins) $(canon_map_txts) $(canon_list_txts) $(canon_ordinal_to_signature_txts) $(orbit_map_txts)

ehd_txts := $(foreach k,$(K), canon_maps/EdgeHammingDistance$(k).txt)
alpha_nbe_txts := $(foreach k, $(K), canon_maps/alpha_list_nbe$(k).txt)
alpha_mcmc_txts := $(foreach k, $(K), canon_maps/alpha_list_mcmc$(k).txt)
subcanon_txts := canon_maps/subcanon_map4-3.txt canon_maps/subcanon_map5-4.txt canon_maps/subcanon_map6-5.txt canon_maps/subcanon_map7-6.txt $(if $(EIGHT),canon_maps/subcanon_map8-7.txt) $(if $(SEVEN),canon_maps/subcanon_map7-6.txt)
magic_table_txts := $(foreach k,$(K), orca_jesse_blant_table/UpperToLower$(k).txt)

base: .firsttime $(LIBWAYNE_HOME)/made blant $(canon_map_files) magic_table

.firsttime:
	@echo "This may take 30-60 minutes if EIGHT is not commented out in the Makefile"
	@echo "(You will only see this message once. Pausing 10 seconds...)"
	@sleep 10
	@touch .firsttime

most: base $(alpha_nbe_txts) $(alpha_mcmc_txts) $(ehd_txts) Draw subcanon_maps

tests: test_maps test_blant

all: most tests

canon_maps: $(LIBWAYNE_HOME)/made $(canon_map_files) subcanon_maps

.PHONY: all most test_blant test_maps realclean clean_canon_maps

### Executables ###

fast-canon-map: $(LIBWAYNE_HOME)/made $(SRCDIR)/fast-canon-map.c | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CC) '-std=c99' -O3 -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/fast-canon-map.c $(LIBWAYNE_OPTS)

slow-canon-maps: $(LIBWAYNE_HOME)/made $(SRCDIR)/slow-canon-maps.c | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CC) -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/slow-canon-maps.c $(LIBWAYNE_OPTS)

make-orbit-maps: $(LIBWAYNE_HOME)/made $(SRCDIR)/make-orbit-maps.c | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CC) -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/make-orbit-maps.c $(LIBWAYNE_OPTS)

blant: $(LIBWAYNE_HOME)/made $(OBJS) $(OBJDIR)/convert.o $(OBJDIR)/libblant.o | $(LIBWAYNE_HOME)/MT19937/mt19937.o
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(OBJS) $(OBJDIR)/convert.o $(LIBWAYNE_OPTS) $(LIBWAYNE_HOME)/MT19937/mt19937.o

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(LIBWAYNE_OPTS)

synthetic: $(LIBWAYNE_HOME)/made $(SRCDIR)/synthetic.c $(SRCDIR)/syntheticDS.h $(SRCDIR)/syntheticDS.c | $(OBJDIR)/libblant.o
	$(CC) -c $(SRCDIR)/syntheticDS.c $(SRCDIR)/synthetic.c $(LIBWAYNE_OPTS)
	$(CXX) -o $@ syntheticDS.o $(OBJDIR)/libblant.o synthetic.o $(LIBWAYNE_OPTS)

CC: $(LIBWAYNE_HOME)/made $(SRCDIR)/CC.c $(OBJDIR)/convert.o | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/CC.c $(OBJDIR)/convert.o $(LIBWAYNE_OPTS)

makeEHD: $(OBJDIR)/makeEHD.o
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(OBJDIR)/makeEHD.o $(LIBWAYNE_OPTS)

compute-alphas-NBE: $(LIBWAYNE_HOME)/made $(SRCDIR)/compute-alphas-NBE.c | $(OBJDIR)/libblant.o
	$(CC) -Wall -O3 -o $@ $(SRCDIR)/compute-alphas-NBE.c $(OBJDIR)/libblant.o $(LIBWAYNE_OPTS)

compute-alphas-MCMC: $(LIBWAYNE_HOME)/made $(SRCDIR)/compute-alphas-MCMC.c | $(OBJDIR)/libblant.o
	$(CC) -Wall -O3 -o $@ $(SRCDIR)/compute-alphas-MCMC.c $(OBJDIR)/libblant.o $(LIBWAYNE_OPTS)

Draw: Draw/graphette2dot

Draw/graphette2dot: $(LIBWAYNE_HOME)/made Draw/DrawGraphette.cpp Draw/Graphette.cpp Draw/Graphette.h Draw/graphette2dotutils.cpp Draw/graphette2dotutils.h  | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CXX) -std=c++11 Draw/DrawGraphette.cpp Draw/graphette2dotutils.cpp Draw/Graphette.cpp $(OBJDIR)/libblant.o -o $@ $(LIBWAYNE_OPTS)

make-subcanon-maps: $(LIBWAYNE_HOME)/made $(SRCDIR)/make-subcanon-maps.c | $(OBJDIR)/libblant.o
	$(CC) -Wall -o $@ $(SRCDIR)/make-subcanon-maps.c $(OBJDIR)/libblant.o $(LIBWAYNE_OPTS)

make-orca-jesse-blant-table: $(LIBWAYNE_HOME)/made $(SRCDIR)/magictable.cpp | $(OBJDIR)/libblant.o
	$(CXX) -std=c++11 -Wall -o $@ $(SRCDIR)/magictable.cpp $(OBJDIR)/libblant.o $(LIBWAYNE_OPTS)

### Object Files/Prereqs ###

$(OBJDIR)/convert.o: $(SRCDIR)/convert.cpp
	@mkdir -p $(dir $@)
	$(CXX) -std=c++11 -c $(SRCDIR)/convert.cpp -o $@

$(LIBWAYNE_HOME)/MT19937/mt19937.o: $(LIBWAYNE_HOME)/made
	cd $(LIBWAYNE_HOME)/MT19937 && $(MAKE)

$(OBJDIR)/libblant.o: $(LIBWAYNE_HOME)/made $(SRCDIR)/libblant.c
	@mkdir -p $(dir $@)
	$(CC) -c $(SRCDIR)/libblant.c $(LIBWAYNE_OPTS) -o $@

$(OBJDIR)/makeEHD.o: $(LIBWAYNE_HOME)/made $(SRCDIR)/makeEHD.c | $(OBJDIR)/libblant.o
	@mkdir -p $(dir $@)
	$(CC) -c $(SRCDIR)/makeEHD.c $(LIBWAYNE_OPTS) -o $@

$(LIBWAYNE_HOME)/made:
	cd $(LIBWAYNE_HOME) && $(MAKE) all

### Generated File Recipes ###

canon_maps/orbit_map%.txt: make-orbit-maps | canon_maps/canon_list%.txt
	./make-orbit-maps $* > $@

canon_maps/canon_map%.bin canon_maps/perm_map%.bin: $(LIBWAYNE_HOME)/made $(SRCDIR)/create-bin-data.c | $(OBJDIR)/libblant.o $(SRCDIR)/blant.h canon_maps/canon_list%.txt canon_maps/canon_map%.txt
	$(CC) '-std=c99' "-Dkk=$*" "-DkString=\"$*\"" -o create-bin-data$* $(SRCDIR)/libblant.c $(SRCDIR)/create-bin-data.c $(LIBWAYNE_OPTS)
	./create-bin-data$*
	/bin/rm -f create-bin-data$*

canon_maps/canon_map%.txt canon_maps/canon_list%.txt canon_maps/canon-ordinal-to-signature%.txt: fast-canon-map 
	mkdir -p canon_maps
	./fast-canon-map $* | cut -f2- | tee canon_maps/canon_map$*.txt | awk '!seen[$$1]{seen[$$1]=1;numEdges[n]=$$4;connected[n]=$$3;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d %d %d\n", map[i],connected[i],numEdges[i]}' | tee canon_maps/canon_list$*.txt | awk 'NR>1{print NR-2, $$1}' > canon_maps/canon-ordinal-to-signature$*.txt

canon_maps/EdgeHammingDistance%.txt: makeEHD | canon_maps/canon_list%.txt canon_maps/canon_map%.bin
	@if [ -f canon_maps.correct/EdgeHammingDistance$*.txt.xz ]; then echo "EdgeHammingDistance8.txt takes weeks to generate, and 7 can't be done on a 32-bit machine; uncompressing instead"; unxz < canon_maps.correct/EdgeHammingDistance$*.txt.xz > $@ && touch $@; else ./makeEHD $* > $@; fi

canon_maps/alpha_list_nbe%.txt: compute-alphas-NBE canon_maps/canon_list%.txt
	./compute-alphas-NBE $* > $@

canon_maps/alpha_list_mcmc%.txt: compute-alphas-MCMC | canon_maps/canon_list%.txt
	@if [ -f canon_maps.correct/alpha_list_mcmc$*.txt ]; then echo "computing MCMC alphas for k=$* takes days, so just copy it"; cp -p canon_maps.correct/alpha_list_mcmc$*.txt canon_maps/ && touch $@; else ./compute-alphas-MCMC $* > $@; fi

.INTERMEDIATE: .created-magic-tables .created-subcanon-maps
subcanon_maps: $(subcanon_txts) ;
$(subcanon_txts): .created-subcanon-maps
.created-subcanon-maps: make-subcanon-maps | $(canon_list_txts) $(canon_map_bins)
	# only do it for k > 3 since it's 4-3, 5-4, etc.
	for k in $(K); do if [ $$k -gt 3 ]; then ./make-subcanon-maps $$k > canon_maps/subcanon_map$$k-$$(($$k-1)).txt; fi; done
		
magic_table: $(magic_table_txts) ;  	
$(magic_table_txts): .created-magic-tables
.created-magic-tables: make-orca-jesse-blant-table | $(canon_list_txts) $(canon_map_bins)
	./make-orca-jesse-blant-table $(if $(EIGHT),8,$(if $(SEVEN),7,6))

### Testing ###

blant-sanity: $(LIBWAYNE_HOME)/made $(SRCDIR)/blant-sanity.c
	$(CC) -o $@ $(SRCDIR)/blant-sanity.c $(LIBWAYNE_OPTS)
	
test_blant: blant blant-sanity $(canon_map_bins) test_sanity test_freq test_GDV

test_sanity:
	# First run blant-sanity for various values of k
	for k in $(K); do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity check indexing for k=$$k; ./blant -u -s NBE -mi -n 100000 -k $$k networks/syeast.el | sort -n | ./blant-sanity $$k 100000 networks/syeast.el; fi; done

test_freq:
	# Test to see that for k=6, the most frequent 10 graphlets in syeast appear in the expected order in frequency
	# Need 10 million samples to ensure with high probability we get the same graphlets.
	# We then sort them because the top 10 are a pretty stable set but their order is not.
	# The -t option tests parallelism, attemting to run multiple threads simultaneously.
	for k in $(K); do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity checking frequency of graphlets in networks/syeast.el for "k=$$k"; ./blant -s NBE -mf -t $(CORES) -n 10000000 -k $$k networks/syeast.el | sort -nr | awk '$$1{print $$2}' | head | sort -n | diff -b - testing/syeast.top10freq.k$$k.txt; fi; done

test_GDV:
	echo 'testing Graphlet (not orbit) Degree Vectors'
	for k in $(K); do export k; /bin/echo -n "$$k: "; ./blant -u -s NBE -t $(CORES) -mg -n 10000000 -k $$k networks/syeast.el | sort -n | cut -d' ' -f2- |bash -c "paste - <(unxz < testing/syeast.gdv.k$$k.txt.xz)" | $(LIBWAYNE_HOME)/bin/hawk '{cols=NF/2;for(i=1;i<=cols;i++)if($$i>1000&&$$(cols+i)>1000)printf "%.9f\n", 1-MIN($$i,$$(cols+i))/MAX($$i,$$(cols+i))}' | $(LIBWAYNE_HOME)/bin/stats | sed -e 's/#/num/' -e 's/var.*//' | $(LIBWAYNE_HOME)/bin/named-next-col '{if(num<1000 || mean>.005*'$$k' || max>0.2 || stdDev>0.005*'$$k'){printf "BEYOND TOLERANCE:\n%s\n",$$0;exit(1);}else print $$0 }' || break; done

test_maps: blant blant-sanity
	ls canon_maps.correct/ | egrep -v '$(if $(SEVEN),,7|)$(if $(EIGHT),,8|)README|\.xz' | awk '{printf "cmp canon_maps.correct/%s canon_maps/%s\n",$$1,$$1}' | sh

### Cleaning ###

clean:
	/bin/rm -f *.[oa] blant canon-sift fast-canon-map make-orbit-maps compute-alphas-MCMC compute-alphas-NBE makeEHD make-orca-jesse-blant-table Draw/graphette2dot blant-sanity make-subcanon-maps
	/bin/rm -f $(OBJDIR)/*

realclean: clean # also clean all canonical data and libwayne
	cd $(LIBWAYNE_HOME); $(MAKE) clean
	/bin/rm -f canon_maps/* .firsttime

clean_canon_maps:
	/bin/rm -f canon_maps/*[3-7].* # don't remove 8 since it takes a few minutes to create
	/bin/rm -f orca_jesse_blant_table/UpperToLower*.txt
