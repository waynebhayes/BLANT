# Some architectures, eg CYGWIN 32-bit and MacOS("Darwin") need an 80MB stack.
STACKSIZE=$(shell arch | awk '/CYGWIN/{print "-Wl,--stack,83886080"}/Darwin/{print "-Wl,-stack_size -Wl,0x5000000"}')
LIBWAYNE=-O3 -I ./libwayne/include -L libwayne -lwayne -lm $(STACKSIZE) # -static OPTIMIZED
#LIBWAYNE=-O0 -I ./libwayne/include -L libwayne -lwayne-g  -lm -ggdb $(STACKSIZE) # for debugging
#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-pg -lm -pg   # for profiling

CC=gcc
CXX=g++

### Generated File Lists ###
EIGHT := #8
K := 3 4 5 6 7 $(EIGHT)
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
subcanon_txts := canon_maps/subcanon_map4-3.txt canon_maps/subcanon_map5-4.txt canon_maps/subcanon_map6-5.txt canon_maps/subcanon_map7-6.txt $(if $(EIGHT),canon_maps/subcanon_map8-7.txt)
magic_table_txts := $(foreach k,$(K), orca_jesse_blant_table/UpperToLower$(k).txt)

most: libwayne/made blant $(canon_map_files) magic_table $(alpha_nbe_txts) $(alpha_mcmc_txts) $(ehd_txts) Draw subcanon_maps

all: most test_maps test_blant

.PHONY: all most test_blant test_maps realclean clean_canon_maps

### Executables ###

fast-canon-map: libwayne/made fast-canon-map.c | blant.h libblant.o
	$(CC) '-std=c99' -O3 -o $@ libblant.o fast-canon-map.c $(LIBWAYNE)

slow-canon-maps: libwayne/made slow-canon-maps.c | blant.h libblant.o
	$(CC) -o $@ libblant.o slow-canon-maps.c $(LIBWAYNE)

make-orbit-maps: libwayne/made make-orbit-maps.c | blant.h libblant.o
	$(CC) -o $@ libblant.o make-orbit-maps.c $(LIBWAYNE)

blant: libwayne/made blant.c blant.h convert.o libblant.o | libwayne/MT19937/mt19937.o
	$(CC) -c blant.c $(LIBWAYNE)
	$(CXX) -o $@ libblant.o blant.o convert.o $(LIBWAYNE) libwayne/MT19937/mt19937.o

synthetic: libwayne/made synthetic.c syntheticDS.h syntheticDS.c | libblant.o
	$(CC) -c syntheticDS.c synthetic.c $(LIBWAYNE)
	$(CXX) -o $@ syntheticDS.o libblant.o synthetic.o $(LIBWAYNE)

CC: libwayne/made CC.c convert.o | blant.h libblant.o
	$(CXX) -o $@ libblant.o CC.c convert.o $(LIBWAYNE)

makeEHD: libwayne/made makeEHD.c | libblant.o
	$(CC) -c makeEHD.c $(LIBWAYNE)
	$(CXX) -o $@ libblant.o makeEHD.o $(LIBWAYNE)

compute-alphas-NBE: libwayne/made compute-alphas-NBE.c | libblant.o
	$(CC) -Wall -O3 -o $@ compute-alphas-NBE.c libblant.o $(LIBWAYNE)

compute-alphas-MCMC: libwayne/made compute-alphas-MCMC.c | libblant.o
	$(CC) -Wall -O3 -o $@ compute-alphas-MCMC.c libblant.o $(LIBWAYNE)

Draw: Draw/graphette2dot

Draw/graphette2dot: libwayne/made Draw/DrawGraphette.cpp Draw/Graphette.cpp Draw/Graphette.h Draw/graphette2dotutils.cpp Draw/graphette2dotutils.h  | blant.h libblant.o
	$(CXX) -std=c++11 Draw/DrawGraphette.cpp Draw/graphette2dotutils.cpp Draw/Graphette.cpp libblant.o -o $@ $(LIBWAYNE)

make-subcanon-maps: libwayne/made make-subcanon-maps.c | libblant.o
	$(CC) -Wall -o $@ make-subcanon-maps.c libblant.o $(LIBWAYNE)

make-orca-jesse-blant-table: libwayne/made magictable.cpp | libblant.o
	$(CXX) -std=c++11 -Wall -o $@ magictable.cpp libblant.o $(LIBWAYNE)

### Object Files/Prereqs ###

convert.o: convert.cpp
	$(CXX) -std=c++11 -c convert.cpp

libwayne/MT19937/mt19937.o: libwayne/made
	cd libwayne/MT19937 && $(MAKE)

libblant.o: libwayne/made libblant.c
	$(CC) -c libblant.c $(LIBWAYNE)

libwayne/made:
	cd libwayne && $(MAKE) all

### Generated File Recipes ###

canon_maps/orbit_map%.txt: make-orbit-maps | canon_maps/canon_list%.txt
	./make-orbit-maps $* > $@

canon_maps/canon_map%.bin canon_maps/perm_map%.bin: libwayne/made create-bin-data.c | libblant.o blant.h canon_maps/canon_list%.txt canon_maps/canon_map%.txt
	$(CC) '-std=c99' "-Dkk=$*" "-DkString=\"$*\"" -o create-bin-data$* libblant.c create-bin-data.c $(LIBWAYNE)
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
	./make-orca-jesse-blant-table $(if $(EIGHT),8,7)

### Testing ###

blant-sanity: libwayne/made blant-sanity.c
	$(CC) -o $@ blant-sanity.c $(LIBWAYNE)
	
test_blant: blant blant-sanity $(canon_map_bins) test_sanity test_freq test_GDV regression

test_sanity:
	# First run blant-sanity for various values of k
	for k in $(K); do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity check indexing for k=$$k; ./blant -s NBE -mi -n 100000 -k $$k networks/syeast.el | sort -n | ./blant-sanity $$k 100000 networks/syeast.el; fi; done

test_freq:
	# Test to see that for k=6, the most frequent 10 graphlets in syeast appear in the expected order in frequency
	# Need 10 million samples to ensure with high probability we get the same graphlets.
	# We then sort them because the top 10 are a pretty stable set but their order is not.
	# The -t option tests parallelism, attemting to run multiple threads simultaneously.
	for k in $(K); do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity checking frequency of graphlets in networks/syeast.el for "k=$$k"; ./blant -s NBE -mf -t 4 -n 10000000 -k $$k networks/syeast.el | sort -nr | awk '$$1{print $$2}' | head | sort -n | diff -b - testing/syeast.top10freq.k$$k.txt; fi; done

test_GDV:
	echo 'testing Graphlet (not orbit) Degree Vectors'
	for k in $(K); do export k; /bin/echo -n "$$k: "; ./blant -s NBE -t 4 -mg -n 10000000 -k $$k networks/syeast.el | sort -n | cut -d' ' -f2- |bash -c "paste - <(unxz < testing/syeast.gdv.k$$k.txt.xz)" | ./libwayne/bin/hawk '{cols=NF/2;for(i=1;i<=cols;i++)if($$i>1000&&$$(cols+i)>1000)printf "%.9f\n", 1-MIN($$i,$$(cols+i))/MAX($$i,$$(cols+i))}' | ./libwayne/bin/stats | sed -e 's/#/num/' -e 's/var.*//' | ./libwayne/bin/named-next-col '{if(num<1000 || mean>.005*'$$k' || max>0.2 || stdDev>0.005*'$$k'){printf "BEYOND TOLERANCE:\n%s\n",$$0;exit(1);}else print $$0 }' || break; done

test_maps: blant blant-sanity
	ls canon_maps.correct/ | egrep -v '$(if $(EIGHT),,8|)README|\.xz' | awk '{printf "cmp canon_maps.correct/%s canon_maps/%s\n",$$1,$$1}' | sh

regression:
	./regression-test-all.sh

### Cleaning ###

clean:
	/bin/rm -f *.[oa] blant canon-sift fast-canon-map make-orbit-maps compute-alphas-MCMC compute-alphas-NBE makeEHD make-orca-jesse-blant-table Draw/graphette2dot blant-sanity make-subcanon-maps

realclean: clean # also clean all canonical data and libwayne
	cd libwayne; $(MAKE) clean
	/bin/rm -f canon_maps/*

clean_canon_maps:
	/bin/rm -f canon_maps/*[3-7].* # don't remove 8 since it takes a few minutes to create
	/bin/rm -f orca_jesse_blant_table/UpperToLower*.txt
