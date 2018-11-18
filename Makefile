LIBWAYNE=-O3 -I ./libwayne/include -L libwayne -lwayne    -lm # -static OPTIMIZED
#LIBWAYNE=-O0 -I ./libwayne/include -L libwayne -lwayne-g  -lm -ggdb # for debugging
#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-pg -lm -pg   # for profiling

most: libwayne/made blant canon_maps compute-alphas-MCMC compute-alphas-NBE magic_table draw

all: most canon_map8 test_blant test_maps

test_blant:
	# First run blant-sanity for various values of k
	for k in 3 4 5 6 7 8; do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity check indexing for k=$$k; ./blant -s NBE -mi -n 100000 -k $$k networks/syeast.el | sort -n | ./blant-sanity $$k networks/syeast.el; fi; done
	# Test to see that for k=6, the most frequent 10 graphlets in syeast appear in the expected order in frequency
	# Need 10 million samples to ensure with high probability we get the same graphlets.
	# We then sort them because the top 10 are a pretty stable set but their order is not.
	# The -t option tests parallelism, attemting to run multiple threads simultaneously.
	for k in 3 4 5 6 7 8; do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity checking frequency of graphlets in networks/syeast.el for "k=$$k"; ./blant -s NBE -mf -t 4 -n 10000000 -k $$k networks/syeast.el | sort -nr | awk '$$1{print $$2}' | head | sort -n | diff -b - testing/syeast.top10freq.k$$k.txt; fi; done
	echo 'testing Graphlet (not orbit) Degree Vectors'
	for k in 3 4 5 6 7 8; do export k; echo -n "$$k: "; ./blant -s NBE -t 4 -mg -n 10000000 -k $$k networks/syeast.el | bash -c "paste - <(unxz < testing/syeast.gdv.k$$k.txt.xz)" | ./libwayne/bin/hawk '{cols=NF/2;for(i=1;i<=cols;i++)if($$i>1000&&$$(cols+i)>1000)printf "%.9f\n", 1-MIN($$i,$$(cols+i))/MAX($$i,$$(cols+i))}' | ./libwayne/bin/stats | sed -e 's/#/num/' -e 's/var.*//' | ./libwayne/bin/named-next-col '{if(num<1000 || mean>.005*'$$k' || max>0.2 || stdDev>0.005*'$$k'){printf "BEYOND TOLERANCE:\n%s\n",$$0;exit(1);}else print $$0 }' || break; done

test_maps:
	ls canon_maps/ | egrep -v 'README|graphlet_list' | awk '{printf "cmp canon_maps.3-7/%s canon_maps/%s\n",$$1,$$1}' | sh

canon_maps: blant.h fast-canon-map libblant.c libwayne/made subcanon_maps
	mkdir -p canon_maps
	for k in 3 4 5 6 7 $$EIGHT; do if [ ! -f canon_maps/canon_map$$k.bin ]; then ./fast-canon-map $$k | cut -f2- | tee canon_maps/canon_map$$k.txt | awk '!seen[$$1]{seen[$$1]=1;numEdges[n]=$$4;connected[n]=$$3;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d %d %d\n", map[i],connected[i],numEdges[i]}' | tee canon_maps/canon_list$$k.txt | awk 'NR>1{print NR-2, $$1}' > canon_maps/canon-ordinal-to-signature$$k.txt; ./make-orbit-maps $$k > canon_maps/orbit_map$$k.txt; gcc "-Dkk=$$k" "-DkString=\"$$k\"" -o create-bin-data$$k libblant.c create-bin-data.c $(LIBWAYNE); ./create-bin-data$$k; /bin/rm -f create-bin-data$$k; fi; done

canon_map8:
	echo "Warning: this will take a few minutes to a few hours, depending upon your machine"
	make 'EIGHT=8' canon_maps
	make 'EIGHT=8' subcanon_maps

fast-canon-map: fast-canon-map.c blant.h canon-sift.c libblant.c make-orbit-maps
	gcc '-std=c99' -O2 -o fast-canon-map libblant.c fast-canon-map.c $(LIBWAYNE)

slow-canon-maps: slow-canon-maps.c blant.h libblant.c
	gcc -o slow-canon-maps libblant.c slow-canon-maps.c $(LIBWAYNE)

make-orbit-maps: make-orbit-maps.c blant.h canon-sift.c libblant.c
	gcc -o make-orbit-maps libblant.c make-orbit-maps.c $(LIBWAYNE)

blant: libwayne/made blant.c blant.h libblant.c convert.cpp libwayne/MT19937/mt19937.o
	gcc -c libblant.c blant.c $(LIBWAYNE)
	g++ -std=c++11 -c convert.cpp
	g++ -o blant libblant.o blant.o convert.o $(LIBWAYNE) libwayne/MT19937/mt19937.o
	gcc -o blant-sanity blant-sanity.c $(LIBWAYNE)

synthetic: libwayne/made synthetic.c syntheticDS.h syntheticDS.c blant.h libblant.c
	gcc -c syntheticDS.c libblant.c synthetic.c $(LIBWAYNE)
	g++ -o synthetic syntheticDS.o libblant.o synthetic.o $(LIBWAYNE)

libwayne/MT19937/mt19937.o:
	(cd libwayne/MT19937 && make)

compute-alphas-NBE: #compute-alphas-NBE.c
	gcc -Wall -O2 -o compute-alphas-NBE compute-alphas-NBE.c libblant.o $(LIBWAYNE)
	for k in 3 4 5 6 7; do if [ -f canon_maps/canon_list$$k.txt -a ! -f canon_maps/alpha_list_nbe$$k.txt ]; then ./compute-alphas-NBE $$k; fi; done
	
compute-alphas-MCMC: #libblant.c blant.h compute-alphas-MCMC.c canon_maps/alpha_list7.txt
	gcc -Wall -O2 -o compute-alphas-MCMC compute-alphas-MCMC.c libblant.o $(LIBWAYNE)
	if [ "$$EIGHT" = 8 ]; then SEVEN=7; fi; for k in 3 4 5 6 $$SEVEN; do if [ -f canon_maps/canon_list$$k.txt -a ! -f canon_maps/alpha_list$$k.txt ]; then ./compute-alphas-MCMC $$k; fi; done

CC: libwayne/made CC.c blant.h libblant.c convert.cpp
	gcc -c libblant.c CC.c $(LIBWAYNE)
	g++ -std=c++11 -c convert.cpp
	g++ -o CC libblant.o CC.o convert.o $(LIBWAYNE)

libwayne/made:
	cd libwayne; make all

subcanon_maps: #libwayne/made make-subcanon-maps.c blant.h libblant.c
	mkdir -p canon_maps
	gcc -Wall -o make-subcanon-maps make-subcanon-maps.c libblant.c $(LIBWAYNE)
	for k in 4 5 6 7 $$EIGHT; do if [ -f canon_maps/canon_map$$k.bin -a -f canon_maps/canon_list$$k.txt ]; then  ./make-subcanon-maps $$k > canon_maps/subcanon_map$$k-$$((k-1)).txt; fi; done;
	/bin/rm -f make-subcanon-maps # it's not useful after this

magic_table: magictable.cpp
	g++ -std=c++11 -Wall -o make-orca-jesse-blant-table magictable.cpp libblant.o $(LIBWAYNE)
	./make-orca-jesse-blant-table 7

draw: Draw/DrawGraphette.cpp Draw/graphette2dotutils.h
	g++ -std=c++11 Draw/DrawGraphette.cpp -o Draw/graphette2dot

clean:
	/bin/rm -f *.[oa] blant canon-sift fast-canon-map compute-alphas-MCMC compute-alphas-NBE canon_maps/*[3-7]*

realclean: clean
	cd libwayne; make clean
	/bin/rm -f canon_maps/*

clean_canon_maps:
	/bin/rm -f canon_maps/*[3-7].* # don't remove 8 since it takes a few minutes to create
	/bin/rm -f orca_jesse_blant_table/UpperToLower*.txt
