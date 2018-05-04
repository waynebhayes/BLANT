LIBWAYNE=-O3 -I ./libwayne/include -L libwayne -lwayne    -lm # -static OPTIMIZED
#LIBWAYNE=-O0 -I ./libwayne/include -L libwayne -lwayne-g  -lm -ggdb # for debugging
#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-pg -lm -pg   # for profiling

all: canon_maps blant test_blant magic_table

test_blant:
	# First run blant-sanity for various values of k
	for k in 3 4 5 6 7 8; do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity check indexing for k=$$k; ./blant -mi -s 100000 -k $$k networks/syeast.el | sort -n | ./blant-sanity $$k networks/syeast.el; fi; done
	# Test to see that for k=6, the most frequent 10 graphlets in syeast appear in the expected order in frequency
	# Need 10 million samples to ensure with high probability we get the same graphlets.
	# We then sort them because the top 10 are a pretty stable set but their order is not.
	# The -t option tests parallelism, attemting to run multiple threads simultaneously.
	for k in 3 4 5 6 7 8; do if [ -f canon_maps/canon_map$$k.bin ]; then echo sanity checking frequency of graphlets in networks/syeast.el for "k=$$k"; ./blant -mf -t 4 -s 10000000 -k $$k networks/syeast.el | sort -nr | awk '$$1{print $$2}' | head | sort -n | diff -b - testing/syeast.top10freq.k$$k.txt; fi; done
	echo 'testing Graphlet (not orbit) Degree Vectors'
	for k in 3 4 5 6 7 8; do export k; echo -n "$$k: "; ./blant -t 4 -mg -s 10000000 -k $$k networks/syeast.el | bash -c "paste - <(unxz < testing/syeast.gdv.k$$k.txt.xz)" | ./libwayne/bin/hawk '{cols=NF/2;for(i=1;i<=cols;i++)if($$i>1000&&$$(cols+i)>1000)printf "%.9f\n", 1-MIN($$i,$$(cols+i))/MAX($$i,$$(cols+i))}' | ./libwayne/bin/stats | sed -e 's/#/num/' -e 's/var.*//' | ./libwayne/bin/named-next-col '{if(num<1000 || mean>.005*'$$k' || max>0.2 || stdDev>0.005*'$$k'){printf "BEYOND TOLERANCE:\n%s\n",$$0;exit(1);}else print $$0 }' || break; done

canon_maps: libwayne canon_maps/canon_map6.txt blant.h test_maps subcanon_maps

test_maps:
	ls canon_maps.3-6 | fgrep -v README | awk '{printf "cmp canon_maps.3-6/%s canon_maps/%s\n",$$1,$$1}' | sh

canon_map7 canon_map8:
	@echo "Making the k=7 canon_map takes a few hours"
	@echo "Making the k=8 canon_map takes many months of CPU time."
	@echo "See the k8 directory for instructions on how to do it using parallelism."
	@echo "Or see http://www.ics.uci.edu/~wayne/blant to just download the files."

canon_maps/canon_map6.txt: blant.h make-canon-maps libblant.c create-canon-map.c
	mkdir -p canon_maps
	for i in 3 4 5 6; do ./make-canon-maps $$i 1 0 | cut -f2- | tee canon_maps/canon_map$$i.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list$$i.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature$$i.txt; ./make-orbit-maps $$i > canon_maps/orbit_map$$i.txt; gcc "-Dkk=$$i" "-DkString=\"$$i\"" -o create-bin-data libblant.c create-bin-data.c $(LIBWAYNE); ./create-bin-data; done
	#for i in 3 4 5 6 7; do gcc '-std=c99' "-Dk=$$i" -o create-canon-map create-canon-map.c; ./create-canon-map | tee canon_maps/canon_map$$i.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list$$i.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature$$i.txt; ./make-orbit-maps $$i > canon_maps/orbit_map$$i.txt; gcc "-Dkk=$$i" "-DkString=\"$$i\"" -o create-bin-data libblant.c create-bin-data.c $(LIBWAYNE); ./create-bin-data; done
	/bin/rm -f create-bin-data ./create-canon-map # it's not useful after this

canon_maps/canon_map7.txt: blant.h make-canon-maps libblant.c
	echo "Warning: this will take an hour or more depending machine speed"
	echo "You can also look in the k8 directory for the k=8 script and test it for k=7"
	./make-canon-maps 7 1 0 | cut -f2- | tee canon_maps/canon_map7.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list7.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature7.txt
	gcc "-Dk=7" "-DkString=\"7\"" -o create-bin-data libblant.c create-bin-data.c; ./create-bin-data

make-canon-maps: make-canon-maps.c blant.h canon-sift.c libblant.c make-orbit-maps
	gcc -o make-canon-maps libblant.c make-canon-maps.c $(LIBWAYNE)
	gcc -o canon-sift libblant.c canon-sift.c  $(LIBWAYNE)

make-orbit-maps: make-orbit-maps.c blant.h canon-sift.c libblant.c
	gcc -o make-orbit-maps libblant.c make-orbit-maps.c $(LIBWAYNE)

blant: libwayne blant.c blant.h libblant.c convert.cpp
	gcc -c libblant.c blant.c $(LIBWAYNE)
	g++ -std=c++11 -c convert.cpp
	g++ -o blant libblant.o blant.o convert.o $(LIBWAYNE)
	gcc -o blant-sanity blant-sanity.c $(LIBWAYNE)

libwayne: libwayne/libwayne.a libwayne/libwayne-g.a

libwayne/libwayne.a:
	cd libwayne; make opt_clean; make opt

libwayne/libwayne-g.a:
	cd libwayne; make debug_clean; make debug

subcanon_maps: libwayne make-subcanon-maps.c blant.h libblant.c
	mkdir -p canon_maps
	gcc -Wall -o make-subcanon-maps make-subcanon-maps.c libblant.c $(LIBWAYNE)
	for i in 4 5 6 7 8; do if [ -f canon_maps/canon_map$$i.bin -a -f canon_maps/canon_list$$i.txt ]; then  ./make-subcanon-maps $$i > canon_maps/subcanon_map$$i-$$((i-1)).txt; fi; done;
	/bin/rm -f make-subcanon-maps # it's not useful after this

magic_table: magictable.cpp
	g++ -std=c++11 -Wall -o make-orca-jesse-blant-table magictable.cpp libblant.o $(LIBWAYNE)
	./make-orca-jesse-blant-table 7

clean:
	/bin/rm -f *.[oa] blant make-canon-maps canon-sift
	/bin/rm -f canon_maps/*[3-6].* # don't remove 7 or 8 unless you REALLY want to since they take long to create
	#cd libwayne; make clean
