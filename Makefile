LIBWAYNE=-O3 -I ./libwayne/include -L libwayne -lwayne    -lm # -static OPTIMIZED
#LIBWAYNE=-O0 -I ./libwayne/include -L libwayne -lwayne-g  -lm -ggdb # for debugging
#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-pg -lm -pg   # for profiling

all: canon_maps blant test_blant magic_table draw

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

canon_maps: libwayne/made canon_maps/canon_map7.txt blant.h test_maps subcanon_maps

test_maps:
	ls canon_maps.3-7 | fgrep -v README | awk '{printf "cmp canon_maps.3-7/%s canon_maps/%s\n",$$1,$$1}' | sh

canon_maps/canon_map7.txt: blant.h make-canon-maps libblant.c create-canon-map
	mkdir -p canon_maps
	for i in 3 4 5 6 7; do ./create-canon-map $$i | cut -f2- | tee canon_maps/canon_map$$i.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list$$i.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature$$i.txt; ./make-orbit-maps $$i > canon_maps/orbit_map$$i.txt; gcc "-Dkk=$$i" "-DkString=\"$$i\"" -o create-bin-data libblant.c create-bin-data.c $(LIBWAYNE); ./create-bin-data; done
	/bin/rm -f create-bin-data # it's not useful after this

create-canon-map: create-canon-map.c blant.h canon-sift.c libblant.c make-orbit-maps
	gcc '-std=c99' -O2 -o create-canon-map libblant.c create-canon-map.c $(LIBWAYNE)

canon_map8: blant.h libblant.c create-canon-map
	echo "Warning: this will take a few minutes"
	./create-canon-map 8 | cut -f2- | tee canon_maps/canon_map8.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list8.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature8.txt
	gcc "-Dkk=8" "-DkString=\"8\"" -o create-bin-data libblant.c create-bin-data.c $(LIBWAYNE); ./create-bin-data

make-orbit-maps: make-orbit-maps.c blant.h canon-sift.c libblant.c
	gcc -o make-orbit-maps libblant.c make-orbit-maps.c $(LIBWAYNE)

blant: libwayne/made blant.c blant.h libblant.c convert.cpp
	gcc -c libblant.c blant.c $(LIBWAYNE)
	g++ -std=c++11 -c convert.cpp
	g++ -o blant libblant.o blant.o convert.o $(LIBWAYNE)
	gcc -o blant-sanity blant-sanity.c $(LIBWAYNE)

libwayne/made:
	cd libwayne; make all

subcanon_maps: libwayne/made make-subcanon-maps.c blant.h libblant.c
	mkdir -p canon_maps
	gcc -Wall -o make-subcanon-maps make-subcanon-maps.c libblant.c $(LIBWAYNE)
	for i in 4 5 6 7 8; do if [ -f canon_maps/canon_map$$i.bin -a -f canon_maps/canon_list$$i.txt ]; then  ./make-subcanon-maps $$i > canon_maps/subcanon_map$$i-$$((i-1)).txt; fi; done;
	/bin/rm -f make-subcanon-maps # it's not useful after this

magic_table: magictable.cpp
	g++ -std=c++11 -Wall -o make-orca-jesse-blant-table magictable.cpp libblant.o $(LIBWAYNE)
	./make-orca-jesse-blant-table 7

draw: Draw/DrawGraphette.cpp Draw/graphette2dotutils.h
	g++ -std=c++11 Draw/DrawGraphette.cpp -o Draw/graphette2dot

clean:
	/bin/rm -f *.[oa] blant canon-sift create-canon-map
	/bin/rm -f canon_maps/*[3-7].* # don't remove 8 since it takes a few minutes to create
	#cd libwayne; make clean
