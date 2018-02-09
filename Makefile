 LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne    -lm # -static OPTIMIZED
#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-g  -lm -ggdb # for debugging
#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-pg -lm -pg   # for profiling

all: canon_maps blant test_blant

test_blant:
	# First run blant-sanity for various values of k
	for k in 3 4 5 6 7 8; do if [ -f canon_maps/canon_map$$k.bin ]; then echo running sanity check for k=$$k; ./blant $$k 100000 syeast.el | sort -n | ./blant-sanity $$k syeast.el; fi; done
	# Test to see that for k=6, the most frequent 10 graphlets in syeast appear in the expected order in frequency
	# Need 10 million samples to ensure with high probability we get the same graphlets.
	# We then sort them because the top 10 are a pretty stable set but their order is not.
	# The -2 also tests parallelism, attemting to run 2 threads simultaneously.
	# NOTE THIS WILL FAIL UNLESS YOU SET BOTH LOWER_TRIANGLE AND PERMS_CAN2NON TO 1 IN blant.h.
	for k in 3 4 5 6 7 8; do if [ -f canon_maps/canon_map$$k.bin ]; then echo checking frequency of graphlets in syeast.el for "k=$$k"; ./blant -2 $$k 10000000 syeast.el | awk '{++count[$$1]}END{for(i in count) print count[i],i}' | sort -nr | head | awk '{print $$2}' | sort -n | diff -b - blant.k$$k.syeast.out; fi; done

canon_maps: libwayne canon_maps/canon_map6.txt blant.h test_maps subcanon_maps

test_maps:
	ls canon_maps.3-6 | fgrep -v README | awk '{printf "cmp canon_maps.3-6/%s canon_maps/%s\n",$$1,$$1}' | sh

canon_map7: libwayne canon_maps/canon_map7.txt blant.h

canon_map8:
	echo "Making the k=8 canon_map takes many months of CPU time."
	echo "See the k8 directory for instructions on how to do it using parallelism."
	echo "Or see http://www.ics.uci.edu/~wayne/blant to just download the files."

canon_maps/canon_map6.txt: blant.h make-canon-maps
	mkdir -p canon_maps
	for i in 3 4 5 6; do ./make-canon-maps $$i 1 0 | cut -f2- | tee canon_maps/canon_map$$i.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list$$i.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature$$i.txt; gcc -O2 "-Dk=$$i" "-DkString=\"$$i\"" -o create-bin-data create-bin-data.c; ./create-bin-data; done
	/bin/rm -f create-bin-data # it's not useful after this

canon_maps/canon_map7.txt: blant.h make-canon-maps
	echo "Warning: this will take an hour or more depending machine speed"
	echo "You can also look in the k8 directory for the k=8 script and test it for k=7"
	./make-canon-maps 7 1 0 | cut -f2- | tee canon_maps/canon_map7.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list7.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature7.txt
	gcc -O2 "-Dk=7" "-DkString=\"7\"" -o create-bin-data create-bin-data.c; ./create-bin-data

make-canon-maps: make-canon-maps.c blant.h canon-sift.c
	gcc -O2 -o make-canon-maps make-canon-maps.c $(LIBWAYNE)
	gcc -O2 -o canon-sift canon-sift.c $(LIBWAYNE)

blant: libwayne blant.c blant.h
	gcc -O2 -o blant blant.c $(LIBWAYNE)
	gcc -O2 -o blant-sanity blant-sanity.c $(LIBWAYNE)

libwayne: libwayne/libwayne.a

libwayne/libwayne.a:
	cd libwayne; make opt

subcanon_maps: libwayne make-subcanon-maps.c blant.h
	mkdir -p canon_maps
	gcc -O2 -Wall -o make-subcanon-maps make-subcanon-maps.c $(LIBWAYNE)
	for i in  4 5 6 7 8; do if [ -f canon_maps/canon_map$$((i-1)).bin -a -f canon_maps/canon_list$$i.txt ]; then  ./make-subcanon-maps $$i > canon_maps/subcanon_map$$i-$$((i-1)).txt; fi; done;
	/bin/rm -f make-subcanon-maps # it's not useful after this


clean:
	/bin/rm -f *.[oa] blant make-canon-maps canon-sift
	/bin/rm -f canon_maps/*[3-6].* # don't remove 7 or 8 unless you REALLY want to since they take long to create
	# cd libwayne; make clean
