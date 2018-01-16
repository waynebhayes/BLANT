#LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne-g -lm -ggdb  # for debugging
LIBWAYNE=-I ./libwayne/include -L libwayne -lwayne -lm # for debugging

all: canon_maps blant test_blant

test_blant:
	# Test to see that for k=6, the most frequent 10 graphlets in syeast appear in the expected order in frequency
	# Need 10 million samples to ensure with high probability we get the same graphlets.
	# We then sort them because the top 10 are a pretty stable set but their order is not.
	# The -4 also tests parallelism, attemting to run 4 threads simultaneously.
	./blant -4 6 10000000 syeast.el | awk '{print $$1}' | sort | uniq -c | sort -nr | head | awk '{print $$2}' | sort -n | diff - blant.k6.syeast.out

canon_maps: libwayne canon_maps/canon_map6.txt blant.h test_maps

test_maps:
	ls canon_maps.3-6 | awk '{printf "cmp canon_maps.3-6/%s canon_maps/%s\n",$$1,$$1}' | sh

canon_map7: libwayne canon_maps/canon_map7.txt blant.h

canon_map8:
	echo "Making the k=8 canon_map takes many months of CPU time."
	echo "See the k8 directory for instructions on how to do it using parallelism."
	echo "Or see http://www.ics.uci.edu/~wayne/blant to just download the files."

canon_maps/canon_map6.txt: blant.h make-canon-maps
	mkdir -p canon_maps
	for i in 3 4 5 6; do ./make-canon-maps $$i 1 0 | tee canon_maps/canon_map$$i.txt | awk '!seen[$$1]{seen[$$1]=1;map[n++]=$$1}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list$$i.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature$$i.txt; gcc -O2 "-Dk=$$i" "-DkString=\"$$i\"" -o create-bin-data create-bin-data.c; ./create-bin-data; done
	/bin/rm -f create-bin-data # it's not useful after this

canon_maps/canon_map7.txt: blant.h make-canon-maps
	echo "Warning: this will take an hour or more depending machine speed"
	./make-canon-maps 7 1 0 | tee canon_maps/canon_map7.txt | awk '!seen[$$2]{seen[$$2]=1;map[n++]=$$2}END{print n;for(i=0;i<n;i++)printf "%d ", map[i]; print ""}' | tee canon_maps/canon_list7.txt | awk 'NR==2{for(i=1;i<=NF;i++) print i-1, $$i}' > canon_maps/canon-ordinal-to-signature7.txt

make-canon-maps: make-canon-maps.c blant.h
	gcc -O2 -o make-canon-maps make-canon-maps.c $(LIBWAYNE)
	gcc -O2 -o canon-sift canon-sift.c $(LIBWAYNE)

blant: libwayne blant.c blant.h
	gcc -O2 -o blant blant.c $(LIBWAYNE)

libwayne: libwayne/libwayne.a

libwayne/libwayne.a:
	cd libwayne; make opt

clean:
	/bin/rm -f *.[oa] blant make-canon-maps canon-sift
	/bin/rm -f canon_maps/*[3-6].* # don't remove 7 or 8 unless you REALLY want to since they take long to create
	# cd libwayne; make clean
