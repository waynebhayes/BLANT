#!/bin/bash

k=5
G=$k.el

for i in {3..5}
do
	ARGS="-k$i -mi -n3 regression-tests/selfLoops/$G"
	./blant $ARGS 1> /dev/null
	echo "Testing blant with self-loops for k=$i"
done

echo "Done testing self-loops"