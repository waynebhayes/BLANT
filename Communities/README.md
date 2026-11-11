Communities

To compile first do 

export LIBWAYNE_HOME={your path to BLANT starting with ~}/libwayne

Ex: 
export LIBWAYNE_HOME=~/BLANT/libwayne

Then in the communities directory

make

in your command line to compile

To run 

./a.out {Graph to run on} 

Ex:
./a.out ./communities.in

OPTIONS
-s Choose a scoring function
-p Provide a partition to start from 
-r Random seed to use (default: a fancy seed, printed as "Random seed N"), so a run can be reproduced exactly
-t Run this many independent anneals in parallel (one process each, seeds N, N+1, ...) and keep the best;
   0 means one per CPU. Prints the best run's output followed by a "Best of" line listing every run's score.
   Wall time is about that of a single run, as long as there are enough cores.

Ex:
./a.out ./communities.in -t 8 -r 100
