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
