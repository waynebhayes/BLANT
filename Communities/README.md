Communities

To compile first do 

export LIBWAYNE_HOME={your path to BLANT starting with ~}/libwayne

Ex: 
export LIBWAYNE_HOME=~/BLANT/libwayne

Then in the communities directory

../libwayne/bin/wgcc ./communities.c

OR

make


To run 

./a.out {Graph to run on} {Optional partition to start from}
