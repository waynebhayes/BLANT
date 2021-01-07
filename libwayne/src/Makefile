# Makefile for "Wayne's Little DSA Library"
# (DSA == Data Structures and Algorithms)
# Feel free to change, modify, or burn these sources, but if you modify them
# please don't distribute your changes without a clear indication that you've done so.
# If you think your change is spiffy, send it to me and maybe I'll include it
# 
# Wayne Hayes, wayne@csri.utoronto.ca (preffered), or wayne@csri.toronto.edu

OBJS=llfile.o stream48.o longlong.o sets.o smallgraph-transitive.o misc.o dverk.o rkd78.o lsode.o ddriv2.o bsode.o ldbsode.o rk4.o rk4s.o rk12.o rk23.o stack.o event.o heap.o linked-list.o stats.o queue.o compressedInt.o Oalloc.o variable_leapfrog.o leapfrog.o htree.o bintree.o eigen.o mem-debug.o smallgraph.o tinygraph.o graph.o combin.o matvec.o sorts.o heun_euler.o multisets.o dynarray.o hashmap.o #qrkd78.o iqrkd78.o

INCLUDE=-I../include
#LIB=$(HOME)/lib/libwayne.a

CFLAGS=$(INCLUDE)

$(LIBOUT): $(OBJS)
	ar rcvu $(LIBOUT) $(OBJS)

event.c: ../include/event.h
