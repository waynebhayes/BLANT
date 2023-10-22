ESCAPE_HOME := .

OBJECTS := Graph.o GraphIO.o 

TARGETS := libescape.a

all: tests

libescape.a : $(OBJECTS)
	ar cruv $@ $^

include common.mk

tests: libescape.a
	$(MAKE) -C tests

cleantests:
	$(MAKE) -C tests clean cleandep

.PHONY: tests


