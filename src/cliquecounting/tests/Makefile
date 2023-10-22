ESCAPE_HOME := ../

TARGETS := test_cliques_turan_shadow test_cliques_color_coding test_cliques_brute_force \
            test_cliques_edge_sampling test_cliques_graft

OBJECTS := $(TARGETS:%=%.o)

include ../common.mk

INCLUDES += -I..
LDFLAGS += -L..

test_% : test_%.o ../libescape.a
	$(COMPILE_AND_LINK)
