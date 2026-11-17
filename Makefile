# Plain 'make' builds 'base': libwayne, blant, and the graphlet tables in canon_maps/.
#
# Nothing here goes by file dates alone. A program is recompiled only when the contents of
# what it is built from change (its sources, the headers they read, its flags, libwayne),
# and a table is remade only when the code that makes it changes. So a 'touch' or a
# 'git checkout' recompiles nothing, and a 'make clean' remakes no undirected table (it
# deletes the directed ones, which are quick), which matters because the k=8 tables can take
# an hour.

# Number of cores to use when invoking parallelism
# (temporarily set to 1 since I broke threading. :-( )
CORES ?= 2
ifndef PAUSE
	PAUSE := 1
endif
# When DYNAMIC_CANON_MAP=1, skip generating canon maps and all related files (canon_map/canon_list/
# orbit_map/alpha_list/subcanon_map/magic tables), and skip all tests that depend on them.
# Only libwayne and the blant binary itself get built.
DYNAMIC_CANON_MAP ?= 0
# Uncomment either of these to remove them (removing 7 implies removing 8)
MAX_K := 8
EIGHT := 8
SEVEN := 7
ifdef NO8
    MAX_K := 7
    EIGHT :=
endif
ifdef NO7
    MAX_K := 6
    SEVEN :=
    EIGHT := # can't have 8 without 7
endif

# to make the prediction version that agrees with the regression tests
ifdef PRED_REG
    PRED_REG_OPT := -DINTERNAL_DEG_WEIGHTS=0 -DDEG_ORDER_MUST_AGREE=1
endif

ifdef DEBUG
    ifdef PROFILE
	SPEED=-O0 -ggdb -pg
	LIB_OPT=-pg-g
    else
	SPEED=-O0 -ggdb
	LIB_OPT=-g
    endif
else
    ifdef PROFILE
	SPEED=-O3 -pg
	LIB_OPT=-pg
    else
	SPEED=-O3 #-DNDEBUG
	LIB_OPT= #-nd # NDEBUG
    endif
endif

ifdef STATIC
    STATIC_LINK=-static
endif

#GCC= $(shell $(CC) -v 2>&1 | awk '/gcc/{++gcc}{V=$$3}END{if(gcc && (V ~ /[0-9]\.[0-9]\.[0-9]*/))print "$(UNAME).gcc"V; else exit 1}')
#GCC_VER=$(shell echo $(UNAME) $(HOME) | awk '/Darwin/&&/Users.wayne/{V="-6"}END{if(V)print V;else{printf "using default gcc: " > "/dev/null"; exit 1}}')
GCC=gcc$(GCC_VER) $(SPEED) $(NDEBUG)
CXX=g++$(GCC_VER) $(SPEED) $(NDEBUG)

# Some systems, eg CYGWIN 32-bit and MacOS("Darwin") need a very large stack.
export LIBWAYNE_HOME=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/libwayne
UNAME=$(shell uname -a | awk '{if(/CYGWIN/){V="CYGWIN"}else if(/Darwin/){if(/arm64/)V="arm64";else V="Darwin"}else if(/Linux/){V="Linux"}}END{if(V){print V;exit}else{print "unknown OS" > "/dev/stderr"; exit 1}}')
#Old values for Cygwin 32-bit: 83886080, 0x5000000
STACKSIZE=$(shell ($(GCC) -v 2>/dev/null; uname -a) | awk '/CYGWIN/{print "-Wl,--stack,83886080"}/gcc-/{actualGCC=1}/Darwin/{print "-Wl,-stack_size -Wl,0x20000000"}')
# 32-bit machine
BITS=$(shell (uname -a) | awk '/i686/{print 32}/x86_64/||[Ll]inux{print 64}')
# GCC and CXX already hold $(SPEED) and $(NDEBUG), so they are not repeated here.
CC=$(GCC) -Wno-misleading-indentation -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow $(PG)
LIBWAYNE_COMP=-I $(LIBWAYNE_HOME)/include
LIBWAYNE_LINK=-L $(LIBWAYNE_HOME) -lwayne$(LIB_OPT) -lm -lpthread $(STACKSIZE) $(STATIC_LINK)
LIBWAYNE_BOTH=$(LIBWAYNE_COMP) $(LIBWAYNE_LINK)

# Name of BLANT source directory
SRCDIR = src

# Put all c files in SRCDIR below. blant-predict.o has its own rule, further down.
BLANT_SRCS = blant.c \
			 blant-window.c \
			 blant-output.c \
			 blant-utils.c \
			 blant-sampling.c \
			 blant-synth-graph.c \
			 importance.c \
			 odv.c \
			 blant-pthreads.c

OBJDIR = _objs
BLANT_CANON_DIR = canon_maps
UND_CANON_DIR := $(BLANT_CANON_DIR)
DIR_CANON_DIR := $(BLANT_CANON_DIR)/directed

# directed canonical values; start with the full range 3..6
# and then remove entries when NO8 or NO7 are set.  This keeps the
# format identical to the `K` list used for undirected data.
K_DIR := 3 4 5 6
ifdef NO8
	# don't build k=6 if undirected 8 is disabled
	K_DIR := $(filter-out 6,$(K_DIR))
endif
ifdef NO_D6
	K_DIR := $(filter-out 6,$(K_DIR))
endif

ifdef NO7
	# NO7 implies NO8; drop both 5 and 6
	K_DIR := $(filter-out 5 6,$(K_DIR))
endif
ifdef NO_DIRECTED
	# skip all directed canon generation
	K_DIR :=
endif
ifneq ("$(wildcard $(SRCDIR)/EdgePredict/blant-predict.c)","")
    BLANT_PREDICT_SRC = $(SRCDIR)/EdgePredict/blant-predict.c
else
    BLANT_PREDICT_SRC = $(SRCDIR)/blant-predict.c
endif


### Generated File Lists ###

# these variables serve only to help in the creation of the generated file lists variables
K := 3 4 5 6 $(SEVEN) $(EIGHT)
ifdef ONLY_DIRECTED
    K :=
endif
alpha_sampling_methods := NBE EBE MCMC
alpha_txts_ud := $(foreach method,$(alpha_sampling_methods),$(UND_CANON_DIR)/alpha_list_$(method))
canon_txt_ud := $(UND_CANON_DIR)/canon_map $(UND_CANON_DIR)/canon_list $(UND_CANON_DIR)/canon-ordinal-to-signature $(UND_CANON_DIR)/orbit_map $(alpha_txts_ud)
canon_txt_dir := $(DIR_CANON_DIR)/canon_map $(DIR_CANON_DIR)/canon_list $(DIR_CANON_DIR)/canon-ordinal-to-signature $(DIR_CANON_DIR)/orbit_map
canon_bin_ud := $(UND_CANON_DIR)/canon_map $(UND_CANON_DIR)/perm_map
canon_bin_dir := $(DIR_CANON_DIR)/canon_map $(DIR_CANON_DIR)/perm_map

# actual generated file lists variables
canon_all_ud := $(foreach k, $(K), $(addsuffix $(k).txt, $(canon_txt_ud)) $(addsuffix $(k).bin, $(canon_bin_ud)))
canon_all_dir := $(foreach k, $(K_DIR), $(addsuffix $(k).txt, $(canon_txt_dir)) $(addsuffix $(k).bin, $(canon_bin_dir)))
canon_all := $(canon_all_ud) $(canon_all_dir)
# subcanon maps go from each k down to k-1, so there is none for k=3.
prev_4 := 3
prev_5 := 4
prev_6 := 5
prev_7 := 6
prev_8 := 7
subcanon_txts := $(foreach k,$(filter-out 3,$(K)),$(UND_CANON_DIR)/subcanon_map$(k)-$(prev_$(k)).txt)
subcanon_txts_dir := $(foreach k,$(filter-out 3,$(K_DIR)),$(DIR_CANON_DIR)/subcanon_map$(k)-$(prev_$(k)).txt)
magic_table_txts := $(foreach k,$(K), orca_jesse_blant_table/UpperToLower$(k).txt)
# The tracked files that make-orca-jesse-blant-table reads.
magic_table_inputs := $(filter-out orca_jesse_blant_table/UpperToLower%,$(wildcard orca_jesse_blant_table/*.txt))

# ehd takes up too much space and isn't used anywhere yet
#ehd_txts := $(foreach k,$(K), $(BLANT_CANON_DIR)/EdgeHammingDistance$(k).txt)

ifeq ($(DYNAMIC_CANON_MAP),1)
    base: ./.notpristine show-gcc-ver libwayne blant
else ifdef ONLY_DIRECTED
    base: ./.notpristine show-gcc-ver libwayne blant $(canon_all)
else
    base: ./.notpristine show-gcc-ver libwayne blant $(canon_all) magic_table test_all
endif

show-gcc-ver:
	gcc$(GCC_VER) -v

./.notpristine:
	@echo '************ READ THIS. REALLY. WE MEAN IT. READ IT AT LEAST ONCE **************'
	@echo "If you haven't already, you should read the README at"
	@echo "	https://github.com/waynebhayes/BLANT#readme"
	@echo "BLANT can sample graphlets of up to k=8 nodes. The lookup table for k=8 can take"
	@echo "up to an hour to generate, but is needed if you want BLANT-seed to work, and so"
	@echo "we make it by default; set NO8=1 to turn it off."
	@echo "The best way to start the very first time is to run the following command:"
	@echo "    ./regression-test-all.sh -make"
	@echo "This may take an hour or more but performs a full battery of tests."
	@echo "The fastest way to get started is to skip k=8 graphlets:"
	@echo "    PAUSE=0 NO8=1 make base"
	@echo "which will make everything needed to get started sampling up to k=7 graphlets".
	@echo "To skip cleaning and re-making libwayne, set NO_CLEAN_LIBWAYNE=1"
	@echo "You will only see this message once on a 'pristine' repo. Pausing $(PAUSE) seconds..."
	@echo '****************************************'
	sleep $(PAUSE)
	@touch .notpristine

most: base Draw
ifneq ($(DYNAMIC_CANON_MAP),1)
most: sub$(BLANT_CANON_DIR)
endif

ifeq ($(DYNAMIC_CANON_MAP),1)
    test_all:
else
    test_all: $(BLANT_CANON_DIR)/test_index_mode $(BLANT_CANON_DIR)/check_maps test_fast
endif

all: most test_all

$(BLANT_CANON_DIR): base $(canon_all) sub$(BLANT_CANON_DIR)

# If a recipe fails, delete the file it was making, so a half-made file never looks finished.
.DELETE_ON_ERROR:

.PHONY: all base most test_all test_fast pristine realclean clean clean_$(BLANT_CANON_DIR) $(BLANT_CANON_DIR) \
	sub$(BLANT_CANON_DIR) magic_table tables show-gcc-ver libwayne Draw FORCE
FORCE:

### libwayne ###

# BLANT only ever runs make in libwayne. The libwayne pinned here has no incremental build:
# its 'make all' deletes every .a file, remakes them all and runs its tests. So it is run once,
# and again only when the contents of a libwayne source file change or one of the .a files is missing.
LW = libwayne
LIBWAYNE_LIBS = $(addprefix $(LW)/,libwayne.a libwayne-g.a libwayne-pg.a libwayne-pg-g.a libwayne-nd.a)
# Dot files, so 'make clean' (rm -rf $(OBJDIR)/*) keeps them.
LIBWAYNE_STAMP = $(OBJDIR)/.libwayne-built
LIBWAYNE_SUM = $(OBJDIR)/.libwayne.sum

libwayne: $(LIBWAYNE_STAMP)

# Checksum of libwayne's sources. intSizes.h is left out, because libwayne makes it.
# LC_ALL=C keeps the file order the same everywhere.
libwayne_files = ls $(LW)/Makefile $(LW)/src/*.c $(LW)/include/*.h | grep -v intSizes.h | LC_ALL=C sort
$(LIBWAYNE_SUM): $(LW)/Makefile FORCE
	mkdir -p $(OBJDIR)
	$(libwayne_files) | xargs cat | cksum > $@.new
	$(call install_sum,$@,$(LIBWAYNE_STAMP))

# This calls plain 'make', not $$(MAKE), on purpose: with $$(MAKE), even 'make -n' would really
# run libwayne's 'make all', which starts by deleting its .a files.
$(LIBWAYNE_STAMP): $(LIBWAYNE_SUM) $(if $(filter-out $(wildcard $(LIBWAYNE_LIBS)),$(LIBWAYNE_LIBS)),FORCE)
	cd $(LW) && make all
	touch $@

$(LIBWAYNE_LIBS): $(LIBWAYNE_STAMP) ;

$(LW)/Makefile:
	echo "Hmm, submodule libwayne doesn't seem to exist; getting it now"
	git submodule init libwayne
	git submodule update libwayne
	(cd libwayne && git checkout master && git pull)

### How each file is compiled ###

# Every compiler call is one "unit". A unit is named by the file it makes, relative to
# $(OBJDIR) for an object (eg blant.o) or to the top for a program (eg Draw/graphette2dot).
# For a unit U:
#   unit_src_U    its source files      (default for an object X.o: $(SRCDIR)/X.c)
#   unit_cc_U     its compiler          (default: $(CC))
#   unit_flags_U  its extra flags       (default: none)
unit_src = $(or $(unit_src_$(1)),$(SRCDIR)/$(patsubst %.o,%.c,$(1)))
unit_cmd = $(or $(unit_cc_$(1)),$(CC)) $(unit_flags_$(1)) $(LIBWAYNE_COMP)
# $(call unit_files,U): every file U is built from, one per line: its sources and every
# header they read, as the compiler finds them (-MM only lists them; it compiles nothing).
unit_files = $(call unit_cmd,$(1)) -MM $(call unit_src,$(1)) | sed -e 's/^[^:]*://' -e 's/\\$$//' | tr ' ' '\n' | grep .

# Each unit U has a checksum file, $(OBJDIR)/U.sum, of its compile command and the contents of
# every file it is built from. It is rewritten only when that checksum changes, and U depends
# on it instead of on its source files. So U is remade when its code really changes, and a
# 'touch' or a 'git checkout' that leaves the contents alone remakes nothing.
#
# $(call install_sum,SUM,OUT): move SUM.new over SUM, but only if it is different.
# make 3.81 (the macOS one) only sees whole seconds. So if SUM changes in the same second as SUM
# or OUT last did, make would not see the change, and OUT would wrongly be kept. In that rare
# case this first waits for the next second (never more than 3), so the new SUM is surely newer.
install_sum = if cmp -s $(1).new $(1); then rm -f $(1).new; else \
	perl -e '$$m = 0; for (@ARGV) { my $$t = (stat)[9]; $$m = $$t if defined $$t && $$t > $$m } for (1..3) { last if time > $$m; sleep 1 }' $(1) $(2); \
	mv -f $(1).new $(1); touch $(1); fi
unit_out = $(if $(filter %.o,$(1)),$(OBJDIR)/$(1),$(1))

# blant's objects. blant-predict.o has its own recipe, further down.
BLANT_OBJ_UNITS := $(BLANT_SRCS:.c=.o)
$(foreach u,$(BLANT_OBJ_UNITS),$(eval unit_flags_$(u) := -D"DYNAMIC_CANON_MAP=$(DYNAMIC_CANON_MAP)"))
# If src/EdgePredict/ exists, this checksum covers only its blant-predict.c and the headers that
# reads; its own sub-make looks after the rest of that folder.
unit_src_blant-predict.o := $(BLANT_PREDICT_SRC)
unit_flags_blant-predict.o := $(PRED_REG_OPT)
# Built here, not with libwayne's own C++ Makefile, which goes by dates.
unit_src_mt19937.o := $(LW)/C++/mt19937.cpp
unit_cc_mt19937.o = $(CXX)
unit_flags_mt19937.o := -std=c++11

# Programs made of one source file plus libblant.o, each with the extra flags it has always had.
# Most of them make the tables in canon_maps/.
PROG_FLAGS_fast-canon-map := -std=gnu11 -O3
PROG_FLAGS_slow-canon-maps :=
PROG_FLAGS_make-orbit-maps :=
PROG_FLAGS_compute-alphas-NBE := -Wall -O3
PROG_FLAGS_compute-alphas-EBE := -Wall -O3
PROG_FLAGS_compute-alphas-MCMC := -Wall -O3
# Currently unused, was the old method for calculating MCMC
PROG_FLAGS_compute-alphas-MCMC-slow := -Wall -O3
PROG_FLAGS_make-subcanon-maps := -Wall
PROGS := fast-canon-map slow-canon-maps make-orbit-maps compute-alphas-NBE compute-alphas-EBE compute-alphas-MCMC compute-alphas-MCMC-slow make-subcanon-maps
$(foreach p,$(PROGS),$(eval unit_flags_$(p).o := $(PROG_FLAGS_$(p))))

# create-bin-data is compiled once for each k, with k built in.
bin_flags = -std=gnu11 -Dkk=$(1) "-DkString=\"$(1)\""
CREATE_BIN_DATA := $(addprefix create-bin-data,$(sort $(K) $(K_DIR)))
$(foreach p,$(CREATE_BIN_DATA),$(eval unit_src_$(p).o := $(SRCDIR)/create-bin-data.c))
$(foreach p,$(CREATE_BIN_DATA),$(eval unit_flags_$(p).o := $(call bin_flags,$(patsubst create-bin-data%,%,$(p)))))

PROG_FLAGS_make-orca-jesse-blant-table := -Wall -std=gnu++11
unit_src_magictable.o := $(SRCDIR)/magictable.cpp
unit_cc_magictable.o = $(CXX)
unit_flags_magictable.o := $(PROG_FLAGS_make-orca-jesse-blant-table)

# Programs compiled straight from their sources in one call.
unit_src_l-k-func-sort-test := $(SRCDIR)/l-k-func-sort-test.c $(SRCDIR)/blant-utils.c
unit_flags_l-k-func-sort-test := -std=gnu11 -DDEBUG_ATTEMPTS
unit_src_l-k-perm-test := $(SRCDIR)/l-k-perm-test.c $(SRCDIR)/blant-utils.c
unit_flags_l-k-perm-test := -std=gnu11 -DDEBUG_ATTEMPTS
unit_src_Draw/graphette2dot := Draw/DrawGraphette.cpp Draw/graphette2dotutils.cpp Draw/Graphette.cpp
unit_cc_Draw/graphette2dot = $(CXX)
unit_flags_Draw/graphette2dot := -std=gnu++11
unit_src_cluster-similarity-graph := $(SRCDIR)/cluster-similarity-graph.c
unit_flags_cluster-similarity-graph := $(SPEED) -Wall

OBJ_UNITS := $(BLANT_OBJ_UNITS) blant-predict.o libblant.o mt19937.o $(PROGS:=.o) $(CREATE_BIN_DATA:=.o) \
	magictable.o makeEHD.o blant-sanity.o synthetic.o syntheticDS.o
PROG_UNITS := l-k-func-sort-test l-k-perm-test Draw/graphette2dot cluster-similarity-graph
UNIT_SUMS := $(addprefix $(OBJDIR)/,$(addsuffix .sum,$(OBJ_UNITS) $(PROG_UNITS)))
# Never delete a checksum file as a throwaway in-between file.
.SECONDARY: $(UNIT_SUMS) $(LIBWAYNE_SUM)

# libwayne makes intSizes.h, which -MM has to find, so it comes first.
$(UNIT_SUMS): $(OBJDIR)/%.sum: FORCE | $(LIBWAYNE_STAMP)
	mkdir -p $(dir $@)
	{ echo '$(subst $(LIBWAYNE_HOME),LIBWAYNE_HOME,$(call unit_cmd,$*) $(call unit_src,$*))'; \
	  $(call unit_files,$*) | LC_ALL=C sort -u | xargs cat; } | cksum > $@.new
	$(call install_sum,$@,$(call unit_out,$*))

### Object Files ###

OBJS = $(addprefix $(OBJDIR)/,$(BLANT_OBJ_UNITS) blant-predict.o)

$(addprefix $(OBJDIR)/,$(filter-out blant-predict.o,$(OBJ_UNITS))): $(OBJDIR)/%: $(OBJDIR)/%.sum
	$(call unit_cmd,$*) -c -o $@ $(call unit_src,$*)

$(OBJDIR)/blant-predict.o: $(OBJDIR)/blant-predict.o.sum
	if [ -f $(SRCDIR)/EdgePredict/Makefile ]; then (CC="$(CC) $(PRED_REG_OPT) $(LIBWAYNE_COMP)"; export CC; OBJDIR="$(OBJDIR)"; export OBJDIR; cd $(SRCDIR)/EdgePredict && $(MAKE)); else $(call unit_cmd,blant-predict.o) -c -o $@ $(call unit_src,blant-predict.o); fi

### Executables ###

# A program is linked again only when one of its objects or libwayne is really remade.
blant: $(OBJS) $(OBJDIR)/libblant.o $(OBJDIR)/mt19937.o $(LIBWAYNE_STAMP)
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(OBJS) $(OBJDIR)/mt19937.o $(LIBWAYNE_LINK)

$(PROGS): %: $(OBJDIR)/%.o $(OBJDIR)/libblant.o $(LIBWAYNE_STAMP)
	$(CC) $(PROG_FLAGS_$@) -o $@ $(OBJDIR)/libblant.o $(OBJDIR)/$@.o $(LIBWAYNE_LINK)

$(CREATE_BIN_DATA): create-bin-data%: $(OBJDIR)/create-bin-data%.o $(OBJDIR)/libblant.o $(LIBWAYNE_STAMP)
	$(CC) $(call bin_flags,$*) -o $@ $(OBJDIR)/libblant.o $< $(LIBWAYNE_LINK)

make-orca-jesse-blant-table: $(OBJDIR)/magictable.o $(OBJDIR)/libblant.o $(LIBWAYNE_STAMP)
	$(CXX) $(PROG_FLAGS_make-orca-jesse-blant-table) -o $@ $(OBJDIR)/magictable.o $(OBJDIR)/libblant.o $(LIBWAYNE_LINK)

# Programs compiled in one call, plus libblant.o:
#   l-k-func-sort-test  debug harness: runs L_K_Func_Sort over the canonicals read from a canon
#                       list file (default canon_maps/canon_list{<k>}.txt) and does nothing else.
#   l-k-perm-test       debug harness: for each canonical given on the command line, runs
#                       L_K_Func_Sort on every decimal in its isomorphism class (all k! relabelings)
#                       and reports whether the per-call attempt count varies across the class.
#   Draw/graphette2dot  draws a graphlet for graphviz.
# The two debug harnesses are NOT part of the normal build (base/most/all); build them only on
# explicit request, eg 'make l-k-func-sort-test'.
$(filter-out cluster-similarity-graph,$(PROG_UNITS)): %: $(OBJDIR)/%.sum $(OBJDIR)/libblant.o $(LIBWAYNE_STAMP)
	$(call unit_cmd,$@) -o $@ $(call unit_src,$@) $(OBJDIR)/libblant.o $(LIBWAYNE_LINK)

# Uses only libwayne's headers, so it links nothing in.
cluster-similarity-graph: $(OBJDIR)/cluster-similarity-graph.sum $(LIBWAYNE_STAMP)
	$(call unit_cmd,$@) -o $@ $(call unit_src,$@)

synthetic: $(OBJDIR)/synthetic.o $(OBJDIR)/syntheticDS.o $(OBJDIR)/libblant.o $(LIBWAYNE_STAMP)
	$(CXX) -o $@ $(OBJDIR)/syntheticDS.o $(OBJDIR)/libblant.o $(OBJDIR)/synthetic.o $(LIBWAYNE_LINK)

makeEHD: $(OBJDIR)/makeEHD.o $(OBJDIR)/libblant.o $(LIBWAYNE_STAMP)
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(OBJDIR)/makeEHD.o $(LIBWAYNE_LINK)

Draw: Draw/graphette2dot

blant-sanity: $(OBJDIR)/blant-sanity.o $(LIBWAYNE_STAMP)
	$(CC) -o $@ $(OBJDIR)/blant-sanity.o $(LIBWAYNE_LINK)

### Generated File Recipes ###

# How the tables stay correct without being remade for nothing:
# 1. Each program that makes tables has a checksum file, $(STAMPS)/<program>.sum, of everything
#    it is built from: its source files and the headers they read, libwayne's source files,
#    and the compile flags. It is rewritten only when that checksum changes, so a table that
#    depends on it is remade only when that code really changed. Like the unit checksums (see
#    "How each file is compiled"), it never goes by the program's date, so a 'make clean' that
#    relinks the program remakes no table.
# 2. Each step, e.g. "make the k=5 canon maps", writes its tables to *.new files and moves each
#    one over the old table only if it is different. So an unchanged table keeps its old date,
#    and the steps after it are not remade. The step then touches its own stamp,
#    $(STAMPS)/<step>.done, which is what the tables themselves depend on.
# 3. If any table of a step is missing, that step runs again.
STAMPS = $(UND_CANON_DIR)/.stamps
DSTAMPS = $(DIR_CANON_DIR)/.stamps

TABLE_PROGS = fast-canon-map make-orbit-maps compute-alphas-NBE compute-alphas-EBE compute-alphas-MCMC make-subcanon-maps make-orca-jesse-blant-table $(CREATE_BIN_DATA)
prog_flags = $(if $(filter create-bin-data%,$(1)),$(call bin_flags,$(patsubst create-bin-data%,%,$(1))),$(PROG_FLAGS_$(1)))
prog_obj = $(if $(filter make-orca-jesse-blant-table,$(1)),magictable,$(1))

# The libwayne path is replaced by a fixed word, so moving the repo does not remake the tables.
# intSizes.h is left out, because libwayne makes it. LC_ALL=C keeps the file order the same everywhere.
$(TABLE_PROGS:%=$(STAMPS)/%.sum): $(STAMPS)/%.sum: % FORCE
	mkdir -p $(STAMPS)
	{ echo '$(subst $(LIBWAYNE_HOME),LIBWAYNE_HOME,$(CC) | $(CXX) | $(LIBWAYNE_BOTH) | $(call prog_flags,$*))'; \
	  { $(call unit_files,$(call prog_obj,$*).o); $(call unit_files,libblant.o); } | LC_ALL=C sort -u | xargs cat; \
	  $(libwayne_files) | xargs cat; } | cksum > $@.new
	$(call install_sum,$@,$(STAMPS)/*.done $(DSTAMPS)/*.done)

# $(call install_new,FILES): move each FILE.new over FILE, but only if it is different.
install_new = for f in $(1); do if cmp -s $$f.new $$f; then rm -f $$f.new; else mv -f $$f.new $$f; fi; done
# For programs that write their tables in place: save each old table first, then put it back
# (with its old date) if the new one is the same, or if the program failed.
save_old = for f in $(1); do if [ -f $$f ]; then cp -p $$f $$f.old; fi; done
keep_old_if_same = for f in $(1); do if [ -f $$f.old ] && cmp -s $$f $$f.old; then mv -f $$f.old $$f; else rm -f $$f.old; fi; done
put_back_old = for f in $(1); do if [ -f $$f.old ]; then mv -f $$f.old $$f; fi; done

# canon_map, canon_list, and canon-...-signature are all made together, because they all come from fast-canon-map.
# If there are several lines for the same canonical, canon_list keeps the one with the most fields.
canon_list_awk = '{if(!($$1 in seen)){seen[$$1]=$$0}else if(NF>split(seen[$$1],tmp,"\t")){seen[$$1]=$$0}}END{n=asorti(seen,sorted,"@ind_num_asc");print n;for(i=1;i<=n;i++)print seen[sorted[i]]}'
# $(call make_canon,DIR,k,ARGS)
define make_canon
	mkdir -p $(1)
	./fast-canon-map $(2) $(3) > $(1)/canon_map$(2).txt.new
	gawk -F '\t' $(canon_list_awk) $(1)/canon_map$(2).txt.new > $(1)/canon_list$(2).txt.sorted
	cut -f1,3- $(1)/canon_list$(2).txt.sorted > $(1)/canon_list$(2).txt.new
	rm -f $(1)/canon_list$(2).txt.sorted
	awk 'NR>1{print NR-2, $$1}' $(1)/canon_list$(2).txt.new > $(1)/canon-ordinal-to-signature$(2).txt.new
	$(call install_new,$(1)/canon_map$(2).txt $(1)/canon_list$(2).txt $(1)/canon-ordinal-to-signature$(2).txt)
endef

$(STAMPS)/canon%.done: $(STAMPS)/fast-canon-map.sum | fast-canon-map
	$(call make_canon,$(UND_CANON_DIR),$*,)
	touch $@
$(UND_CANON_DIR)/canon_map%.txt $(UND_CANON_DIR)/canon_list%.txt $(UND_CANON_DIR)/canon-ordinal-to-signature%.txt: $(STAMPS)/canon%.done ;

# directed counterpart (only k values listed in K_DIR)
$(DSTAMPS)/canon%.done: $(STAMPS)/fast-canon-map.sum | fast-canon-map
	$(call make_canon,$(DIR_CANON_DIR),$*,directed)
	mkdir -p $(DSTAMPS)
	touch $@
$(DIR_CANON_DIR)/canon_map%.txt $(DIR_CANON_DIR)/canon_list%.txt $(DIR_CANON_DIR)/canon-ordinal-to-signature%.txt: $(DSTAMPS)/canon%.done ;

# create-bin-data<k> writes canon_map<k>.bin and perm_map<k>.bin in place.
# $(call make_bin,DIR,k,ARGS)
define make_bin
	$(call save_old,$(1)/canon_map$(2).bin $(1)/perm_map$(2).bin)
	./create-bin-data$(2) $(3) || { $(call put_back_old,$(1)/canon_map$(2).bin $(1)/perm_map$(2).bin); exit 1; }
	$(call keep_old_if_same,$(1)/canon_map$(2).bin $(1)/perm_map$(2).bin)
endef

$(STAMPS)/bin%.done: $(STAMPS)/create-bin-data%.sum $(UND_CANON_DIR)/canon_list%.txt $(UND_CANON_DIR)/canon_map%.txt | create-bin-data%
	$(call make_bin,$(UND_CANON_DIR),$*,)
	touch $@
$(UND_CANON_DIR)/canon_map%.bin $(UND_CANON_DIR)/perm_map%.bin: $(STAMPS)/bin%.done ;

$(DSTAMPS)/bin%.done: $(STAMPS)/create-bin-data%.sum $(DIR_CANON_DIR)/canon_list%.txt $(DIR_CANON_DIR)/canon_map%.txt | create-bin-data%
	$(call make_bin,$(DIR_CANON_DIR),$*,directed)
	touch $@
$(DIR_CANON_DIR)/canon_map%.bin $(DIR_CANON_DIR)/perm_map%.bin: $(DSTAMPS)/bin%.done ;

$(STAMPS)/orbit%.done: $(STAMPS)/make-orbit-maps.sum $(UND_CANON_DIR)/canon_map%.txt $(UND_CANON_DIR)/canon_list%.txt $(UND_CANON_DIR)/canon_map%.bin | make-orbit-maps
	./make-orbit-maps $* > $(UND_CANON_DIR)/orbit_map$*.txt.new
	$(call install_new,$(UND_CANON_DIR)/orbit_map$*.txt)
	touch $@
$(UND_CANON_DIR)/orbit_map%.txt: $(STAMPS)/orbit%.done ;

# directed orbit maps
$(DSTAMPS)/orbit%.done: $(STAMPS)/make-orbit-maps.sum $(DIR_CANON_DIR)/canon_map%.txt $(DIR_CANON_DIR)/canon_list%.txt $(DIR_CANON_DIR)/canon_map%.bin | make-orbit-maps
	./make-orbit-maps $* directed > $(DIR_CANON_DIR)/orbit_map$*.txt.new
	$(call install_new,$(DIR_CANON_DIR)/orbit_map$*.txt)
	touch $@
$(DIR_CANON_DIR)/orbit_map%.txt: $(DSTAMPS)/orbit%.done ;

# canon-upper.sh adds lower-case names (eg alpha_list_nbe5.txt) for the alpha lists.
# $(call make_alpha,METHOD,k)
define make_alpha
	./compute-alphas-$(1) $(2) > $(UND_CANON_DIR)/alpha_list_$(1)$(2).txt.new
	$(call install_new,$(UND_CANON_DIR)/alpha_list_$(1)$(2).txt)
	./canon-upper.sh
endef

$(STAMPS)/alpha-NBE%.done: $(STAMPS)/compute-alphas-NBE.sum $(UND_CANON_DIR)/canon_list%.txt | compute-alphas-NBE
	$(call make_alpha,NBE,$*)
	touch $@
$(STAMPS)/alpha-EBE%.done: $(STAMPS)/compute-alphas-EBE.sum $(UND_CANON_DIR)/canon_list%.txt | compute-alphas-EBE
	$(call make_alpha,EBE,$*)
	touch $@
$(STAMPS)/alpha-MCMC%.done: $(STAMPS)/compute-alphas-MCMC.sum $(UND_CANON_DIR)/canon_list%.txt | compute-alphas-MCMC
	$(call make_alpha,MCMC,$*)
	touch $@
$(UND_CANON_DIR)/alpha_list_NBE%.txt: $(STAMPS)/alpha-NBE%.done ;
$(UND_CANON_DIR)/alpha_list_EBE%.txt: $(STAMPS)/alpha-EBE%.done ;
$(UND_CANON_DIR)/alpha_list_MCMC%.txt: $(STAMPS)/alpha-MCMC%.done ;

sub$(BLANT_CANON_DIR): $(subcanon_txts) $(subcanon_txts_dir)

# undirected and directed subcanon maps, all in one step
$(STAMPS)/subcanon.done: $(STAMPS)/make-subcanon-maps.sum $(canon_all) | make-subcanon-maps
	for k in $(filter-out 3,$(K)); do ./make-subcanon-maps $$k > $(UND_CANON_DIR)/subcanon_map$$k-$$(($$k-1)).txt.new || exit 1; done
	for k in $(filter-out 3,$(K_DIR)); do ./make-subcanon-maps $$k directed > $(DIR_CANON_DIR)/subcanon_map$$k-$$(($$k-1)).txt.new || exit 1; done
	$(call install_new,$(subcanon_txts) $(subcanon_txts_dir))
	touch $@
$(subcanon_txts) $(subcanon_txts_dir): $(STAMPS)/subcanon.done ;

# make-orca-jesse-blant-table writes every UpperToLower<k>.txt in place, in one run.
magic_table: $(magic_table_txts)
$(STAMPS)/magic.done: $(STAMPS)/make-orca-jesse-blant-table.sum $(canon_all) $(magic_table_inputs) | make-orca-jesse-blant-table
	$(call save_old,$(magic_table_txts))
	./make-orca-jesse-blant-table $(if $(EIGHT),8,$(if $(SEVEN),7,6)) || { $(call put_back_old,$(magic_table_txts)); exit 1; }
	ls $(magic_table_txts) > /dev/null
	$(call keep_old_if_same,$(magic_table_txts))
	touch $@
$(magic_table_txts): $(STAMPS)/magic.done ;

# Every step's stamp, by name, so make never treats one as a throwaway in-between file.
TABLE_STAMPS := \
	$(foreach k,$(K),$(addprefix $(STAMPS)/,canon$(k).done bin$(k).done orbit$(k).done alpha-NBE$(k).done alpha-EBE$(k).done alpha-MCMC$(k).done)) \
	$(foreach k,$(K_DIR),$(addprefix $(DSTAMPS)/,canon$(k).done bin$(k).done orbit$(k).done)) \
	$(STAMPS)/subcanon.done $(if $(K),$(STAMPS)/magic.done)
tables: $(TABLE_STAMPS)

# A step must run again if any of its tables is missing, even when its stamp looks up to date.
# $(call rerun_if_missing,STAMP,TABLES)
rerun_if_missing = $(if $(filter-out $(wildcard $(2)),$(2)),$(eval $(1): FORCE))
$(foreach k,$(K), \
	$(call rerun_if_missing,$(STAMPS)/canon$(k).done,$(addsuffix $(k).txt,$(addprefix $(UND_CANON_DIR)/,canon_map canon_list canon-ordinal-to-signature))) \
	$(call rerun_if_missing,$(STAMPS)/bin$(k).done,$(addsuffix $(k).bin,$(canon_bin_ud))) \
	$(call rerun_if_missing,$(STAMPS)/orbit$(k).done,$(UND_CANON_DIR)/orbit_map$(k).txt) \
	$(foreach m,$(alpha_sampling_methods),$(call rerun_if_missing,$(STAMPS)/alpha-$(m)$(k).done,$(UND_CANON_DIR)/alpha_list_$(m)$(k).txt)))
$(foreach k,$(K_DIR), \
	$(call rerun_if_missing,$(DSTAMPS)/canon$(k).done,$(addsuffix $(k).txt,$(addprefix $(DIR_CANON_DIR)/,canon_map canon_list canon-ordinal-to-signature))) \
	$(call rerun_if_missing,$(DSTAMPS)/bin$(k).done,$(addsuffix $(k).bin,$(canon_bin_dir))) \
	$(call rerun_if_missing,$(DSTAMPS)/orbit$(k).done,$(DIR_CANON_DIR)/orbit_map$(k).txt))
$(call rerun_if_missing,$(STAMPS)/subcanon.done,$(subcanon_txts) $(subcanon_txts_dir))
$(call rerun_if_missing,$(STAMPS)/magic.done,$(magic_table_txts))

# Currently unused target
$(BLANT_CANON_DIR)/EdgeHammingDistance%.txt: makeEHD | $(BLANT_CANON_DIR)/canon_list%.txt $(BLANT_CANON_DIR)/canon_map%.bin
	if [ ! -f $(BLANT_CANON_DIR).correct/EdgeHammingDistance$*.txt.xz ]; then ./makeEHD $* > $@; cmp $(BLANT_CANON_DIR).correct/EdgeHammingDistance$*.txt $@; else echo "EdgeHammingDistance8.txt takes weeks to generate; uncompressing instead"; unxz < $(BLANT_CANON_DIR).correct/EdgeHammingDistance$*.txt.xz > $@ && touch $@; fi
	#(cd $(BLANT_CANON_DIR).correct && ls EdgeHammingDistance$*.txt*) | awk '{printf "cmp $(BLANT_CANON_DIR).correct/%s $(BLANT_CANON_DIR)/%s\n",$$1,$$1}' | sh

### Testing ###

# Each test below fails if any one of its checks fails, not just the last one.
test_stamp: blant blant-sanity $(canon_all) $(subcanon_txts)
	@echo Touching test_stamp so $(BLANT_CANON_DIR)/check_maps and $(BLANT_CANON_DIR)/test_index_mode tests only occur if the $(BLANT_CANON_DIR) are changed.
	touch test_stamp

test_fast: blant blant-sanity
	# Run blant sanity test only for MCMC tables. If this fails, chances are EVERYTHING ELSE is wrong. These tests will run every time base is made
	for k in $(K); do if [ -f $(BLANT_CANON_DIR)/canon_map$$k.bin ]; then echo FAST basic sanity for ONLY MCMC with k=$$k; ./blant -q -s MCMC -mi -n 100000 -k $$k networks/syeast.el | sort -n | ./blant-sanity $$k 100000 networks/syeast.el || exit 1; fi; done

# The stamp is touched only after every check passes, so a failed test runs again next time.
$(BLANT_CANON_DIR)/test_index_mode: test_stamp
	# First run blant-sanity for various values of k
	for S in NBE MCMC EBE; do for k in $(K); do if [ -f $(BLANT_CANON_DIR)/canon_map$$k.bin ]; then echo basic sanity check sampling method $$S indexing k=$$k; ./blant -q -s $$S -mi -n 100000 -k $$k networks/syeast.el | sort -n | ./blant-sanity $$k 100000 networks/syeast.el || exit 1; fi; done; done
	touch $@

# Compares every table that has a copy in $(BLANT_CANON_DIR).correct/. A difference only warns,
# because some of those copies are known to be out of date; the check then runs again next time.
$(BLANT_CANON_DIR)/check_maps: test_stamp
	differ=; for f in `ls $(BLANT_CANON_DIR).correct/ | egrep -v 'canon_list2|$(if $(SEVEN),,7|)$(if $(EIGHT),,8|)README|\.[gx]z|EdgeHamming'`; do \
	    [ -d $(BLANT_CANON_DIR).correct/$$f ] && continue; \
	    cmp $(BLANT_CANON_DIR).correct/$$f $(BLANT_CANON_DIR)/$$f || differ="$$differ $$f"; \
	done; \
	if [ -n "$$differ" ]; then echo "WARNING: check_maps: these differ from $(BLANT_CANON_DIR).correct/:$$differ" >&2; else touch $@; fi

### Cleaning ###

# On a disk that ignores case (macOS), 'synthetic' is the Synthetic/ folder; leave that alone.
clean:
	/bin/rm -f *.[oa] blant create-bin-data3 create-bin-data4 create-bin-data5 create-bin-data6 create-bin-data7 create-bin-data8 canon-sift fast-canon-map slow-canon-maps make-orbit-maps compute-alphas-MCMC-slow compute-alphas-MCMC compute-alphas-NBE compute-alphas-EBE make-orca-jesse-blant-table Draw/graphette2dot blant-sanity make-subcanon-maps makeEHD cluster-similarity-graph test_stamp l-k-func-sort-test l-k-perm-test $(BLANT_CANON_DIR)/check_maps $(BLANT_CANON_DIR)/test_index_mode
	[ ! -f synthetic ] || /bin/rm -f synthetic
	/bin/rm -rf $(OBJDIR)/*
	/bin/rm -rf $(DIR_CANON_DIR)/* || true

realclean:
	echo "'realclean' is now called 'pristine'; try again"
	false

pristine: clean clean_$(BLANT_CANON_DIR)
ifndef NO_CLEAN_LIBWAYNE
	cd $(LIBWAYNE_HOME); $(MAKE) clean
endif
	[ ! -d $(BLANT_CANON_DIR) ] || find $(BLANT_CANON_DIR) -mindepth 1 -not -type d -delete; /bin/rm -f .notpristine .firsttime # .firsttime is the old name but remove it anyway
	@echo "Finding all python crap and removing it... this may take awhile..." >/dev/null
	./scripts/delete-python-shit.sh $(UNAME)

clean_$(BLANT_CANON_DIR):
	/bin/rm -f $(BLANT_CANON_DIR)/*[3-7].* # don't remove 8 since it takes too long to create
	/bin/rm -f $(DIR_CANON_DIR)/*[3-6].* || true # directed only goes through k=6
	/bin/rm -f orca_jesse_blant_table/UpperToLower*.txt
