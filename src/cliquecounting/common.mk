CC       := gcc
INCLUDES := -I $(ESCAPE_HOME)
DEFINES  := 
CFLAGS   := -Wall -std=c++11 -g -O3 #-O3 -Werror
LDFLAGS  := -L $(ESCAPE_HOME)
LDLIBS   := -lescape -lstdc++ -lm


DEPDIR   := .deps
DEPENDS  := $(OBJECTS:%.o=$(DEPDIR)/%.dep)
MAKEDEP  = $(CC) -MM $(INCLUDES) $(DEFINES) $(CFLAGS) 


all : dep $(TARGETS)


dep : $(DEPENDS)


$(DEPDIR):
	mkdir -p $(DEPDIR)


$(DEPDIR)/%.dep : %.cpp | $(DEPDIR)
	$(MAKEDEP) $< -o /dev/stdout | sed 's,:, $(DEPDIR)/$*.dep:,' > $@

COMPILE = $(CC) -c $(DEFINES) $(INCLUDES) $(CFLAGS) -o $@ $<
COMPILE_AND_LINK = $(CC) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS)

%.o : %.cpp Makefile $(ESCAPE_HOME)/common.mk
	$(COMPILE)


clean:
	rm -f $(OBJECTS) $(TARGETS)


cleandep:
	rm -f $(DEPENDS)
	rmdir $(DEPDIR)


debug:
	@echo abc$(filter clean%,$(MAKECMDGOALS))

ifeq (abc$(filter clean%,$(MAKECMDGOALS)),abc)
#ifneq ($(MAKECMDGOALS),cleandep)
-include $(DEPENDS)
endif
