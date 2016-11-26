def: all

# dirs
BDIR = bin
SDIR = src
ODIR = obj

# tools
CC = gcc
CPPC = g++
RM = rm
MKDIR = mkdir
# default parameters
CFLAGS_COMMON = -std=c++0x -Wno-unused-variable
CFLAGS_DBG = $(CFLAGS_COMMON) -ggdb -Wall -O0
CFLAGS_OPT = $(CFLAGS_COMMON) -O4
LFLAGS_DBG = 
LFLAGS_OPT = -s

# targets

# parameters
NWS_BDIR = $(BDIR)
NWS_SDIR = $(SDIR)
NWS_ODIR = $(ODIR)/nws
NWS_IDIRS = 
NWS_OBJ =
NWS_OBJ += $(NWS_ODIR)/main.o
NWS_OBJ += $(NWS_ODIR)/timer.o
NWS_OBJ += $(NWS_ODIR)/nwsimplex_misc.o
NWS_OBJ += $(NWS_ODIR)/dimacs.o
NWS_OBJ += $(NWS_ODIR)/solver.o
NWS_CFLAGS = $(CFLAGS_DBG)
NWS_BIN = $(BDIR)/nwsimplex
NWS_BIN_OBJ = $(NWS_OBJ)
NWS_LDIRS = 
NWS_LIBS = 
NWS_LFLAGS = $(LFLAGS_DBG)
# rules
$(NWS_ODIR): 
	$(MKDIR) -p $(NWS_ODIR)
$(NWS_ODIR)/%.o: $(NWS_SDIR)/%.cpp | $(NWS_ODIR)
	$(CPPC) $(NWS_CFLAGS) -MMD -c $< -o $@ $(NWS_IDIRS)
-include $(NWS_ODIR)/*.d
$(NWS_BIN): $(NWS_BIN_OBJ) | $(BDIR)
	$(CPPC) -o $(NWS_BIN) $(NWS_BIN_OBJ) $(NWS_LFLAGS) $(NWS_LIBS) $(NWS_LDIRS)


# general targets
$(BDIR):
	$(MKDIR) -p $(BDIR)
all: $(NWS_BIN)
clean: cleanobj cleanbin
cleanobj:
	$(RM) $(NWS_OBJ) -f
	$(RM) $(ODIR) -rf
cleanbin:
	$(RM) $(NWS_BIN) -f
