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
NMFS_BDIR = $(BDIR)
NMFS_SDIR = $(SDIR)
NMFS_ODIR = $(ODIR)/nmfs
NMFS_IDIRS = 
NMFS_OBJ =
NMFS_OBJ += $(NMFS_ODIR)/main.o
NMFS_OBJ += $(NMFS_ODIR)/timer.o
NMFS_OBJ += $(NMFS_ODIR)/nmfs_util.o
NMFS_OBJ += $(NMFS_ODIR)/dimacs.o
NMFS_OBJ += $(NMFS_ODIR)/solver.o
NMFS_CFLAGS = $(CFLAGS_OPT)
NMFS_BIN = $(BDIR)/nmfs
NMFS_BIN_OBJ = $(NMFS_OBJ)
NMFS_LDIRS = 
NMFS_LIBS = 
NMFS_LFLAGS = $(LFLAGS_OPT)
# rules
$(NMFS_ODIR): 
	$(MKDIR) -p $(NMFS_ODIR)
$(NMFS_ODIR)/%.o: $(NMFS_SDIR)/%.cpp | $(NMFS_ODIR)
	$(CPPC) $(NMFS_CFLAGS) -MMD -c $< -o $@ $(NMFS_IDIRS)
-include $(NMFS_ODIR)/*.d
$(NMFS_BIN): $(NMFS_BIN_OBJ) | $(BDIR)
	$(CPPC) -o $(NMFS_BIN) $(NMFS_BIN_OBJ) $(NMFS_LFLAGS) $(NMFS_LIBS) $(NMFS_LDIRS)


# general targets
$(BDIR):
	$(MKDIR) -p $(BDIR)
all: $(NMFS_BIN)
clean: cleanobj cleanbin
cleanobj:
	$(RM) $(NMFS_OBJ) -f
	$(RM) $(ODIR) -rf
cleanbin:
	$(RM) $(NMFS_BIN) -f
