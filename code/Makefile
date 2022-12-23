#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for full recourse kidney exchange problem
#@author  Christopher Hojny

#-----------------------------------------------------------------------------
# path
#-----------------------------------------------------------------------------

SCIPDIR         =       $(SCIP_PATH)

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP (need to do this twice, once to
# find the correct binary, then, after getting the correct flags from the
# binary (which is necessary since the ZIMPL flags differ from the default
# if compiled with the SCIP Optsuite instead of SCIP), we need to set the
# compile flags, e.g., for the ZIMPL library, which is again done in make.project
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
SCIPVERSION			:=$(shell $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION) -v | sed -e 's/$$/@/')
override ARCH		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ARCH=\([^@]*\).*/\1/')
override EXPRINT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* EXPRINT=\([^@]*\).*/\1/')
override GAMS		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GAMS=\([^@]*\).*/\1/')
override GMP		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GMP=\([^@]*\).*/\1/')
override SYM		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')
override IPOPT		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPT=\([^@]*\).*/\1/')
override IPOPTOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPTOPT=\([^@]*\).*/\1/')
override LPSCHECK	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSCHECK=\([^@]*\).*/\1/')
override LPSOPT 	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSOPT=\([^@]*\).*/\1/')
override NOBLKBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKBUFMEM=\([^@]*\).*/\1/')
override NOBLKMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKMEM=\([^@]*\).*/\1/')
override NOBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBUFMEM=\([^@]*\).*/\1/')
override PARASCIP	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* PARASCIP=\([^@]*\).*/\1/')
override READLINE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* READLINE=\([^@]*\).*/\1/')
override SANITIZE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SANITIZE=\([^@]*\).*/\1/')
override ZIMPL		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPL=\([^@]*\).*/\1/')
override ZIMPLOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPLOPT=\([^@]*\).*/\1/')
override ZLIB		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZLIB=\([^@]*\).*/\1/')
include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# default settings for test runs

VERSION		=	1.0
CONTINUE	=	false
LOCK		=	false
TEST    	=
SETTINGS    	=   	default
NATTACKS	=	2
TIME		=	3600
NODES		=	2100000000
MEM 		=	24000
ONLYPRE		=	false
DISPFREQ	=	10000
TABUITERS	=	1000000
FEASTOL		=	default
PERMUTE         =       0
SETCUTOFF	=	false
CLIENTTMPDIR    =	/home/mcs001/s145344/mytmp/
PPN		=	-1
QUEUE		=	mcs.default.q

# symmetry settings for CPLEX - can be -1 (automatic), 0 (off), 1-5 (increasingly aggressive)
SYMM		=   	5
CPLEXBIN	= 	/Users/bin/workspace/cplex121_0.darwin/cplex121/bin/x86_darwin8_gcc4.0/cplex
GUROBIBIN	= 	gurobi.sh
THREADS		=	1
NOWAITCLUSTER	=	1

GETSTARGITHASH	=	$(shell git describe --always --dirty  | sed 's/^.*-g//')

USRCXXFLAGS	+=      -DSTARGITHASH=\"$(GETSTARGITHASH)\"


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	fullrecoursekidney
CMAINOBJ	=	auxiliaryStructures.o \
			kidneyParams.o \
			kidneyPlugins.o \
			probdata_master_kidney.o \
			probdata_benders.o \
			probdata_benders_subproblem.o \
			probdata_glorie.o \
			problem_master_kidneyexchange.o \
			problem_benders.o \
			problem_glorie.o \
			solveMasterProblem.o \
			probdata_benders_picef.o \
			probdata_benders_subproblem_picef.o \
			probdata_master_kidney_picef.o \
			problem_master_kidneyexchange_picef.o \
			problem_benders_picef.o \
			solveMasterPICEFProblem.o



CXXMAINOBJ	=	main.o \
			read_graph.o \
			find_graphstructures.o

MAINSRC		=	$(addprefix $(SRCDIR)/,$(CMAINOBJ:.o=.c))
MAINSRC		+=	$(addprefix $(SRCDIR)/,$(CXXMAINOBJ:.o=.cpp))

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(CMAINOBJ))
MAINOBJFILES	+=	$(addprefix $(OBJDIR)/,$(CXXMAINOBJ))

#-----------------------------------------------------------------------------
# External libraries
#-----------------------------------------------------------------------------

# FLAGS		+=	-DBMS_NOBLOCKMEM -DSCIP_NOBUFFERMEM
FLAGS		+=
LDFLAGS		+=

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)


.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -I$(SCIPDIR) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE $$i; \
			done'
.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

.PHONY: doc
doc:
		cd doc; $(DOXY) kidney.dxy

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: test
test:
		cd check; \
	cd check; \
	$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(ONLYPRE) $(NATTACKS);

testcluster:
	cd check; \
	$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) \
	$(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(NATTACKS) \
	$(DISPFREQ) $(CONTINUE) $(VERSION) \
	$(LPS) $(QUEUE) $(PPN) $(CLIENTTMPDIR) $(NOWAITCLUSTER) $(PERMUTE);

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h src/*.cpp $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

-include	$(MAINOBJFILES:.o=.d)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) $(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(DFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) $(DFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
