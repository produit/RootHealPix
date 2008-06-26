# $Id: Makefile,v 1.2 2008/06/26 00:51:38 oxon Exp $
# Author: Akira Okumura 2008/06/20

###############################################################################
#  Copyright (C) 2008-, Akira Okumura                                         #
#  All rights reserved.                                                       #
###############################################################################

include $(ROOTSYS)/test/Makefile.arch

NAME	:=	RootHealPix
DEPEND	:=	

SRCDIR	:=	src
INCDIR	:=	include

DICT	:=	$(NAME)Dict
DICTS	:=	$(SRCDIR)/$(NAME)Dict.$(SrcSuf)
DICTI	:=	$(SRCDIR)/$(NAME)Dict.h
DICTO	:=	$(SRCDIR)/$(NAME)Dict.$(ObjSuf)

INCS	:=	$(filter-out $(INCDIR)/LinkDef.h,$(wildcard $(INCDIR)/*.h))
SRCS	:=	$(filter-out $(SRCDIR)/$(DICT).%,$(wildcard $(SRCDIR)/*.$(SrcSuf)))
OBJS	:=	$(patsubst %.$(SrcSuf),%.$(ObjSuf),$(SRCS)) $(DICTO)

ifeq ($(PLATFORM),macosx)
LIB	=	lib$(NAME).$(DllSuf)
LIB_SYMBOLIC=	$(subst .$(DllSuf),.so,$(LIB))
endif
ifeq ($(PLATFORM),win32)
LIB	=	lib$(NAME).lib
else
LIB	=	lib$(NAME).$(DllSuf)
endif

RMAP	=	lib$(NAME).rootmap

EXTLIBS	=	-lcfitsio

.SUFFIXES:	.$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:		all clean doc html

all:		$(RMAP)

$(LIB):		$(OBJS)
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		echo $(LDFLAGS)
		echo $(SOFLAGS)
		echo $(EXTLIBS)
		$(LD) $(SOFLAGS) $(EXTLIBS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(LIB_SYMBOLIC)
else
		echo $(LDFLAGS)
		echo $(SOFLAGS)
		echo $(EXTLIBS)
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $(EXTLIBS) $^ \
		$(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(EVENTLIB)
                $(LD) $(EXTLIBS) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) $(OutPutOpt) $@
else
		$(LD) $(EXTLIBS) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
		@echo "$@ done"


$(SRCDIR)/%.$(ObjSuf):	$(SRCDIR)/%.$(SrcSuf) $(INCDIR)/%.h
		@echo "Compiling" $<
		$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(DICTS):	$(INCS) $(INCDIR)/LinkDef.h
		@echo "Generating dictionary ..."
		rootcint -f $@ -c -p $^

$(DICTO):	$(DICTS)
		@echo "Compiling" $<
		$(CXX) $(CXXFLAGS) -I. -I/usr/X11/include -c $< -o $@

$(RMAP):	$(LIB) $(INCDIR)/LinkDef.h
		rlibmap -f -o $@ -l $(LIB) -d $(DEPEND) -c $(INCDIR)/LinkDef.h

doc:	all htmldoc

htmldoc:
	sh mkhtml.sh

clean:
		rm -rf $(LIB) $(LIB_SYMBOLIC) $(OBJS) $(DICTI) $(DICTS) $(DICTO) htmldoc
