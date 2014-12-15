AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libAnaAlignmentProd.la

libAnaAlignmentProd_la_LIBADD = \
  -L$(OFFLINE_MAIN)/lib -L$(libdir) -lphool -lSubsysReco @ROOTLIBS@

noinst_HEADERS = AnaAlignmentProdLinkDef.h

include_HEADERS = AnaAlignmentProd.h

libAnaAlignmentProd_la_SOURCES = \
  AnaAlignmentProd.C \
  AnaAlignmentProd_Dict.C

AnaAlignmentProd_Dict.C: AnaAlignmentProd.h AnaAlignmentProdLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

