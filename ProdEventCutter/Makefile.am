AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libProdEventCutter.la

libProdEventCutter_la_LIBADD = \
  -L$(OFFLINE_MAIN)/lib -L$(libdir) -lphool -lSubsysReco @ROOTLIBS@

noinst_HEADERS = ProdEventCutterLinkDef.h

include_HEADERS = ProdEventCutter.h

libProdEventCutter_la_SOURCES = \
  ProdEventCutter.C \
  ProdEventCutter_Dict.C

ProdEventCutter_Dict.C: ProdEventCutter.h ProdEventCutterLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

