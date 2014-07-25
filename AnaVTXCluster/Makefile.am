AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include 

lib_LTLIBRARIES = libAnaVTXCluster.la   

libAnaVTXCluster_la_LIBADD = \
  -L$(OFFLINE_MAIN)/lib -L$(libdir) -lphool -lSubsysReco @ROOTLIBS@ 

noinst_HEADERS = AnaVTXClusterLinkDef.h 

include_HEADERS = AnaVTXCluster.h 

libAnaVTXCluster_la_SOURCES = \
  AnaVTXCluster.C \
  AnaVTXCluster_Dict.C 

AnaVTXCluster_Dict.C: AnaVTXCluster.h AnaVTXClusterLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

