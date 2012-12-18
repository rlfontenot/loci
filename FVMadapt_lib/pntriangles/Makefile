#include ../chem.conf
#LOCI_BASE=/usr/tmp/qxue/Loci
include $(LOCI_BASE)/Loci.conf

#COPTLESS=-g
COPT = $(COPTLESS)

INCLUDES = -I../src/include

LOCAL_LIBS =

JUNK = *~ core ti_files ii_files rii_files 

GEOMOBJS = compute_geom.o coeff.o geom_output.o  set_bnodes_index.o

OBJS = adjust_pos.o surface_delta.o mesh_deform.o 

default: computegeom adjustpos getgeom viscousFix


FRC:


computegeom: $(GEOMOBJS)
	$(LD) -o computegeom $(GEOMOBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
adjustpos: $(OBJS)
	$(LD) -o adjustpos $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
getgeom: get_geom.o
	$(LD) -o getgeom get_geom.o 
viscousFix: viscousFix.o
	$(LD) -o viscousFix viscousFix.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


clean:
	rm -fr $(GEOMOBJS) $(OBJS) $(JUNK) make.depend make.depend.bak *.o $(OTHER_APPS)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

distclean:
	rm -f $(DEPEND_FILES)
	rm -f computegeom adjustpos getgeom viscousFix
	rm -fr $(GEOMOBJS) $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak *.o $(OTHER_APPS)

DEPEND_FILES=$(subst .o,.d,$(GEOMOBJS) $(OBJS)) $(subst .o,.d,$(EOSGEOMOBJS)$(EOSOBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
