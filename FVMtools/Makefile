include $(LOCI_BASE)/Loci.conf

COPT = $(COPTLESS)
#COPT = -g -O0 -DDEBUG -DBOUNDS_CHECK

# List the name of your compiled program here
TARGET = extract

OTHER_APPS = cfd++2xdr cobalt2xdr solidMesh2xdr plot3d2xdr pb xdr2cobalt xdr2vog vog2xdr ugrid2vog

# List your object files here
OBJS =  extract.o extract_ensight.o extract_2dgv.o extract_ascii.o extract_tecplot.o extract_fv.o



LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE) $(TARGET) $(OTHER_APPS)

all: $(TARGET) $(OTHER_APPS)


$(TARGET): $(OBJS) 
	$(LD) -o $(TARGET) $(OBJS) $(TEC_LIBS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
FRC : 

xdr2cobalt: xdr2cobalt.o
	$(LD) -o xdr2cobalt xdr2cobalt.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

xdr2vog: xdr2vog.o vogtools.o
	$(LD) -o xdr2vog xdr2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ugrid2vog: ugrid2vog.o vogtools.o
	$(LD) -o ugrid2vog ugrid2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vog2xdr: vog2xdr.o
	$(LD) -o vog2xdr vog2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cobalt2xdr: cobalt2xdr.o
	$(LD) -o cobalt2xdr cobalt2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cfd++2xdr: cfd++2xdr.o
	$(LD) -o cfd++2xdr cfd++2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

plot3d2xdr: plot3d2xdr.o
	$(LD) -o plot3d2xdr plot3d2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

solidMesh2xdr: solidMesh2xdr.o
	$(LD) -o solidMesh2xdr solidMesh2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

pb: pb.o predicates.o
	$(LD) -o pb pb.o predicates.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

predicates.o:
	$(CC) $(COPT) -c $(DEFINES) predicates.c

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(OTHER_APPS) $(JUNK) *.o

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) $(OTHER_APPS) *.o

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
