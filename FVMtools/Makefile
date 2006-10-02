include $(LOCI_BASE)/Loci.conf

COPT = $(COPTLESS)

# List the name of your compiled program here
TARGET = extract

OTHER_APPS = cobalt2xdr solidMesh2xdr solidMesh2xdrMPI plot3d2xdr pb xdr2cobalt xdropt 

# List your object files here
OBJS =  fv_extract.o tec_extract.o extract.o extract_2dgv.o ensight_extract.o


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

cobalt2xdr: cobalt2xdr.o
	$(LD) -o cobalt2xdr cobalt2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

plot3d2xdr: plot3d2xdr.o
	$(LD) -o plot3d2xdr plot3d2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

solidMesh2xdr: solidMesh2xdr.o
	$(LD) -o solidMesh2xdr solidMesh2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

solidMesh2xdrMPI: solidMesh2xdrMPI.o
	$(LD) -o solidMesh2xdrMPI solidMesh2xdrMPI.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

xdropt: xdropt.o
	$(LD) -o xdropt xdropt.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

pb: pb.o predicates.o
	$(LD) -o pb pb.o predicates.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
predicates.o:
	$(CC) $(COPT) -c $(DEFINES) predicates.c

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK)

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
