# Set up LOCI_BASE to point to the location where Loci is installed
# on your system
LOCI_BASE = /ccs/mssl/general/loci/debug/Loci-sun4-gcc-rel-1-0-beta-6

# Include Loci.conf to get all of the rules for compiling Loci programs
include $(LOCI_BASE)/Loci.conf


TARGET = entities containers relations

# List your object files here

EX1_OBJS = entities.o
EX2_OBJS = containers.o
EX3_OBJS = relations.o
EX4_OBJS = helpers.o
EX5_OBJS = lnklst.o
EX6_OBJS = life.o

OBJS = $(EX1_OBJS) $(EX2_OBJS) $(EX3_OBJS) $(EX4_OBJS) $(EX5_OBJS) $(EX6_OBJS)

# Junk files created by compilers that may need to be cleaned up
JUNK = *~ *.ii core ti_files 

default:
	$(MAKE) $(TARGET)

# To compile a Loci program, use the .cc or .C prefix on your Loci source 
# modules and link using a command similar to below.
entities: $(EX1_OBJS)
	$(LD) -o entities $(EX1_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

containers: $(EX2_OBJS)
	$(LD) -o containers $(EX2_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

relations: $(EX3_OBJS)
	$(LD) -o relations $(EX3_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

helpers: $(EX4_OBJS)
	$(LD) -o helpers $(EX4_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

lnklst: $(EX5_OBJS)
	$(LD) -o lnklst $(EX5_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

life:	$(EX6_OBJS)
	$(LD) -o life $(EX6_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

esoteric: esoteric.o
	$(LD) -o esoteric esoteric.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak

make.depend:
	$(MAKEDEPEND) $(CFLAGS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

depend: make.depend

include make.depend
