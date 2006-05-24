LOCI_BASE=/usr/local/Loci-Linux-x86_64-mpic++-rel-3-0-p8/
include $(LOCI_BASE)/Loci.conf
COPT = -g

# List the name of your compiled program here
TARGET = lpp


# List your object files here
OBJS =  lpp.o main.o


LOCAL_LIBS =

JUNK = *~ core ti_files ii_files rii_files 

default:
	$(MAKE) $(TARGET) 

all: $(TARGET) 


$(TARGET): $(OBJS) 
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
FRC : 

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak *.o $(OTHER_APPS)
