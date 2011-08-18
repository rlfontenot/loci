
include $(LOCI_BASE)/Loci.conf

INCLUDES=-I.

COPT = $(CC_LIB_FLAGS) $(COPTLESS)

OBJS = \
	b64.o blas.o csr.o error.o evaluate.o fkvkselect.o fs.o getopt.o \
	gkregex.o htable.o io.o itemsets.o mcore.o memory.o omp.o pdb.o \
	pqueue.o random.o seq.o sort.o string.o timers.o tokenizer.o \
	util.o

libgk.a: $(OBJS) 
	$(AR) $(ARFLAGS) libgk.a $(OBJS) 
	$(RANLIB) libgk.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libgk.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libgk.so $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)


distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
