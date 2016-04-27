
default:
	$(MAKE) -C GKLib LOCI_BASE="$(LOCI_BASE)" libgk.so
	$(MAKE) -C METISLib LOCI_BASE="$(LOCI_BASE)" libmetis.so
	$(MAKE) -C ParMETISLib LOCI_BASE="$(LOCI_BASE)" libparmetis.so 


static:
	$(MAKE) -C GKLib LOCI_BASE="$(LOCI_BASE)" libgk.a
	$(MAKE) -C METISLib LOCI_BASE="$(LOCI_BASE)" libmetis.a 
	$(MAKE) -C ParMETISLib LOCI_BASE="$(LOCI_BASE)" libparmetis.a 


clean:
	$(MAKE) -C GKLib LOCI_BASE="$(LOCI_BASE)" clean 
	$(MAKE) -C METISLib LOCI_BASE="$(LOCI_BASE)" clean  
	$(MAKE) -C ParMETISLib LOCI_BASE="$(LOCI_BASE)" clean  

distclean:
	$(MAKE) -C GKLib LOCI_BASE="$(LOCI_BASE)" distclean 
	$(MAKE) -C METISLib LOCI_BASE="$(LOCI_BASE)" distclean  
	$(MAKE) -C ParMETISLib LOCI_BASE="$(LOCI_BASE)" distclean  

.PHONY: default static clean distclean