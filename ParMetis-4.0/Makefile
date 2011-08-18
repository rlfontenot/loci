
default:
	(cd GKLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libgk.so) 
	(cd METISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libmetis.so) 
	(cd ParMETISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libparmetis.so) 


static:
	(cd GKLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libgk.a) 
	(cd METISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libmetis.a) 
	(cd ParMETISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libparmetis.a) 


clean:
	(cd GKLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" clean ) 
	(cd METISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" clean ) 
	(cd ParMETISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" clean ) 

distclean:
	(cd GKLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" distclean ) 
	(cd METISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" distclean ) 
	(cd ParMETISLib ; $(MAKE) LOCI_BASE="$(LOCI_BASE)" distclean ) 



