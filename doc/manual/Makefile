#----------------------------------------------------------------------------
# Makefile for the Manual of cpmd.x (plane wave electronic calculation)
#----------------------------------------------------------------------------
#
SHELL=/bin/sh
HTMLDIR=manual
LTXFILE=manual
#----------------------------------------------------------------------------
# Generate the manual
#----------------------------------------------------------------------------
all:         pdf

clean:      
	     rm -f $(LTXFILE).aux $(LTXFILE).ilg $(LTXFILE).log $(LTXFILE).idx $(LTXFILE).ind $(LTXFILE).out $(LTXFILE).toc 
	     @$(MAKE) cppsource IS_HTML=0

realclean:   
	     @$(MAKE) clean
	     rm -f $(LTXFILE).dvi $(LTXFILE).pdf $(LTXFILE).ps
	     rm -rf $(HTMLDIR)

ps:          $(LTXFILE).ps
	     @$(MAKE) clean
	     @rm -f $(LTXFILE).dvi

pdf:         $(LTXFILE).pdf
	     @$(MAKE) clean
	     @rm -f $(LTXFILE).dvi

html:        $(LTXFILE).html
	     @$(MAKE) clean
	     @rm -f $(LTXFILE).dvi

cppsource:
ifeq ($(IS_HTML),1)
	     @awk '{ if (match($$0,"%  HTML TAG -  DO NOT REMOVE  %")!=0) {print "\\usepackage{html}"} else {print} }' manual.tex > manual.tmp 
	     @grep -i "\\usepackage{html}" manual.tmp > /dev/null ; if [ $$? -eq 0 ] ; then mv manual.tmp manual.tex ; else echo "Inconsistency! Bug! Report!"; exit 9; fi
else
	     $(shell awk '{ if (match($$0,"\\usepackage{html}")!=0) {print "%  HTML TAG -  DO NOT REMOVE  %"} else {print} }' manual.tex > manual.tmp )
	     @grep -i "%  HTML TAG -  DO NOT REMOVE  %" manual.tmp > /dev/null ; if [ $$? -eq 0 ] ; then mv manual.tmp manual.tex ; else echo "Inconsistency! Bug! Report!"; exit 9; fi
endif

manual.dvi:  
	     @$(MAKE) cppsource 
	     rm -f $(LTXFILE).dvi
	     latex $(LTXFILE).tex
	     latex $(LTXFILE).tex
	     makeindex $(LTXFILE)
	     latex $(LTXFILE).tex
	     latex $(LTXFILE).tex

manual.pdf:  
	     @$(MAKE) cppsource IS_HTML=0
	     rm -f $(LTXFILE).pdf
	     pdflatex $(LTXFILE).tex
	     pdflatex $(LTXFILE).tex
	     makeindex $(LTXFILE)
	     pdflatex $(LTXFILE).tex
	     pdflatex $(LTXFILE).tex

manual.html: 
	     @$(MAKE) manual.dvi IS_HTML=1
	     rm -rf $(HTMLDIR)
	     mkdir $(HTMLDIR)
	     latex2html -dir $(HTMLDIR) $(LTXFILE).tex
	     latex2html -dir $(HTMLDIR) $(LTXFILE).tex

manual.ps:   
	     @$(MAKE) manual.dvi IS_HTML=0
	     rm -f $(LTXFILE).ps
	     dvips -o $@ $(LTXFILE).dvi

