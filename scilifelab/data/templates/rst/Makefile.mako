# Makefile for rst documentation
#

# You can set these variables from the command line.
OPTS    = -s ${stylefile}
RST2PDF = rst2pdf

.PHONY: latexpdf

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo " "
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  file.pdf   to make a given pdf file"

RSTFILES := $(wildcard *.rst)
TRG	= $(RSTFILES:%.rst=%.pdf)

latexpdf:
	@echo "Running LaTeX files through pdflatex..."
	@echo $(TRG)
	$(MAKE) $(TRG)

%%.pdf: %.rst
	$(RST2PDF) $(OPTS) $<
	@echo
	@echo "Build finished."

