MANUAL_NAME := user_guide

all:
	latex $(MANUAL_NAME).tex
	latex $(MANUAL_NAME).tex
	bibtex $(MANUAL_NAME)
	latex $(MANUAL_NAME).tex
	dvips $(MANUAL_NAME).dvi
	ps2pdf $(MANUAL_NAME).ps

clean:
	rm -f *.toc *.aux *.dvi *.log

veryclean: 
	rm -f *.pdf
