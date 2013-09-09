MANUAL_NAME := user_guide 

all:
	latex $(MANUAL_NAME).tex
	latex $(MANUAL_NAME).tex
	bibtex $(MANUAL_NAME)
	latex $(MANUAL_NAME)/tex
	dvipdfm $(MANUAL_NAME)

clean:
	rm -f *.toc *.aux *.dvi *.log

veryclean: 
	rm -f *.pdf
