default : optados

all : optados docs

optados : 
	make -C ./src all

# Utility targets
.PHONY: clean veryclean

install:
	make -C ./src install

docs:
	make -C ./documents all

clean:
	make -C ./src clean 
	make -C ./documents clean


veryclean:
	make  -C ./src veryclean
	make -C ./documents veryclean

dist:
	 @(tar cf - \
                ./src/*.?90 \
                ./examples/*/*.cell \
                ./examples/*/*.param \
                ./examples/*/*.odi \
                ./documents/*.tex \
                ./documents/*.pdf \
		./documents/Makefile \
		./documents/THINGS_TO_DO.txt \
                ./Makefile \
		./make.sys \
		COPYING \ 
	| gzip -c > \
		./optados.tar.gz)

