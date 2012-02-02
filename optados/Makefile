VER := 1.0.0

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

dist: optados-$(VER).tar.gz

CURRENT_DIR := $(shell basename ${PWD})

optados-$(VER).tar.gz: 
	 cd .. && tar -czf $(CURRENT_DIR)/optados-$(VER).tar.gz \
                $(CURRENT_DIR)/src/*.?90 \
                $(CURRENT_DIR)/src/Makefile \
                $(CURRENT_DIR)/examples/*/*.cell \
                $(CURRENT_DIR)/examples/*/*.param \
                $(CURRENT_DIR)/examples/*/*.odi \
                $(CURRENT_DIR)/documents/*.tex \
                $(CURRENT_DIR)/documents/*.pdf \
		$(CURRENT_DIR)/documents/Makefile \
		$(CURRENT_DIR)/documents/THINGS_TO_DO.txt \
                $(CURRENT_DIR)/Makefile \
		$(CURRENT_DIR)/make.system \
		$(CURRENT_DIR)/COPYING  

