include make.system

VER := 1.2.380

default : optados

all : optados docs tools

tools: 
	cd ./src && $(MAKE) tools

optados: 
	cd ./src && $(MAKE) optados

# Utility targets
.PHONY: ford clean veryclean 

# These have to be phony as they're recursive makes
.PHONY: tools optados

install:
	 cd ./src && $(MAKE) install

docs:
	cd ./documents && $(MAKE) all

clean:
	cd ./src && $(MAKE) clean 
	cd ./documents && $(MAKE) clean
	cd ./python  && $(MAKE) clean
	cd ./test-suite && $(PYTHON_VER) clean_tests

ford: 
	ford -d ../src -o ../auto_doc ./ford/project.md

python: optados
	cd ./python && $(MAKE) python

veryclean:
	cd ./src && $(MAKE) veryclean
	cd ./documents && $(MAKE) veryclean
	cd ./python &&  $(MAKE) veryclean
	cd ./test-suite && $(PYTHON_VER) clean_tests -i

dist: optados-$(VER).tar.gz

CURRENT_DIR := $(shell basename ${PWD})

optados-$(VER).tar.gz: 
	 cd .. && tar -czf $(CURRENT_DIR)/optados-$(VER).tar.gz \
                $(CURRENT_DIR)/src/*.?90 \
                $(CURRENT_DIR)/src/Makefile \
                $(CURRENT_DIR)/examples/*/*.cell \
                $(CURRENT_DIR)/examples/*/*.param \
                $(CURRENT_DIR)/examples/*/*.recpot \
		$(CURRENT_DIR)/examples/*/*.odi \
                $(CURRENT_DIR)/documents/*.tex \
                $(CURRENT_DIR)/documents/*.bib \
                $(CURRENT_DIR)/documents/*.pdf \
		$(CURRENT_DIR)/documents/images/*.eps \
		$(CURRENT_DIR)/documents/Makefile \
		$(CURRENT_DIR)/README.INSTALL \
		$(CURRENT_DIR)/CHANGE.LOG \
		$(CURRENT_DIR)/Makefile \
		$(CURRENT_DIR)/make.system \
		$(CURRENT_DIR)/COPYING \
                $(CURRENT_DIR)/tools/*.f90 

