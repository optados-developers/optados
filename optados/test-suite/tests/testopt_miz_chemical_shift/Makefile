all: $(patsubst %.cell.bz2, %.cell, $(wildcard *.cell.bz2)) $(patsubst %.castep.bz2, %.castep, $(wildcard *.castep.bz2)) $(patsubst %.param.bz2, %.param, $(wildcard *.param.bz2))

%.cell: %.cell.bz2
	$(eval SEEDNAME:=$(patsubst %.cell.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).cell

%.param: %.param.bz2
	$(eval SEEDNAME:=$(patsubst %.param.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).param

%.castep: %.castep.bz2
	$(eval SEEDNAME:=$(patsubst %.castep.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).castep

.PHONY: all
