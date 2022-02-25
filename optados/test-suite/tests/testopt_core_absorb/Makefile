all: $(patsubst %.elnes_fmt.bz2, %.elnes_bin, $(wildcard *.elnes_fmt.bz2)) $(patsubst %.dome_fmt.bz2, %.dome_bin, $(wildcard *.dome_fmt.bz2)) $(patsubst %.bands.bz2, %.bands, $(wildcard *.bands.bz2))

%.elnes_bin:  %.elnes_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.elnes_fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).elnes_fmt && ../../../od2od -i elnes_fmt -o elnes_bin  $(SEEDNAME) && rm $(SEEDNAME).elnes_fmt

%.dome_bin:  %.dome_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.dome_fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).dome_fmt && ../../../od2od -i dome_fmt -o dome_bin  $(SEEDNAME) && rm $(SEEDNAME).dome_fmt

%.bands: %.bands.bz2
	$(eval SEEDNAME:=$(patsubst %.bands.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).bands

.PHONY: all
