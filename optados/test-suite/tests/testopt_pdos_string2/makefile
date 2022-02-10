all: $(patsubst %.dome_fmt.bz2, %.dome_bin, $(wildcard *.dome_fmt.bz2)) $(patsubst %.bands.bz2, %.bands, $(wildcard *.bands.bz2)) $(patsubst %.pdos_fmt.bz2, %.pdos_bin, $(wildcard *.pdos_fmt.bz2))

%.dome_bin:  %.dome_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.dome_fmt.bz2, %, $<))
	echo "DOME"	
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).dome_fmt && ../../../od2od -i dome_fmt -o dome_bin  $(SEEDNAME) && rm $(SEEDNAME).dome_fmt

%.pdos_bin:  %.pdos_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.pdos_fmt.bz2, %, $<))
	echo "PDOS"
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).pdos_fmt && ../../../od2od -i pdos_fmt -o pdos_bin  $(SEEDNAME) && rm $(SEEDNAME).pdos_fmt

%.bands: %.bands.bz2
	$(eval SEEDNAME:=$(patsubst %.bands.bz2, %, $<))
	echo "BANDS"
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).bands

.PHONY: all
