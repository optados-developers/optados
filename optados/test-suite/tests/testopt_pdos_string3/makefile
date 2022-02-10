all: $(patsubst %.ome_fmt.bz2, %.ome_bin, $(wildcard *.ome_fmt.bz2)) $(patsubst %.bands.bz2, %.bands, $(wildcard *.bands.bz2)) $(patsubst %.pdos_fmt.bz2, %.pdos_bin, $(wildcard *.pdos_fmt.bz2))

%.ome_bin:  %.ome_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.ome_fmt.bz2, %, $<))
	echo "OME"	
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).ome_fmt && ../../../od2od -i ome_fmt -o ome_bin  $(SEEDNAME) && rm $(SEEDNAME).ome_fmt

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
