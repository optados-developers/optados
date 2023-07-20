all: $(patsubst %.bands.bz2, %.bands, $(wildcard *.bands.bz2)) $(patsubst %.dome_fmt.bz2, %.dome_bin, $(wildcard *.dome_fmt.bz2)) $(patsubst %.fem_fmt.bz2, %.fem_bin, $(wildcard *.fem_fmt.bz2)) $(patsubst %.ome_fmt.bz2, %.ome_bin, $(wildcard *.ome_fmt.bz2)) $(patsubst %.pdos_fmt.bz2, %.pdos_bin, $(wildcard *.pdos_fmt.bz2)) $(patsubst %=C2.odi.bz2, %.odi, $(wildcard *=C2.odi.bz2))

%.dome_bin:  %.dome_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.dome_fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).dome_fmt && ../../../od2od -i dome_fmt -o dome_bin  $(SEEDNAME) && rm $(SEEDNAME).dome_fmt

%.fem_bin:  %.fem_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.fem_fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).fem_fmt && ../../../od2od -i fem_fmt -o fem_bin  $(SEEDNAME) && rm $(SEEDNAME).fem_fmt

%.ome_bin:  %.ome_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.ome_fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).ome_fmt && ../../../od2od -i ome_fmt -o ome_bin  $(SEEDNAME) && rm $(SEEDNAME).ome_fmt

%.pdos_bin:  %.pdos_fmt.bz2 %.bands
	$(eval SEEDNAME:=$(patsubst %.pdos_fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).pdos_fmt && ../../../od2od -i pdos_fmt -o pdos_bin  $(SEEDNAME) && rm $(SEEDNAME).pdos_fmt && rm $(SEEDNAME).log

%.bands: %.bands.bz2
	$(eval SEEDNAME:=$(patsubst %.bands.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).bands

.PHONY: all
