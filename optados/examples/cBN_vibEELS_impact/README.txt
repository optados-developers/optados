########################################
#   Impact vib-EELS example, for cBN   #
########################################

# NOTE: It is assumed that the 'odplot.py' plotting utility (found     
#       in optados/python) is available on PATH.

(1) Run the example with:

        optados.x cBN

    This will produce various .dat files, as well as cBN.odo

(2) Plot the vib-EELS impact intensity on top of the phonon dispersion, with:

        odplot.py impact_intensity -hsp_labels G X UbK G L W X

    (the -hsp_labels argument is optional, but if it is omitted, high-symmetry points will not be labelled).

(3) Plot a 'heatmap' of the vib-EELS impact intensity, with:

        odplot.py impact_heatmap -hsp_labels G X UbK G L W X  -elim 0 200  -qslices 0.15 0.40 0.82

    (-elim is optional, and simply adjusts the energy limits of the plot)

    (-qslices is also optional; if supplied, the heatmap will be vertically sliced at various fractions along the q-path)
