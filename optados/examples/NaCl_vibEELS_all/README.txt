#################################################
#   Aloof & Impact vib-EELS example, for NaCl   #
#################################################

# NOTE: It is assumed that the 'odplot.py' plotting utility (found
#       in optados/python) is available on PATH.

(1) Run the example with:

        optados.x NaCl

    This will produce various .dat files, as well as NaCl.odo

(2) Plot the vib-EELS impact intensity on top of the phonon dispersion, with:

        odplot.py impact_intensity -hsp_labels G X M R G

    (the -hsp_labels argument is optional, but if it is omitted, high-symmetry points will not be labelled).

(3) Plot a 'heatmap' of the vib-EELS impact intensity, with:

        odplot.py impact_heatmap -hsp_labels G X M R G  -elim 0 40  -qslices 0.25 0.55 0.78

    (-elim is optional, and simply adjusts the energy limits of the plot)

    (-qslices is also optional; if supplied, the heatmap will be vertically sliced at various fractions along the q-path)

(4) Plot the (normalized) aloof-loss spectrum, and oscillator strength tensor, with:

        odplot.py aloof_ost

(5) Other interesting quantities can be plotted, such as the Polarizability (alpha):

        odplot.py alpha
