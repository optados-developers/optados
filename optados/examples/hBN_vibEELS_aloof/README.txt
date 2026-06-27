########################################
#   Aloof vib-EELS example, for hBN    #
########################################

# NOTE: It is assumed that the 'odplot.py' plotting utility (found     
#       in optados/python) is available on PATH.

(1) Run the example with:

        optados.x hBN

    This will produce various .dat files, as well as hBN.odo

(2) Plot the (normalized) aloof-loss spectrum, and oscillator strength tensor, with:

        odplot.py aloof_ost

(3) Other interesting quantities can be plotted, such as the Polarizability (alpha):

        odplot.py alpha
