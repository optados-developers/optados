"""
Parser function parse() to parse the .wout output file of Wannier90.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

def parse(fname):
    """
    Open the file, parses it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        for line in f:
            pass
        last_line = line

    elements = last_line.split('   ')

    if len(elements)==5:
        retdict["energy"].append(float(elements[0]))
        retdict["up_dos"].append(float(elements[1]))
        retdict["down_dos"].append(float(elements[2]))
        retdict["up_intdos"].append(float(elements[3]))
        retdict["down_intdos"].append(float(elements[4]))
    elif len(elements)==3:
        retdict["energy"].append(float(elements[0]))
        retdict["dos"].append(float(elements[1]))
        retdict["intdos"].append(float(elements[2]))
    else:
        sys.exit("Error reading the number of entries at end of"+fname)

    ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
