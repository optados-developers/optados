"""
Parser function parse() to parse the .wout output file of Wannier90.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

e_fermi_fb = re.compile("Fermi\ energy\ \(Fixed\ broadening\)\ \:\s*([0-9\.-]+)\s*")
e_fermi_ab = re.compile("Fermi\ energy\ \(Adaptive\ broadening\)\ \:\s*([0-9\.-]+)\s*")
e_fermi_lb = re.compile("Fermi\ energy\ \(Linear\ broadening\)\ \:\s*([0-9\.-]+)\s*")


def parse(fname):
    """
    Open the file, parses it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for lno, l in enumerate(lines):

        # various Omegas (four numbers, typically at the end)
        match = e_fermi_fb.search(l)
        if match:
            retdict["fermi_fb"].append(float(match.groups()[0]))
            continue
        match = e_fermi_ab.search(l)
        if match:
            retdict["fermi_ab"].append(float(match.groups()[0]))
            continue
        match = e_fermi_lb.search(l)
        if match:
            retdict["fermi_lb"].append(float(match.groups()[0]))
            continue

        ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
