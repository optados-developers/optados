"""
Parser function parse() to parse the .odo file of OptaDOS.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

e_fermi_fb = re.compile("Fermi\ energy\ \(Fixed\ broadening\)\ \:\s*([0-9\.-]+)\s*")
e_fermi_ab = re.compile("Fermi\ energy\ \(Adaptive\ broadening\)\ \:\s*([0-9\.-]+)\s*")
e_fermi_lb = re.compile("Fermi\ energy\ \(Linear\ broadening\)\ \:\s*([0-9\.-]+)\s*")

qe_bulk = re.compile("Bulk \s*([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s*")
qe_total = re.compile("Total Quantum Efficiency \(electrons/photon\):\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s*") 
mte = re.compile("Weighted\ Mean\ Transverse\ Energy\ \(eV\)\:\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s*")

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
        match = qe_bulk.search(l)
        if match:
            retdict['qe_bulk'].append(float(match.groups()[0]+'E'+match.groups()[1]))
            continue
        match = qe_total.search(l)
        if match:
            retdict['qe_total'].append(float(match.groups()[0]+'E'+match.groups()[1]))
            continue
        match = mte.search(l)
        if match:
            retdict['mte'].append(float(match.groups()[0]+'E'+match.groups()[1]))
            continue
        ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
