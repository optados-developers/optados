"""
Parser function parse() to parse the .log file of od2od
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

nbands = re.compile("Number\ of\ bands\ :\s*([0-9\.-]+)\s*")

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
        match = nbands.search(l)
        if match:
            retdict["nbands"].append(float(match.groups()[0]))
            continue
        ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
