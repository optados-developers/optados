"""
Parser function parse() to parse the .odo file of OptaDOS.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

error_msg = re.compile("Error: ([A-Za-z]+)")

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

        match = error_msg.search(l)
        if match:
            retdict['error_msg'].append(match.groups()[0])
            continue
        ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
