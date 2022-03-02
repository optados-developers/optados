"""
Parser function parse() to parse general columnar data optput files of OptaDOS and od2od
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

    count = len(open(fname).readlines(  ))

    modulo_line = count//10

    line_counter=0
    with open(fname) as f:

        for line in f:
            # Let's only grab every 10 lines, no point going overboard
            if (line_counter % modulo_line) == 0:
                elements = line.split()
                # If there's only one element it's probably a header so lets ignore it.
                if len(elements)>1:
                    try:
                        # Maybe the header isn't a float if so, let's ignore it, and pick
                        # it up again when we're out of the header
                        retdict["l"+str(line_counter)+"_c2"].append(float(elements[1]))
                    except:
                        pass
                # Let's try another column too, if it's there.
                if len(elements)>4:
                    try:
                        retdict["l"+str(line_counter)+"_c5"].append(float(elements[4]))
                    except:
                        pass
            line_counter=line_counter+1

    ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
