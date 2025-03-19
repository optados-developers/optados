"""
Parser function parse() to parse agr output files of OptaDOS
"""
from __future__ import print_function

import inspect
import sys
import re
from collections import defaultdict

from . import show_output


def parse(fname):
    """
    Open the .agr file, parses it and checks for retirms the values.
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    header_end = 0
    nsets = 0
    ndata = 0  # data in a set

    # Regex pattern for various parts of a xmgrace file
    leg_pat = r'@\s*s\d+\s*legend\s*"([^"]*)'  # legend label
    set_pat = r'@target\s*G0\.s(\d+)'

    with open(fname, 'r') as f:
        for lineno, line in enumerate(f):
            s = line.strip()

            if re.search(leg_pat, s):
                # Count the number of sets with a legend.
                # Here we assume this is the number of sets.
                nsets += 1

            if s.startswith('@target'):
                header_end = lineno
                continue
            if s.startswith('@type'):
                continue

            # Count how many data points are in this set.
            # This should be the same in all sets.
            if header_end > 0:
                if s.startswith('&'):
                    break
                ndata += 1

    count = len(open(fname).readlines()[header_end:])
    modulo_line = 10
    lineno = 0

    # Read from the beginning but now skip the header.
    with open(fname, 'r') as f:
        for line in f:
            if lineno > header_end - 1:  # -1 since this includes target
                break
            lineno += 1

        # Now loop around the sets ensuring that we have valid format for a set
        validset = False
        for iset in range(nsets):
            lineno = 0

            # Check that we have a correctly formatted set.
            if re.search(set_pat, line.strip()):
                validset = True
            else:
                # We have an invalid set so set the current set
                # to a silly value and continue to next set.
                retdict[f'set{iset},col1'] = -1*sys.float_info.max
                retdict[f'set{iset},col2'] = -1*sys.float_info.max

            for line in f:
                # Grab every 10 lines in a set ensuring that we do not exceed the set
                if lineno < ndata and (lineno % modulo_line) == 0 and validset:
                    elements = line.split()
                    try:
                        # Grab data from both columns
                        retdict[f'set{iset},col1'].append(float(elements[0]))
                        retdict[f'set{iset},col2'].append(float(elements[1]))
                    except:
                        pass
                else:
                    # Finished reading, let's advance to start of the next set
                    # looking for the ampersand and then reset validset flag.
                    if line.strip().startswith('@'):
                        validset = False
                        break

                lineno += 1

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
