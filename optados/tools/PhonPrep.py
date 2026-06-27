"""
J. A. J. Whaley-Baldwin, June 2026

Reads partial charges, atomic displacement parameters and Born effective charge tensors from a
CASTEP .phonon file, and prints the result to stdout, in a format compatible for use with the
vib-EELS functionality of OptaDOS.

The output of this script can be placed directly into an OptaDOS .odd file, for use with vib-EELS.

Usage:

  python3 PhonPrep.py castep_file

Note: Original PhonPrep.py written by R. J. Nicholls, which produced .chge_trans and .adf files.

      This script has been rewritten from scratch to instead produce output that is compliant with
      the data blocks within the new OptaDOS .odd file format.
"""

import sys
import numpy as np

castep_file = sys.argv[1]

file = open("./"+castep_file,"r")
lines = file.readlines()

def read_partial_charges():
	count = 1
	read_data = False
	for line in lines:
		if ("Atomic Populations (Mulliken)" in line):
			read_data = True
			continue
		if ("Species" in line):
			continue
		if ("=====" in line):
			if (count > 1):
				break
			continue
		if (read_data):
			line = line.replace("\n","")
			line = [ x for x in line.split(" ") if x ]
			if ( len(line) == 8 ):
				print( "%d  %.6f" % (count,float(line[-1])) )
				count += 1

def read_born_eff_charges():
	bec_tensor = []
	bec_tensors = []
	atoms = []
	idxs = []
	read_data = False
	count = 1
	for line in lines:
		if ("Born Effective Charges" in line):
			read_data = True
			continue
		if ("-----" in line):
			continue
		if (read_data):
			if ( ("=====" in line) and (count > 1) ):
				break
			line = line.replace("\n","")
			line = [ x for x in line.split(" ") if x ]
			if ( len(line) == 5 ):
				atoms.append(line[0])
				idxs.append(int(line[1]))
				if ( len(bec_tensor) > 0 ):
					bec_tensors.append(bec_tensor)
				bec_tensor = []
				bec_tensor.append( [float(line[2]),float(line[3]),float(line[4])] )
				count += 1
				continue
			elif ( len(line) == 3 ):
				bec_tensor.append( [float(line[0]),float(line[1]),float(line[2])] )
				continue

	bec_tensors.append(bec_tensor)

	bec_tensors = np.array(bec_tensors)
	#print(bec_tensors.shape)

	for i,atm in enumerate(atoms):
		print(i+1)
		print( "%.8f  %.8f  %.8f" % (bec_tensors[i][0][0],bec_tensors[i][0][1],bec_tensors[i][0][2]) )
		print( "%.8f  %.8f  %.8f" % (bec_tensors[i][1][0],bec_tensors[i][1][1],bec_tensors[i][1][2]) )
		print( "%.8f  %.8f  %.8f" % (bec_tensors[i][2][0],bec_tensors[i][2][1],bec_tensors[i][2][2]) )
		print("")

def read_adps(target_temp=300):
	print(target_temp)
	adps = []
	count = 1
	read_data = False
	for line in lines:
		if ( ("U11" in line) and ("U22" in line) ):
			read_data = True
			continue
		if ("-----" in line):
			if (count > 1):
				break
			continue
		if (read_data):
			line = line.replace("\n","")
			line = [ x for x in line.split(" ") if x ]
			T = float(line[0])
			if ( abs(T-target_temp) < 0.1 ):
				U11 = float(line[3])
				U22 = float(line[3])
				U33 = float(line[3])
				U23 = float(line[3])
				U31 = float(line[3])
				U12 = float(line[3])
				#adps.append([U11,U22,U33,U23,U31,U12])
				print( "%d  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f" % (count,U11,U22,U33,U23,U31,U12) )
				count += 1
			else:
				continue

################################################################################################################
# Run.

# Parse ADPs.
try:
	print("BEGIN ATOMIC_DISPLACEMENT_PARAMETERS")
	read_adps()
	print("END ATOMIC_DISPLACEMENT_PARAMETERS")
	print("\n")
except:
	print("Could not parse Atomic Displacement Parameters from .castep file")

# Parse partial charges.
try:
	print("BEGIN PARTIAL_CHARGES")
	read_partial_charges()
	print("END PARTIAL_CHARGES")
	print("\n")
except:
	print("Could not parse partial charges from .castep file")

# Parse BEC tensor.
try:
	print("BEGIN BORN_EFF_CHARGE_TENSOR")
	read_born_eff_charges()
	print("END BORN_EFF_CHARGE_TENSOR")
except:
	print("Could not parse Born effective charges from .castep file")
