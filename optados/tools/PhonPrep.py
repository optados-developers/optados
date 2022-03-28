#!/usr/bin/python3

import os
import argparse
import sys
sys.tracebacklimit = 0
############################################################################################
# This module produces a .chge_trans and .adf file. The former contains a single           #
# number per atom in the unit cell describing the static/dynamic charge (from an efield    #
# calculation), whilst the latter contains (for every input k-point and temperature) the   #
# atomic displacement parameters from a thermodynamics calculation.                        #
############################################################################################
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# Fill the charge_dic dictionary containing atom names and their static/dynamic charges

def fill_charge_dic(castep_file,charge_type,outfile): 
  data=[]
  file_in=open(castep_file)
  file_in=file_in.readlines()
  for line in file_in:
    line=line.rstrip("\n").split(" ")
    line=list(filter(None,line))
    data.append(line)
  data=list(filter(None,data))

  # Check it is a .castep file and contains total number of ions in file.
  for line in data:
    if line[0]=="Total" and line[1]=="number" and line[3]=="ions":
      number_of_species = int(line[-1])
      break
  try:
    number_of_species
  except NameError:
    print("Cannot find total number of species, please input a CASTEP output file!")

  # Fill charge_dic with Born, Mulliken, or Hirshfield charges
  if charge_type == "born":
    for line in data:
      if len(line)==3 and line[0]=="Born" and line[1]=="Effective" and line[2]=="Charges":
        Born_index=data.index(line)+2
        break
    try:
      Born_index
    except NameError:
      print("Cannot find Born effective charges, please perform the appropriate calculations (CALCULATE_BORN_CHARGES = TRUE).")

    for atom_ind in range(number_of_species): 
      atom=data[Born_index+atom_ind*3][0]+data[Born_index+atom_ind*3][1]
      row1=[float(i) for i in [data[Born_index+atom_ind*3][2],data[Born_index+atom_ind*3][3],data[Born_index+atom_ind*3][4]]]
      row2=[float(i) for i in [data[Born_index+atom_ind*3+1][0],data[Born_index+atom_ind*3+1][1],data[Born_index+atom_ind*3+1][2]]]
      row3=[float(i) for i in [data[Born_index+atom_ind*3+2][0],data[Born_index+atom_ind*3+2][1],data[Born_index+atom_ind*3+2][2]]]
      bornav = (1./3.)*(row1[0]+row2[1]+row3[2])
      charge_dic["{}".format(atom)]=bornav

  elif charge_type == "mulliken":
    for line in data:
      if len(line)==3 and line[0]=="Atomic" and line[1]=="Populations" and line[2]=="(Mulliken)":
        Mulliken_index = data.index(line)+4
        break
    try:
      Mulliken_index
    except NameError:
      print("Cannot find Mulliken charges, please perform the appropriate calculations.")

    for atom_ind in range(number_of_species):
      atom=data[Mulliken_index+atom_ind][0]+data[Mulliken_index+atom_ind][1]
      mull = data[Mulliken_index+atom_ind][-1]
      charge_dic["{}".format(atom)]=mull

  elif charge_type == "hirshfeld":
    for line in data:
      if len(line)==2 and line[0]=="Hirshfeld" and line[1]=="Analysis":
        Hirsh_index = data.index(line)+4
        break
    try:
      Hirsh_index
    except NameError:
      print("Cannot find Hirshfeld charges, please perform the appropriate calculations. (CALCULATE_HIRSHFELD : TRUE)")

    for atom_ind in range(number_of_species):
      atom=data[Hirsh_index+atom_ind][0]+data[Hirsh_index+atom_ind][1]
      hirsh = data[Hirsh_index+atom_ind][-1]
      charge_dic["{}".format(atom)]=hirsh

  else:
    raise NameError("Please enter a valid type of charge input (mulliken, hirshfeld or born).")

  f=open(outfile,'w')
  f.write("# This file contains the atomic charge analysis (mulliken, hirshfeld or born). The born charge is 1/3*Trace(Born matrix).\n")
  f.write("Number of atoms: {}\n".format(number_of_species))
  f.write("Type of charge: {}\n".format(charge_type))
  for atom in charge_dic.keys():
    f.write("{:>5} {:>5} {:>10}\n".format(mysplit(atom)[0], mysplit(atom)[1], charge_dic[atom]))
  f.close()

# Fill the adf_dic dictionary containing atom names and their thermal factors at all calculated temperatures
def fill_adf_dic(castep_file,outfile,temperature_val):
  data=[]
  file_in=open(castep_file)
  file_in=file_in.readlines()
  for line in file_in:
    line=line.rstrip("\n").split(" ")
    line=list(filter(None,line))
    data.append(line)
  data=list(filter(None,data))

  for line in data:
    if len(line)==8 and line[0]+line[1]+line[2]+line[3]=="Totalnumberofions":
      number_of_species = int(line[-1])
      break
  try:
    number_of_species
  except NameError:
    print("Cannot find total number of species, please input a CASTEP output file!")

  for line in data:
    if len(line)==4 and line[0]=="Atomic" and line[1]=="Displacement" and line[3]=="(A**2)":
      adf_index = data.index(line)+5
      # Uii_labels = data[data.index(line)+3]
      break
  try:
    adf_index
  except NameError:
    print("Cannot find Atomic Displacement Parameters, please perform the appropriate thermodynamics calculation.")

  #Assuming the last try is passed, so thermodynamics calculation is performed:
  for line in data:
    if len(line)==6 and line[0]=="Number" and line[2]=="temperature" and line[3]=="values":
      number_of_temp_vals = int(line[-1])
      break

  temperature_list=[]
  for Uii_line_index in range(number_of_species):
    atom = data[adf_index+Uii_line_index][1]+data[adf_index+Uii_line_index][2]
    adf_dic["{}".format(atom)]={}
    temperature=float(data[adf_index+Uii_line_index][0])
    temperature_list.append(temperature)
    adf_dic["{}".format(atom)]["{}".format(temperature)]=data[adf_index+Uii_line_index][3:]
  for Uii_tempblock_index in range(1,number_of_temp_vals):
    for Uii_atomblock_index in range(number_of_species):
      atom = data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][1]+data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][2]
      temperature=float(data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][0])
      temperature_list.append(temperature)
      adf_dic["{}".format(atom)]["{}".format(temperature)]=data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][3:]

  temperature_list=list(set(temperature_list))
  
  if len(temperature_list)>1 and args.temperature==None:
    raise ValueError("The input thermodynamics file contains {} temperature values ({}), please choose a specific one with the -t option.".format(len(temperature_list),temperature_list))
  elif len(temperature_list)==1 and args.temperature==None:
    f=open(outfile,'w')
    f.write(" # This file contains the atomic displacement parameters or thermal factors used to calculate the debye-waller factor \n")
    f.write(" # Each line contains 6 numbers in the order of U11, U22, U33, U23 U31 U12 (since the matrix is symmetric).\n")
    f.write("Number of atoms: {} \n".format(number_of_species))
    f.write("Temperature: {}\n".format(temperature_list[0]))
    for atom in adf_dic.keys():
      f.write("{:>5} {:>5}".format(mysplit(atom)[0],mysplit(atom)[1] + " "))
      for Uii in range(6):
        f.write("{:>10}".format(adf_dic[atom]["{}".format(temperature_list[0])][Uii]))
      f.write("\n")
    f.close()
  else:
    if float(args.temperature) in temperature_list:
      f=open(outfile,'w')
      f.write(" # This file contains the atomic displacement parameters or thermal factors used to calculate the debye-waller factor \n")
      f.write(" # Each line of 6 numbers is the order of U11, U22, U33, U23 U31 U12, since the matrix is symmetry.\n")
      f.write("Number of atoms: {} \n".format(number_of_species))
      f.write("Temperature: {}\n".format(float(args.temperature)))
      for atom in adf_dic.keys():
        f.write("{:>5}".format(atom + " "))
        for Uii in range(6):
          f.write("{:>10}".format(adf_dic[atom]["{}".format(float(args.temperature))][Uii]))
        f.write("\n")
      f.close()
    else:
      raise ValueError("Temperature specified does not exist in thermodynamics calculation.")

# Function to split text into characters vs numbers (e.g. Atom H1 to H 1)
def mysplit(text):
  atom_name = text.rstrip('01234567890')
  atom_number = text[len(atom_name):]
  return atom_name, atom_number

if __name__ == '__main__':
  # parser = argparse.ArgumentParser(description='This script can produce a .chge_trans from an efield calculation containing \
  # mulliken/born/hirshfeld charges and/or a .adf file from a thermodynamics calcilation containing \
  # atomic displacement parameters/thermal factors.')
  parser = MyParser(description="Produces a .chge_trans from an efield calculation containing \
  mulliken/born/hirshfeld charges and/or produces a .adf file from a thermodynamics calcilation containing \
  atomic displacement parameters/thermal factors.")
  parser.add_argument('-c', '--charge', type=str, nargs=2, help='Type of charge wanted from efield calculation. Please enter \
  two arguments: mulliken, hirshfeld or born followed by the path to the .castep file.')
  parser.add_argument('-a', '--adf', type=str, help='Requesting Atomic displacement factors require only the .castep file path.')
  parser.add_argument('-t', '--temperature', type=float, help='The temperature at which thermodynamics calculations are performed. Please select one if multiple temeratures are considered.')
  args = parser.parse_args()
  charge_dic = {} # Dictionary for charge analysis storage to produce .chge_trans file
  adf_dic = {} #Dictionary for thermal analysis storage to produce .adf file


  if args.charge != None:
    try:
      args.charge[0] in ['born','mulliken','hirshfeld']
    except ValueError:
      print("Charge type {} does not exist, please choose between hirshfeld, mulliken and born.")
    try:
      f=open(args.charge[1],'r')
    except FileNotFoundError:
      print("Cannot open {}, file does not exist.".format(args.charge[1]))
    
    chargetype=args.charge[0]
    castepfile_chgetrans=args.charge[1]
    outfile_chgetrans=castepfile_chgetrans.split(".")[0]+".chge_trans"
    fill_charge_dic(castepfile_chgetrans,chargetype,outfile_chgetrans)
  if args.adf != None:
    try:
      f=open(args.adf,'r')
    except FileNotFoundError:
      print("Cannot open {}, file does not exist.".format(args.charge[1]))
    
    castepfile_adf=args.adf 
    outfile_adf=castepfile_adf.split(".")[0]+".adf"
    fill_adf_dic(castepfile_adf,outfile_adf,args.temperature)

  if args.charge==None and args.adf == None:
    parser.error("No input requests, exiting.")


