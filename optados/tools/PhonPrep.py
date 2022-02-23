import sys
import numpy
############################################################################################
# This module produces a .chge_trans and .adf file. The former contains a single           #
# number per atom in the unit cell describing the static/dynamic charge (from an efield    #
# calculation), whilst the latter contains (for every input k-point and temperature) the   #
# atomic displacement parameters from a thermodynamics calculation.                        #
############################################################################################

outputtype=sys.argv[1] #Type of calculation wanted (-c or -a)
if outputtype=="-c":
  chargetype=sys.argv[2] #Type of charge (mulliken, hirschfield or born) requested by user.
  castepfile=sys.argv[3] #CASTEP file containing efield calculations
  outfile = castepfile.split(".")[0]+".chge_trans"
elif outputtype=="-a":
  castepfile=sys.argv[2] #CASTEP file containing efield calculations
  outfile = castepfile.split(".")[0]+".adf"
else:
  raise NameError("Please choose the type of calculation as -c or -a")

charge_dic = {} # Dictionary for charge analysis storage to produce .chge_trans file
adf_dic = {} #Dictionary for thermal analysis storage to produce .adf file

# Fill the charge_dic dictionary containing atom names and their static/dynamic charges
def fill_charge_dic(castep_file,charge_type): 
  data=[]
  file_in=open(castepfile)
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

    for atom_ind in numpy.arange(number_of_species): 
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

    for atom_ind in numpy.arange(number_of_species):
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

    for atom_ind in numpy.arange(number_of_species):
      atom=data[Hirsh_index+atom_ind][0]+data[Hirsh_index+atom_ind][1]
      hirsh = data[Hirsh_index+atom_ind][-1]
      charge_dic["{}".format(atom)]=hirsh

  else:
    raise NameError("Please enter a valid type of charge input (mulliken, hirshfeld or born")

  f=open(outfile,'w')
  f.write("# This file contains the atomic charge analysis (mulliken, hirshfeld or born). The born charge is 1/3*Trance(Born matrix).\n")
  f.write("Number of atoms: {}\n".format(number_of_species))
  f.write("Type of charge: {}\n".format(charge_type))
  for atom in charge_dic.keys():
    f.write("{} {}\n".format(atom, charge_dic[atom]))
  f.close()


def fill_adf_dic(castep_file):
  data=[]
  file_in=open(castepfile)
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
  for Uii_line_index in numpy.arange(number_of_species):
    atom = data[adf_index+Uii_line_index][1]+data[adf_index+Uii_line_index][2]
    adf_dic["{}".format(atom)]={}
    temperature=float(data[adf_index+Uii_line_index][0])
    temperature_list.append(temperature)
    adf_dic["{}".format(atom)]["{}".format(temperature)]=data[adf_index+Uii_line_index][3:]
  for Uii_tempblock_index in numpy.arange(1,number_of_temp_vals):
    for Uii_atomblock_index in numpy.arange(number_of_species):
      atom = data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][1]+data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][2]
      temperature=float(data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][0])
      temperature_list.append(temperature)
      adf_dic["{}".format(atom)]["{}".format(temperature)]=data[adf_index+Uii_tempblock_index*number_of_species+Uii_atomblock_index][3:]

  temperature_list=list(set(temperature_list))
  f=open(outfile,'w')
  f.write(" # This file contains the atomic displacement parameters or thermal factors used to calculate the debye-waller factor \n")
  f.write(" # Each line of 6 numbers is the order of U11, U22, U33, U23 U31 U12, since the matrix is symmetry.")
  f.write("Number of atoms: {} \n".format(number_of_species))
  f.write("Temperature List: ")
  for temp in temperature_list:
    f.write(" {}".format(temp))
  f.write("\n")
  for temp in temperature_list:
    f.write("Temperature: {} \n".format(temp))
    for atom in adf_dic.keys():
      f.write(atom + " ")
      for Uii in numpy.arange(6):
        f.write("{} ".format(adf_dic[atom]["{}".format(temp)][Uii]))
      f.write("\n")
  f.close()

if outputtype=="-c":
  fill_charge_dic(castepfile,chargetype)
elif outputtype=="-a":
  fill_adf_dic(castepfile)
