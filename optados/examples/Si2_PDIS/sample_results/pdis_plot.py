#!/usr/bin/python
import os,sys
import matplotlib.pyplot as plt
import random
import argparse


def read_pdis(pdisfile):

    pdis_data = []

    # Need to determine what projectors there are the options
    # are species, species_ang, ang, or sites
    proj_num = 0
    
    # scrape out the header data to determine what is being
    # projected
    for i,line in enumerate(pdisfile):
        if line.startswith('#') and 'Projector:' in line:
            proj_num += 1
            j = i+2
            projline = pdisfile[j]
            
            # keep track of what element, site number, and ang channel is in projector
            elem = []
            site = []
            ang = []
            while '+' not in projline:
                strippedline = projline.strip().split('|')[1].strip().split(' ') 
                nospaces = []
                for ele in strippedline:
                    if ele.strip():
                        nospaces.append((ele))
                if nospaces[0] not in elem:
                    elem.append(nospaces[0])
                if nospaces[1] not in site:
                    site.append(nospaces[1])
                if nospaces[2] not in ang:
                    ang.append(nospaces[2])
                j += 1
                projline = pdisfile[j]
            
            if len(elem) > 1:
                elem = []
            if len(site) > 1:
                site = []
            if len(ang) > 1:
                ang = []
                
            # Get Al1s, or Sip, or p, as the name of projector
            proj_string = elem+site+ang
            proj_string = ''.join(proj_string)
            pdis_data.append({'position':proj_num,'name':proj_string,
                              'kpoints':[],'energies':[],'intensities':[]})

    k_point = 0
    k_point_list = []
    for i,line in enumerate(pdisfile):
        if line.startswith('K-point'):
            newline = line.strip().split(' ') 
            nospaces = []
            for ele in newline:
                if ele.strip():
                    nospaces.append(ele)
            k_point_list.append([round(float(point),2) for point in nospaces[2:]])
            for proj in pdis_data:
                proj['energies'].append([])
                proj['intensities'].append([])

            j = i + 1
            # Only read the data for that K-point and then stop
            k = 0 #counter for each energy level in the k-point

            while j < len(pdisfile) and 'K-point' not in pdisfile[j]:
                nextline = pdisfile[j]
                nextline = nextline.strip().split(' ')
                nospaces = []
                for ele in nextline:
                    if ele.strip():
                        nospaces.append((float(ele)))
                

                for proj in pdis_data:
                    proj['energies'][k_point].append(nospaces[0])
                    proj['intensities'][k_point].append(nospaces[proj['position']])
                j += 1
                k += 1

            k_point += 1

    for proj in pdis_data:
        proj['kpoints'] = [i for i in range(len(proj['energies']))]
    return pdis_data,k_point_list
    
def plot_pdis(infile):
    """Plot the projected dispersion from a pdis.dat file 
    """

    # First open the .pdis.dat file in <infile>.pdis.dat
    try:
        pdisfile = open('%s.pdis.dat'%infile).readlines()
    except:
        print(('\n%s.pdis.dat not found in folder'%infile))
        print(('\n...EXITING...'))
        exit()
    
    pdis_data,k_point_list = read_pdis(pdisfile)

    
    # Now plot the resulting data
    fig, ax = plt.subplots(1,1,figsize=(10,8))

    #Here make the list of colors
    #Possible FIXME is for users to add in their own custom colors
    color = [ "red", "blue", "green", "yellow", "purple", "orange", "magenta", "black" ]


    
    for proj in pdis_data:
        label = False
        for x,y,i in zip(proj['kpoints'],proj['energies'],proj['intensities']):
            i = [num*50 for num in i]
            if not label:
                plt.scatter([x]*len(y),y,s=i,color=color[proj['position']-1],label=proj['name'])
                label = True
            else:
                plt.scatter([x]*len(y),y,s=i,color=color[proj['position']-1])
    
    k_point_list_str = []
    for kpoint in k_point_list:
        if kpoint == [0.0,0.0,0.0]:
            k_point_list_str.append(r'$\Gamma$')
        else:
            k_point_list_str.append(None)

    for i,kpoint in enumerate(k_point_list_str):
        if kpoint is not None:
            ax.plot([pdis_data[0]['kpoints'][i],pdis_data[0]['kpoints'][i]],[min(pdis_data[0]['energies'][0]),max(pdis_data[0]['energies'][0])],color='k',linestyle='--')

    ax.tick_params(which='major', length=0)
    ax.set_xticks(pdis_data[0]['kpoints'])
    ax.set_xticklabels(k_point_list_str)
    plt.legend()
    plt.show()



#############
#Main script#
#############
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Generates a quick plot for analysing output from a \
            pdis OptaDOS calculation. 
            Example usage is `pdis_plotter.py -i <filename>\
            ''')


    #Options
    parser.add_argument('-i', '--inputfile', type=str,
                        help='Name of the input file (<input>.pdis.dat only) from OptaDOS')

    args = parser.parse_args()

    if args.inputfile is not None:
        # remove the .pdis.dat ending if the user has supplied it
        if args.inputfile.endswith('.pdis.dat'):
            args.inputfile = args.inputfile.split('.pdis.dat')[0]
    else:
        # look for a .pdis.dat file
        for file in os.listdir():
            if file.endswith('.pdis.dat'):
                args.inputfile = file.split('.pdis.dat')[0]

    plot_pdis(args.inputfile)
