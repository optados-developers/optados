#!/usr/bin/python
import os,sys
import matplotlib
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
            pdis_data.append({'name':proj_string,
                              'x':[],'y':[],'i':[]})

    #FIXME Header is now parsed, but the data still needs to be read in
    # currently have the data as a list so for each line at each k-point
    # y i_0 i_1 --> pdis_data[0]['y'] etc...

    #WARNING FIXME reading this into a list is a bit dangerous as lists are
    #changeable later on...so need to either be extra careful not to move around 
    #items in list or just read in the file and then be done with it
    print(pdis_data)
    return pdis_data
    
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
    
    pdis_data = read_pdis(pdisfile)




#############
#Main script#
#############
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Generates a quick plot for analysing outpu from a \
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
