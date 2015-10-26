import OptaPyDOS as opd
# OptaPyDOS Tools
# ===============
# For doing stuff with OptaDOS in Python
# A J Morris 26/10/2015
od_OptaPyDOS_tools_version = 1.0


def od_setup_from_file(seedname):
    """ Read an odi file, write headers and read minimun mencessary
    castep files then distribute to nodes.
    """
    import OptaPyDOS as opd	
    print " OptaPyDOS Tools -- Setup Calculation from Seedname"
    opd.od_io.stdout = opd.od_io.io_file_unit()
    opd.od_io.stderr = opd.od_io.io_file_unit()
    opd.od_comms.comms_setup()
    opd.od_io.seedname=seedname
    print " Seedname:", opd.od_io.seedname
    opd.od_parameters.param_read()
    opd.od_parameters.param_write_header()
    opd.od_parameters.param_write()
    opd.od_electronic.elec_read_band_energy()
    opd.od_cell.cell_calc_lattice()
    opd.od_cell.cell_report_parameters()
    opd.od_electronic.elec_report_parameters()
    opd.od_parameters.param_dist()
    opd.od_cell.cell_dist()

    
def od_set_efermi():
    """ Deal with the Fermi energy

    Call the internal optados routines and shift the energy scale if it is
    required. Note that the OptaDOS code only does the shifting just
    before it prints.
    """
    import OptaPyDOS as opd	
    opd.od_dos_utils.dos_utils_set_efermi() # Decide which efermi to use
    chemical_potential = opd.od_electronic.efermi

    # If the user wanted the fermi level setting to zero, do it here.
    if opd.od_parameters.set_efermi_zero:
        opd.od_dos_utils.e = [ x-chemical_potential for x in opd.od_dos_utils.e]
        chemical_potential=0

    return chemical_potential
   
