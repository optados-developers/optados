import _OptaPyDOS
import f90wrap.runtime
import logging

class Od_Algorithms(f90wrap.runtime.FortranModule):
    """
    Module od_algorithms
    
    
    Defined at ../src/algorithms.f90 lines 35-345
    
    """
    @staticmethod
    def channel_to_am(no):
        """
        channel_to_am = channel_to_am(no)
        
        
        Defined at ../src/algorithms.f90 lines 52-70
        
        Parameters
        ----------
        no : int
        
        Returns
        -------
        channel_to_am : str
        
        """
        channel_to_am = _OptaPyDOS.f90wrap_channel_to_am(no=no)
        return channel_to_am
    
    @staticmethod
    def gaussian(m, w, x):
        """
        gaussian = gaussian(m, w, x)
        
        
        Defined at ../src/algorithms.f90 lines 73-92
        
        Parameters
        ----------
        m : float
        w : float
        x : float
        
        Returns
        -------
        gaussian : float
        
        =========================================================================
         ** Return value of Gaussian(mean=m,width=w) at position x
         I don't know who's this function originally was, CJP? MIJP?
        =========================================================================
        """
        gaussian = _OptaPyDOS.f90wrap_gaussian(m=m, w=w, x=x)
        return gaussian
    
    @staticmethod
    def heap_sort(num_items, weight):
        """
        heap_sort(num_items, weight)
        
        
        Defined at ../src/algorithms.f90 lines 95-170
        
        Parameters
        ----------
        num_items : int
        weight : float array
        
        =========================================================================
         This subroutine sorts the list of weights into descending order.
         On exit, if present, the array of indexes contains the original index
         of each item.
         This is a heap sort
        -------------------------------------------------------------------------
         Arguments:
           num_items (input) :: The number of items to sort
           weight (in/out) :: The weights of each item. On exit these are
                              sorted into descending order.
        -------------------------------------------------------------------------
         Parent module variables used:
           None
        -------------------------------------------------------------------------
         Modules used:
           None
        -------------------------------------------------------------------------
         Key Internal Variables:
           None
        -------------------------------------------------------------------------
         Necessary conditions:
           None
        -------------------------------------------------------------------------
         Written by Chris Pickard 22nd May 2009
        =========================================================================
        """
        _OptaPyDOS.f90wrap_heap_sort(num_items=num_items, weight=weight)
    
    @staticmethod
    def utility_lowercase(string_bn):
        """
        utility_lowercase = utility_lowercase(string_bn)
        
        
        Defined at ../src/algorithms.f90 lines 173-203
        
        Parameters
        ----------
        string_bn : str
        
        Returns
        -------
        utility_lowercase : str
        
        =========================================================================
         Takes a string and converts to lowercase characters
        =========================================================================
        """
        utility_lowercase = _OptaPyDOS.f90wrap_utility_lowercase(string_bn=string_bn)
        return utility_lowercase
    
    @staticmethod
    def utility_frac_to_cart(frac, cart, real_lat):
        """
        utility_frac_to_cart(frac, cart, real_lat)
        
        
        Defined at ../src/algorithms.f90 lines 206-226
        
        Parameters
        ----------
        frac : float array
        cart : float array
        real_lat : float array
        
        ==================================================================
          Convert from fractional to Cartesian coordinates
        ===================================================================
        """
        _OptaPyDOS.f90wrap_utility_frac_to_cart(frac=frac, cart=cart, real_lat=real_lat)
    
    @staticmethod
    def utility_cart_to_frac(cart, frac, recip_lat):
        """
        utility_cart_to_frac(cart, frac, recip_lat)
        
        
        Defined at ../src/algorithms.f90 lines 229-254
        
        Parameters
        ----------
        cart : float array
        frac : float array
        recip_lat : float array
        
        ==================================================================
          Convert from fractional to Cartesian coordinates
        ===================================================================
        """
        _OptaPyDOS.f90wrap_utility_cart_to_frac(cart=cart, frac=frac, \
            recip_lat=recip_lat)
    
    @staticmethod
    def algorithms_erf(x):
        """
        algorithms_erf = algorithms_erf(x)
        
        
        Defined at ../src/algorithms.f90 lines 256-310
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        algorithms_erf : float
        
        """
        algorithms_erf = _OptaPyDOS.f90wrap_algorithms_erf(x=x)
        return algorithms_erf
    
    _dt_array_initialisers = []
    

od_algorithms = Od_Algorithms()

class Od_Cell(f90wrap.runtime.FortranModule):
    """
    Module od_cell
    
    
    Defined at ../src/cell.f90 lines 35-1124
    
    """
    @staticmethod
    def cell_find_mp_grid(kpoints, num_kpts, kpoint_grid_dim, kpoint_offset=None):
        """
        cell_find_mp_grid(kpoints, num_kpts, kpoint_grid_dim[, kpoint_offset])
        
        
        Defined at ../src/cell.f90 lines 90-323
        
        Parameters
        ----------
        kpoints : float array
        num_kpts : int
        kpoint_grid_dim : int array
        kpoint_offset : float array
        
        =========================================================================
        """
        _OptaPyDOS.f90wrap_cell_find_mp_grid(kpoints=kpoints, num_kpts=num_kpts, \
            kpoint_grid_dim=kpoint_grid_dim, kpoint_offset=kpoint_offset)
    
    @staticmethod
    def cell_get_symmetry():
        """
        cell_get_symmetry()
        
        
        Defined at ../src/cell.f90 lines 330-408
        
        
        =========================================================================
         Read in the cell symmetries
        -------------------------------------------------------------------------
         Arguments: kpoints - an array of kpoints
                    num_kpts - size of the kpoint array
        -------------------------------------------------------------------------
         Returns: kpint_grid_dim - the number of kpoints in each dimension
        -------------------------------------------------------------------------
         Parent module variables used: None
        -------------------------------------------------------------------------
         Modules used:  None
        -------------------------------------------------------------------------
         Key Internal Variables:
         Described below
        -------------------------------------------------------------------------
         Necessary conditions: None
        --------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         JRY, April 2011
        =========================================================================
        """
        _OptaPyDOS.f90wrap_cell_get_symmetry()
    
    @staticmethod
    def cell_read_cell():
        """
        cell_read_cell()
        
        
        Defined at ../src/cell.f90 lines 413-712
        
        
        =========================================================================
        """
        _OptaPyDOS.f90wrap_cell_read_cell()
    
    @staticmethod
    def cell_get_atoms():
        """
        cell_get_atoms()
        
        
        Defined at ../src/cell.f90 lines 717-953
        
        
        =========================================================================
        """
        _OptaPyDOS.f90wrap_cell_get_atoms()
    
    @staticmethod
    def cell_calc_lattice():
        """
        cell_calc_lattice()
        
        
        Defined at ../src/cell.f90 lines 956-1013
        
        
        =========================================================================
         Begin with a real lattice. Convert from bohr. Calculate a reciprocal
         lattice and the volume of the cell
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables used: real_lattice, recip_lattice, cell_volume
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by Andrew Morris from the LinDOS program 11/10/2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_cell_calc_lattice()
    
    @staticmethod
    def cell_report_parameters():
        """
        cell_report_parameters()
        
        
        Defined at ../src/cell.f90 lines 1016-1063
        
        
        =========================================================================
         Begin with a real lattice. Convert from bohr. Calculate a reciprocal
         lattice and the volume of the cell
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables used: real_lattice, recip_lattice, cell_volume
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by J R Yates, modified A J Morris                     Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_cell_report_parameters()
    
    @staticmethod
    def cell_dist():
        """
        cell_dist()
        
        
        Defined at ../src/cell.f90 lines 1065-1124
        
        
        -------------------------------------------------------------------------
        """
        _OptaPyDOS.f90wrap_cell_dist()
    
    @property
    def real_lattice(self):
        """
        Element real_lattice ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__real_lattice(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            real_lattice = self._arrays[array_handle]
        else:
            real_lattice = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__real_lattice)
            self._arrays[array_handle] = real_lattice
        return real_lattice
    
    @real_lattice.setter
    def real_lattice(self, real_lattice):
        self.real_lattice[...] = real_lattice
    
    @property
    def recip_lattice(self):
        """
        Element recip_lattice ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__recip_lattice(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            recip_lattice = self._arrays[array_handle]
        else:
            recip_lattice = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__recip_lattice)
            self._arrays[array_handle] = recip_lattice
        return recip_lattice
    
    @recip_lattice.setter
    def recip_lattice(self, recip_lattice):
        self.recip_lattice[...] = recip_lattice
    
    @property
    def cell_volume(self):
        """
        Element cell_volume ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 45
        
        """
        return _OptaPyDOS.f90wrap_od_cell__get__cell_volume()
    
    @cell_volume.setter
    def cell_volume(self, cell_volume):
        _OptaPyDOS.f90wrap_od_cell__set__cell_volume(cell_volume)
    
    @property
    def kpoint_r(self):
        """
        Element kpoint_r ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__kpoint_r(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kpoint_r = self._arrays[array_handle]
        else:
            kpoint_r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__kpoint_r)
            self._arrays[array_handle] = kpoint_r
        return kpoint_r
    
    @kpoint_r.setter
    def kpoint_r(self, kpoint_r):
        self.kpoint_r[...] = kpoint_r
    
    @property
    def kpoint_r_cart(self):
        """
        Element kpoint_r_cart ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__kpoint_r_cart(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kpoint_r_cart = self._arrays[array_handle]
        else:
            kpoint_r_cart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__kpoint_r_cart)
            self._arrays[array_handle] = kpoint_r_cart
        return kpoint_r_cart
    
    @kpoint_r_cart.setter
    def kpoint_r_cart(self, kpoint_r_cart):
        self.kpoint_r_cart[...] = kpoint_r_cart
    
    @property
    def kpoint_weight(self):
        """
        Element kpoint_weight ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__kpoint_weight(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kpoint_weight = self._arrays[array_handle]
        else:
            kpoint_weight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__kpoint_weight)
            self._arrays[array_handle] = kpoint_weight
        return kpoint_weight
    
    @kpoint_weight.setter
    def kpoint_weight(self, kpoint_weight):
        self.kpoint_weight[...] = kpoint_weight
    
    @property
    def num_kpoints_on_node(self):
        """
        Element num_kpoints_on_node ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__num_kpoints_on_node(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            num_kpoints_on_node = self._arrays[array_handle]
        else:
            num_kpoints_on_node = \
                f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__num_kpoints_on_node)
            self._arrays[array_handle] = num_kpoints_on_node
        return num_kpoints_on_node
    
    @num_kpoints_on_node.setter
    def num_kpoints_on_node(self, num_kpoints_on_node):
        self.num_kpoints_on_node[...] = num_kpoints_on_node
    
    @property
    def nkpoints(self):
        """
        Element nkpoints ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 57
        
        """
        return _OptaPyDOS.f90wrap_od_cell__get__nkpoints()
    
    @nkpoints.setter
    def nkpoints(self, nkpoints):
        _OptaPyDOS.f90wrap_od_cell__set__nkpoints(nkpoints)
    
    @property
    def kpoint_grid_dim(self):
        """
        Element kpoint_grid_dim ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 58
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__kpoint_grid_dim(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kpoint_grid_dim = self._arrays[array_handle]
        else:
            kpoint_grid_dim = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__kpoint_grid_dim)
            self._arrays[array_handle] = kpoint_grid_dim
        return kpoint_grid_dim
    
    @kpoint_grid_dim.setter
    def kpoint_grid_dim(self, kpoint_grid_dim):
        self.kpoint_grid_dim[...] = kpoint_grid_dim
    
    @property
    def num_crystal_symmetry_operations(self):
        """
        Element num_crystal_symmetry_operations ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 62
        
        """
        return _OptaPyDOS.f90wrap_od_cell__get__num_crystal_symmetry_operations()
    
    @num_crystal_symmetry_operations.setter
    def num_crystal_symmetry_operations(self, num_crystal_symmetry_operations):
        \
            _OptaPyDOS.f90wrap_od_cell__set__num_crystal_symmetry_operations(num_crystal_symmetry_operations)
    
    @property
    def crystal_symmetry_disps(self):
        """
        Element crystal_symmetry_disps ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__crystal_symmetry_disps(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            crystal_symmetry_disps = self._arrays[array_handle]
        else:
            crystal_symmetry_disps = \
                f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__crystal_symmetry_disps)
            self._arrays[array_handle] = crystal_symmetry_disps
        return crystal_symmetry_disps
    
    @crystal_symmetry_disps.setter
    def crystal_symmetry_disps(self, crystal_symmetry_disps):
        self.crystal_symmetry_disps[...] = crystal_symmetry_disps
    
    @property
    def crystal_symmetry_operations(self):
        """
        Element crystal_symmetry_operations ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__crystal_symmetry_operations(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            crystal_symmetry_operations = self._arrays[array_handle]
        else:
            crystal_symmetry_operations = \
                f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__crystal_symmetry_operations)
            self._arrays[array_handle] = crystal_symmetry_operations
        return crystal_symmetry_operations
    
    @crystal_symmetry_operations.setter
    def crystal_symmetry_operations(self, crystal_symmetry_operations):
        self.crystal_symmetry_operations[...] = crystal_symmetry_operations
    
    @property
    def atoms_pos_frac(self):
        """
        Element atoms_pos_frac ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 67
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__atoms_pos_frac(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            atoms_pos_frac = self._arrays[array_handle]
        else:
            atoms_pos_frac = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__atoms_pos_frac)
            self._arrays[array_handle] = atoms_pos_frac
        return atoms_pos_frac
    
    @atoms_pos_frac.setter
    def atoms_pos_frac(self, atoms_pos_frac):
        self.atoms_pos_frac[...] = atoms_pos_frac
    
    @property
    def atoms_pos_cart(self):
        """
        Element atoms_pos_cart ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/cell.f90 line 68
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__atoms_pos_cart(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            atoms_pos_cart = self._arrays[array_handle]
        else:
            atoms_pos_cart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__atoms_pos_cart)
            self._arrays[array_handle] = atoms_pos_cart
        return atoms_pos_cart
    
    @atoms_pos_cart.setter
    def atoms_pos_cart(self, atoms_pos_cart):
        self.atoms_pos_cart[...] = atoms_pos_cart
    
    @property
    def atoms_species_num(self):
        """
        Element atoms_species_num ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 69
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__atoms_species_num(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            atoms_species_num = self._arrays[array_handle]
        else:
            atoms_species_num = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__atoms_species_num)
            self._arrays[array_handle] = atoms_species_num
        return atoms_species_num
    
    @atoms_species_num.setter
    def atoms_species_num(self, atoms_species_num):
        self.atoms_species_num[...] = atoms_species_num
    
    @property
    def atoms_label(self):
        """
        Element atoms_label ftype=character(len=maxlen) pytype=str
        
        
        Defined at ../src/cell.f90 line 70
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__atoms_label(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            atoms_label = self._arrays[array_handle]
        else:
            atoms_label = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__atoms_label)
            self._arrays[array_handle] = atoms_label
        return atoms_label
    
    @atoms_label.setter
    def atoms_label(self, atoms_label):
        self.atoms_label[...] = atoms_label
    
    @property
    def atoms_symbol(self):
        """
        Element atoms_symbol ftype=character(len=2) pytype=str
        
        
        Defined at ../src/cell.f90 line 71
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_cell__array__atoms_symbol(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            atoms_symbol = self._arrays[array_handle]
        else:
            atoms_symbol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_cell__array__atoms_symbol)
            self._arrays[array_handle] = atoms_symbol
        return atoms_symbol
    
    @atoms_symbol.setter
    def atoms_symbol(self, atoms_symbol):
        self.atoms_symbol[...] = atoms_symbol
    
    @property
    def num_atoms(self):
        """
        Element num_atoms ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 72
        
        """
        return _OptaPyDOS.f90wrap_od_cell__get__num_atoms()
    
    @num_atoms.setter
    def num_atoms(self, num_atoms):
        _OptaPyDOS.f90wrap_od_cell__set__num_atoms(num_atoms)
    
    @property
    def num_species(self):
        """
        Element num_species ftype=integer pytype=int
        
        
        Defined at ../src/cell.f90 line 73
        
        """
        return _OptaPyDOS.f90wrap_od_cell__get__num_species()
    
    @num_species.setter
    def num_species(self, num_species):
        _OptaPyDOS.f90wrap_od_cell__set__num_species(num_species)
    
    def __str__(self):
        ret = ['<od_cell>{\n']
        ret.append('    real_lattice : ')
        ret.append(repr(self.real_lattice))
        ret.append(',\n    recip_lattice : ')
        ret.append(repr(self.recip_lattice))
        ret.append(',\n    cell_volume : ')
        ret.append(repr(self.cell_volume))
        ret.append(',\n    kpoint_r : ')
        ret.append(repr(self.kpoint_r))
        ret.append(',\n    kpoint_r_cart : ')
        ret.append(repr(self.kpoint_r_cart))
        ret.append(',\n    kpoint_weight : ')
        ret.append(repr(self.kpoint_weight))
        ret.append(',\n    num_kpoints_on_node : ')
        ret.append(repr(self.num_kpoints_on_node))
        ret.append(',\n    nkpoints : ')
        ret.append(repr(self.nkpoints))
        ret.append(',\n    kpoint_grid_dim : ')
        ret.append(repr(self.kpoint_grid_dim))
        ret.append(',\n    num_crystal_symmetry_operations : ')
        ret.append(repr(self.num_crystal_symmetry_operations))
        ret.append(',\n    crystal_symmetry_disps : ')
        ret.append(repr(self.crystal_symmetry_disps))
        ret.append(',\n    crystal_symmetry_operations : ')
        ret.append(repr(self.crystal_symmetry_operations))
        ret.append(',\n    atoms_pos_frac : ')
        ret.append(repr(self.atoms_pos_frac))
        ret.append(',\n    atoms_pos_cart : ')
        ret.append(repr(self.atoms_pos_cart))
        ret.append(',\n    atoms_species_num : ')
        ret.append(repr(self.atoms_species_num))
        ret.append(',\n    atoms_label : ')
        ret.append(repr(self.atoms_label))
        ret.append(',\n    atoms_symbol : ')
        ret.append(repr(self.atoms_symbol))
        ret.append(',\n    num_atoms : ')
        ret.append(repr(self.num_atoms))
        ret.append(',\n    num_species : ')
        ret.append(repr(self.num_species))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_cell = Od_Cell()

class Od_Comms(f90wrap.runtime.FortranModule):
    """
    Module od_comms
    
    
    Defined at ../src/comms.F90 lines 29-665
    
    """
    @staticmethod
    def comms_setup():
        """
        comms_setup()
        
        
        Defined at ../src/comms.F90 lines 87-107
        
        
        """
        _OptaPyDOS.f90wrap_comms_setup()
    
    @staticmethod
    def comms_end():
        """
        comms_end()
        
        
        Defined at ../src/comms.F90 lines 109-122
        
        
        """
        _OptaPyDOS.f90wrap_comms_end()
    
    @staticmethod
    def _comms_bcast_int(array, size_bn):
        """
        _comms_bcast_int(array, size_bn)
        
        
        Defined at ../src/comms.F90 lines 124-145
        
        Parameters
        ----------
        array : int
        size_bn : int
        
        """
        _OptaPyDOS.f90wrap_comms_bcast_int(array=array, size_bn=size_bn)
    
    @staticmethod
    def _comms_bcast_logical(array, size_bn):
        """
        _comms_bcast_logical(array, size_bn)
        
        
        Defined at ../src/comms.F90 lines 170-191
        
        Parameters
        ----------
        array : bool
        size_bn : int
        
        """
        _OptaPyDOS.f90wrap_comms_bcast_logical(array=array, size_bn=size_bn)
    
    @staticmethod
    def _comms_bcast_real(array, size_bn):
        """
        _comms_bcast_real(array, size_bn)
        
        
        Defined at ../src/comms.F90 lines 147-168
        
        Parameters
        ----------
        array : float
        size_bn : int
        
        """
        _OptaPyDOS.f90wrap_comms_bcast_real(array=array, size_bn=size_bn)
    
    @staticmethod
    def _comms_bcast_cmplx(array, size_bn):
        """
        _comms_bcast_cmplx(array, size_bn)
        
        
        Defined at ../src/comms.F90 lines 216-238
        
        Parameters
        ----------
        array : complex
        size_bn : int
        
        """
        _OptaPyDOS.f90wrap_comms_bcast_cmplx(array=array, size_bn=size_bn)
    
    @staticmethod
    def _comms_bcast_char(array, size_bn):
        """
        _comms_bcast_char(array, size_bn)
        
        
        Defined at ../src/comms.F90 lines 193-214
        
        Parameters
        ----------
        array : str
        size_bn : int
        
        """
        _OptaPyDOS.f90wrap_comms_bcast_char(array=array, size_bn=size_bn)
    
    @staticmethod
    def comms_bcast(*args, **kwargs):
        """
        comms_bcast(*args, **kwargs)
        
        
        Defined at ../src/comms.F90 lines 54-60
        
        Overloaded interface containing the following procedures:
          _comms_bcast_int
          _comms_bcast_logical
          _comms_bcast_real
          _comms_bcast_cmplx
          _comms_bcast_char
        
        """
        for proc in [_comms_bcast_int, _comms_bcast_logical, _comms_bcast_real, \
            _comms_bcast_cmplx, _comms_bcast_char]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _comms_send_int(array, size_bn, to):
        """
        _comms_send_int(array, size_bn, to)
        
        
        Defined at ../src/comms.F90 lines 268-292
        
        Parameters
        ----------
        array : int
        size_bn : int
        to : int
        
        """
        _OptaPyDOS.f90wrap_comms_send_int(array=array, size_bn=size_bn, to=to)
    
    @staticmethod
    def _comms_send_logical(array, size_bn, to):
        """
        _comms_send_logical(array, size_bn, to)
        
        
        Defined at ../src/comms.F90 lines 242-266
        
        Parameters
        ----------
        array : bool
        size_bn : int
        to : int
        
        """
        _OptaPyDOS.f90wrap_comms_send_logical(array=array, size_bn=size_bn, to=to)
    
    @staticmethod
    def _comms_send_real(array, size_bn, to):
        """
        _comms_send_real(array, size_bn, to)
        
        
        Defined at ../src/comms.F90 lines 320-344
        
        Parameters
        ----------
        array : float
        size_bn : int
        to : int
        
        """
        _OptaPyDOS.f90wrap_comms_send_real(array=array, size_bn=size_bn, to=to)
    
    @staticmethod
    def _comms_send_cmplx(array, size_bn, to):
        """
        _comms_send_cmplx(array, size_bn, to)
        
        
        Defined at ../src/comms.F90 lines 346-370
        
        Parameters
        ----------
        array : complex
        size_bn : int
        to : int
        
        """
        _OptaPyDOS.f90wrap_comms_send_cmplx(array=array, size_bn=size_bn, to=to)
    
    @staticmethod
    def _comms_send_char(array, size_bn, to):
        """
        _comms_send_char(array, size_bn, to)
        
        
        Defined at ../src/comms.F90 lines 294-318
        
        Parameters
        ----------
        array : str
        size_bn : int
        to : int
        
        """
        _OptaPyDOS.f90wrap_comms_send_char(array=array, size_bn=size_bn, to=to)
    
    @staticmethod
    def comms_send(*args, **kwargs):
        """
        comms_send(*args, **kwargs)
        
        
        Defined at ../src/comms.F90 lines 62-68
        
        Overloaded interface containing the following procedures:
          _comms_send_int
          _comms_send_logical
          _comms_send_real
          _comms_send_cmplx
          _comms_send_char
        
        """
        for proc in [_comms_send_int, _comms_send_logical, _comms_send_real, \
            _comms_send_cmplx, _comms_send_char]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _comms_recv_int(array, size_bn, from_):
        """
        _comms_recv_int(array, size_bn, from_)
        
        
        Defined at ../src/comms.F90 lines 401-426
        
        Parameters
        ----------
        array : int
        size_bn : int
        from_ : int
        
        """
        _OptaPyDOS.f90wrap_comms_recv_int(array=array, size_bn=size_bn, from_=from_)
    
    @staticmethod
    def _comms_recv_logical(array, size_bn, from_):
        """
        _comms_recv_logical(array, size_bn, from_)
        
        
        Defined at ../src/comms.F90 lines 374-399
        
        Parameters
        ----------
        array : bool
        size_bn : int
        from_ : int
        
        """
        _OptaPyDOS.f90wrap_comms_recv_logical(array=array, size_bn=size_bn, from_=from_)
    
    @staticmethod
    def _comms_recv_real(array, size_bn, from_):
        """
        _comms_recv_real(array, size_bn, from_)
        
        
        Defined at ../src/comms.F90 lines 455-480
        
        Parameters
        ----------
        array : float
        size_bn : int
        from_ : int
        
        """
        _OptaPyDOS.f90wrap_comms_recv_real(array=array, size_bn=size_bn, from_=from_)
    
    @staticmethod
    def _comms_recv_cmplx(array, size_bn, from_):
        """
        _comms_recv_cmplx(array, size_bn, from_)
        
        
        Defined at ../src/comms.F90 lines 482-509
        
        Parameters
        ----------
        array : complex
        size_bn : int
        from_ : int
        
        """
        _OptaPyDOS.f90wrap_comms_recv_cmplx(array=array, size_bn=size_bn, from_=from_)
    
    @staticmethod
    def _comms_recv_char(array, size_bn, from_):
        """
        _comms_recv_char(array, size_bn, from_)
        
        
        Defined at ../src/comms.F90 lines 428-453
        
        Parameters
        ----------
        array : str
        size_bn : int
        from_ : int
        
        """
        _OptaPyDOS.f90wrap_comms_recv_char(array=array, size_bn=size_bn, from_=from_)
    
    @staticmethod
    def comms_recv(*args, **kwargs):
        """
        comms_recv(*args, **kwargs)
        
        
        Defined at ../src/comms.F90 lines 70-76
        
        Overloaded interface containing the following procedures:
          _comms_recv_int
          _comms_recv_logical
          _comms_recv_real
          _comms_recv_cmplx
          _comms_recv_char
        
        """
        for proc in [_comms_recv_int, _comms_recv_logical, _comms_recv_real, \
            _comms_recv_cmplx, _comms_recv_char]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _comms_reduce_int(array, size_bn, op):
        """
        _comms_reduce_int(array, size_bn, op)
        
        
        Defined at ../src/comms.F90 lines 528-571
        
        Parameters
        ----------
        array : int
        size_bn : int
        op : str
        
        """
        _OptaPyDOS.f90wrap_comms_reduce_int(array=array, size_bn=size_bn, op=op)
    
    @staticmethod
    def _comms_reduce_real(array, size_bn, op):
        """
        _comms_reduce_real(array, size_bn, op)
        
        
        Defined at ../src/comms.F90 lines 573-620
        
        Parameters
        ----------
        array : float
        size_bn : int
        op : str
        
        """
        _OptaPyDOS.f90wrap_comms_reduce_real(array=array, size_bn=size_bn, op=op)
    
    @staticmethod
    def _comms_reduce_cmplx(array, size_bn, op):
        """
        _comms_reduce_cmplx(array, size_bn, op)
        
        
        Defined at ../src/comms.F90 lines 622-663
        
        Parameters
        ----------
        array : complex
        size_bn : int
        op : str
        
        """
        _OptaPyDOS.f90wrap_comms_reduce_cmplx(array=array, size_bn=size_bn, op=op)
    
    @staticmethod
    def comms_reduce(*args, **kwargs):
        """
        comms_reduce(*args, **kwargs)
        
        
        Defined at ../src/comms.F90 lines 78-83
        
        Overloaded interface containing the following procedures:
          _comms_reduce_int
          _comms_reduce_real
          _comms_reduce_cmplx
        
        """
        for proc in [_comms_reduce_int, _comms_reduce_real, _comms_reduce_cmplx]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def on_root(self):
        """
        Element on_root ftype=logical pytype=bool
        
        
        Defined at ../src/comms.F90 line 41
        
        """
        return _OptaPyDOS.f90wrap_od_comms__get__on_root()
    
    @on_root.setter
    def on_root(self, on_root):
        _OptaPyDOS.f90wrap_od_comms__set__on_root(on_root)
    
    @property
    def num_nodes(self):
        """
        Element num_nodes ftype=integer pytype=int
        
        
        Defined at ../src/comms.F90 line 42
        
        """
        return _OptaPyDOS.f90wrap_od_comms__get__num_nodes()
    
    @num_nodes.setter
    def num_nodes(self, num_nodes):
        _OptaPyDOS.f90wrap_od_comms__set__num_nodes(num_nodes)
    
    @property
    def my_node_id(self):
        """
        Element my_node_id ftype=integer pytype=int
        
        
        Defined at ../src/comms.F90 line 42
        
        """
        return _OptaPyDOS.f90wrap_od_comms__get__my_node_id()
    
    @my_node_id.setter
    def my_node_id(self, my_node_id):
        _OptaPyDOS.f90wrap_od_comms__set__my_node_id(my_node_id)
    
    @property
    def root_id(self):
        """
        Element root_id ftype=integer pytype=int
        
        
        Defined at ../src/comms.F90 line 43
        
        """
        return _OptaPyDOS.f90wrap_od_comms__get__root_id()
    
    def __str__(self):
        ret = ['<od_comms>{\n']
        ret.append('    on_root : ')
        ret.append(repr(self.on_root))
        ret.append(',\n    num_nodes : ')
        ret.append(repr(self.num_nodes))
        ret.append(',\n    my_node_id : ')
        ret.append(repr(self.my_node_id))
        ret.append(',\n    root_id : ')
        ret.append(repr(self.root_id))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_comms = Od_Comms()

class Od_Constants(f90wrap.runtime.FortranModule):
    """
    Module od_constants
    
    
    Defined at ../src/constants.f90 lines 32-65
    
    """
    @property
    def optados_version(self):
        """
        Element optados_version ftype=character(len=6) pytype=str
        
        
        Defined at ../src/constants.f90 line 37
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__optados_version()
    
    @property
    def copyright(self):
        """
        Element copyright ftype=character(len=14) pytype=str
        
        
        Defined at ../src/constants.f90 line 38
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__copyright()
    
    @property
    def dp(self):
        """
        Element dp ftype=integer pytype=int
        
        
        Defined at ../src/constants.f90 line 40
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__dp()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/constants.f90 line 42
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__pi()
    
    @property
    def h2ev(self):
        """
        Element h2ev ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/constants.f90 line 43
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__h2ev()
    
    @property
    def bohr2ang(self):
        """
        Element bohr2ang ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/constants.f90 line 44
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__bohr2ang()
    
    @property
    def inv_sqrt_two_pi(self):
        """
        Element inv_sqrt_two_pi ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/constants.f90 line 45
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__inv_sqrt_two_pi()
    
    @property
    def twopi(self):
        """
        Element twopi ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/constants.f90 line 46
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__twopi()
    
    @property
    def sqrt_two(self):
        """
        Element sqrt_two ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/constants.f90 line 47
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__sqrt_two()
    
    @property
    def cmplx_0(self):
        """
        Element cmplx_0 ftype=complex(dp) pytype=complex
        
        
        Defined at ../src/constants.f90 line 48
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__cmplx_0()
    
    @property
    def cmplx_i(self):
        """
        Element cmplx_i ftype=complex(dp) pytype=complex
        
        
        Defined at ../src/constants.f90 line 49
        
        """
        return _OptaPyDOS.f90wrap_od_constants__get__cmplx_i()
    
    def __str__(self):
        ret = ['<od_constants>{\n']
        ret.append('    optados_version : ')
        ret.append(repr(self.optados_version))
        ret.append(',\n    copyright : ')
        ret.append(repr(self.copyright))
        ret.append(',\n    dp : ')
        ret.append(repr(self.dp))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    h2ev : ')
        ret.append(repr(self.h2ev))
        ret.append(',\n    bohr2ang : ')
        ret.append(repr(self.bohr2ang))
        ret.append(',\n    inv_sqrt_two_pi : ')
        ret.append(repr(self.inv_sqrt_two_pi))
        ret.append(',\n    twopi : ')
        ret.append(repr(self.twopi))
        ret.append(',\n    sqrt_two : ')
        ret.append(repr(self.sqrt_two))
        ret.append(',\n    cmplx_0 : ')
        ret.append(repr(self.cmplx_0))
        ret.append(',\n    cmplx_i : ')
        ret.append(repr(self.cmplx_i))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_constants = Od_Constants()

class Od_Core(f90wrap.runtime.FortranModule):
    """
    Module od_core
    
    
    Defined at ../src/core.f90 lines 26-747
    
    """
    @staticmethod
    def core_calculate():
        """
        core_calculate()
        
        
        Defined at ../src/core.f90 lines 41-82
        
        
        """
        _OptaPyDOS.f90wrap_core_calculate()
    
    @property
    def matrix_weights(self):
        """
        Element matrix_weights ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/core.f90 line 33
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_core__array__matrix_weights(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            matrix_weights = self._arrays[array_handle]
        else:
            matrix_weights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_core__array__matrix_weights)
            self._arrays[array_handle] = matrix_weights
        return matrix_weights
    
    @matrix_weights.setter
    def matrix_weights(self, matrix_weights):
        self.matrix_weights[...] = matrix_weights
    
    @property
    def weighted_dos(self):
        """
        Element weighted_dos ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/core.f90 line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_core__array__weighted_dos(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            weighted_dos = self._arrays[array_handle]
        else:
            weighted_dos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_core__array__weighted_dos)
            self._arrays[array_handle] = weighted_dos
        return weighted_dos
    
    @weighted_dos.setter
    def weighted_dos(self, weighted_dos):
        self.weighted_dos[...] = weighted_dos
    
    @property
    def weighted_dos_broadened(self):
        """
        Element weighted_dos_broadened ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/core.f90 line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_core__array__weighted_dos_broadened(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            weighted_dos_broadened = self._arrays[array_handle]
        else:
            weighted_dos_broadened = \
                f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_core__array__weighted_dos_broadened)
            self._arrays[array_handle] = weighted_dos_broadened
        return weighted_dos_broadened
    
    @weighted_dos_broadened.setter
    def weighted_dos_broadened(self, weighted_dos_broadened):
        self.weighted_dos_broadened[...] = weighted_dos_broadened
    
    def __str__(self):
        ret = ['<od_core>{\n']
        ret.append('    matrix_weights : ')
        ret.append(repr(self.matrix_weights))
        ret.append(',\n    weighted_dos : ')
        ret.append(repr(self.weighted_dos))
        ret.append(',\n    weighted_dos_broadened : ')
        ret.append(repr(self.weighted_dos_broadened))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_core = Od_Core()

class Od_Dos(f90wrap.runtime.FortranModule):
    """
    Module od_dos
    
    
    Defined at ../src/dos.f90 lines 24-323
    
    """
    @staticmethod
    def dos_calculate():
        """
        dos_calculate()
        
        
        Defined at ../src/dos.f90 lines 35-153
        
        
        ===============================================================================
         Main routine in dos module, drives the calculation of density of states for
         both task : dos and also if it is required elsewhere.
        -------------------------------------------------------------------------------
         Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
                    weighted_dos (out)(opt) : Output DOS weigthed by matrix_weights
        -------------------------------------------------------------------------------
         Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
         intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
         efermi_linear, delta_bins, calc_weighted_dos
        -------------------------------------------------------------------------------
         Modules Used: see below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: One of linear, adaptive or fixed must be .true.
        -------------------------------------------------------------------------------
         Known Worries: (1) If more than one of linear, adaptive or fixed are set it
         uses the most complicated method.
         (2) It should be possible to pass optioinal arguments to sub programs as
         optional argumnets without checking whether they are there or not. g95 will
         allow this behaviour. gfotran will not.
        -------------------------------------------------------------------------------
         Written by : A J Morris December 2010
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_dos_calculate()
    
    @staticmethod
    def write_dos(e, dos, intdos, dos_name):
        """
        write_dos(e, dos, intdos, dos_name)
        
        
        Defined at ../src/dos.f90 lines 156-257
        
        Parameters
        ----------
        e : float array
        dos : float array
        intdos : float array
        dos_name : str
        
        ===============================================================================
         This routine receives an energy scale, a density of states and a file name
         and writes out the DOS to disk
        -------------------------------------------------------------------------------
         Arguments: E       (in) : The energy scale
                    dos     (in) : The density of states
                    intdos  (in) : The integrated DOS
                    dos_name(in) : Name of the output file
        -------------------------------------------------------------------------------
         Parent Module Varables Used: None
        -------------------------------------------------------------------------------
         Modules Used: See below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: None
        -------------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------------
         Written by : A J Morris December 2010
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_write_dos(e=e, dos=dos, intdos=intdos, dos_name=dos_name)
    
    @staticmethod
    def write_dos_xmgrace(dos_name, e, dos):
        """
        write_dos_xmgrace(dos_name, e, dos)
        
        
        Defined at ../src/dos.f90 lines 260-323
        
        Parameters
        ----------
        dos_name : str
        e : float array
        dos : float array
        
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_write_dos_xmgrace(dos_name=dos_name, e=e, dos=dos)
    
    _dt_array_initialisers = []
    

od_dos = Od_Dos()

class Od_Dos_Utils(f90wrap.runtime.FortranModule):
    """
    Module od_dos_utils
    
    
    Defined at ../src/dos_utils.f90 lines 34-1863
    
    """
    @staticmethod
    def dos_utils_calculate():
        """
        dos_utils_calculate()
        
        
        Defined at ../src/dos_utils.f90 lines 92-297
        
        
        ===============================================================================
         Main routine in dos module, drives the calculation of density of states for
         both task : dos and also if it is required elsewhere.
        -------------------------------------------------------------------------------
         Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
                    weighted_dos (out)(opt) : Output DOS weigthed by matrix_weights
        -------------------------------------------------------------------------------
         Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
         intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
         efermi_linear, delta_bins, calc_weighted_dos
        -------------------------------------------------------------------------------
         Modules Used: see below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: One of linear, adaptive or fixed must be .true.
        -------------------------------------------------------------------------------
         Known Worries: (1) If more than one of linear, adaptive or fixed are set it
         uses the most complicated method.
         (2) It should be possible to pass optioinal arguments to sub programs as
         optional argumnets without checking whether they are there or not. g95 will
         allow this behaviour. gfotran will not.
        -------------------------------------------------------------------------------
         Written by : A J Morris December 2010
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_dos_utils_calculate()
    
    @staticmethod
    def dos_utils_set_efermi():
        """
        dos_utils_set_efermi()
        
        
        Defined at ../src/dos_utils.f90 lines 300-393
        
        
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_dos_utils_set_efermi()
    
    @staticmethod
    def dos_utils_compute_dos_at_efermi():
        """
        dos_utils_compute_dos_at_efermi()
        
        
        Defined at ../src/dos_utils.f90 lines 396-460
        
        
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_dos_utils_compute_dos_at_efermi()
    
    @staticmethod
    def dos_utils_compute_bandgap():
        """
        dos_utils_compute_bandgap()
        
        
        Defined at ../src/dos_utils.f90 lines 463-749
        
        
        ===============================================================================
         Modified from LINDOS -- AJM 3rd June 2011
         Rewritten 31/1/12 AJM
        """
        _OptaPyDOS.f90wrap_dos_utils_compute_bandgap()
    
    @staticmethod
    def dos_utils_compute_band_energies():
        """
        dos_utils_compute_band_energies()
        
        
        Defined at ../src/dos_utils.f90 lines 856-929
        
        
        ===============================================================================
         High-level subroutine to compute band energies of the DOS calculated.
         Calculates using the band_energies directly and compares with the
         function calc_band_energies which does the low level computation on the DOS.
        -------------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------------
         Parent Module Varables Used: E,dos_fixed,dos_adaptive,dos_linear
        -------------------------------------------------------------------------------
         Modules Used: See below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: E must be set
        -------------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------------
         Written by : A J Morris December 2010
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_dos_utils_compute_band_energies()
    
    @staticmethod
    def dos_utils_deallocate():
        """
        dos_utils_deallocate()
        
        
        Defined at ../src/dos_utils.f90 lines 1113-1146
        
        
        """
        _OptaPyDOS.f90wrap_dos_utils_deallocate()
    
    @staticmethod
    def doslin_sub_cell_corners(grad, step, energy, eigenv):
        """
        doslin_sub_cell_corners(grad, step, energy, eigenv)
        
        
        Defined at ../src/dos_utils.f90 lines 1333-1402
        
        Parameters
        ----------
        grad : float array
        step : float array
        energy : float
        eigenv : float array
        
        ===============================================================================
         A helper subroutine for calculated_dos, which is used for the linear
         broadening method. This routine extrapolates the energies at the corner of the
         sub cells by using the gradient at the centre of the cell
        -------------------------------------------------------------------------------
         Arguments: grad (in) : The Gradient of the band at the centre of a sub-cell
                    step   (in) : The directions to the edge of the sub_cell
                    energy (in) : The Band energy at the centre of the sub cell
                    EigenV (out): The Energies at the corners of the sub-cell
        -------------------------------------------------------------------------------
         Parent Module Varables Used: None
        -------------------------------------------------------------------------------
         Modules Used: See below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: None
        -------------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------------
         Written by : A J Morris December 2010 Heavliy modified from LinDOS
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_doslin_sub_cell_corners(grad=grad, step=step, energy=energy, \
            eigenv=eigenv)
    
    @staticmethod
    def doslin(e0, e1, e2, e3, e4, e):
        """
        doslin, int_bn = doslin(e0, e1, e2, e3, e4, e)
        
        
        Defined at ../src/dos_utils.f90 lines 1405-1542
        
        Parameters
        ----------
        e0 : float
        e1 : float
        e2 : float
        e3 : float
        e4 : float
        e : float
        
        Returns
        -------
        doslin : float
        int_bn : float
        
        ===============================================================================
         Return the DoS contribution for a linear band portion and a cubic cell
        -------------------------------------------------------------------------------
         Arguments: e0 (in) : Energy at centre of sub cell
           e1,e2,e3,e4 (in) : Energies of the four lowest corners of the sub cell
                    e  (in) : Energy at which DOS is evaluated
                  int  (out): Integrated DOS contribution for energy, E)
         (The function itself returns the DOS couribution for energy, E)
        -------------------------------------------------------------------------------
         Parent Module Varables Used: None
        -------------------------------------------------------------------------------
         Modules Used: See below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: None
        -------------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------------
         Written by : C J Pickard. From LinDOS. Extra Comments A J Morris Sept 2010
        ===============================================================================
        """
        doslin, int_bn = _OptaPyDOS.f90wrap_doslin(e0=e0, e1=e1, e2=e2, e3=e3, e4=e4, \
            e=e)
        return doslin, int_bn
    
    @staticmethod
    def dos_utils_calculate_at_e(energy, dos_at_e, weighted_dos_at_e=None):
        """
        dos_utils_calculate_at_e(energy, dos_at_e[, weighted_dos_at_e])
        
        
        Defined at ../src/dos_utils.f90 lines 1545-1670
        
        Parameters
        ----------
        energy : float
        dos_at_e : float array
        weighted_dos_at_e : float array
        
        ===============================================================================
         Main routine in dos module, drives the calculation of density of states for
         both task : dos and also if it is required elsewhere.
        -------------------------------------------------------------------------------
         Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
                    weighted_dos (out)(opt) : Output DOS weigthed by matrix_weights
        -------------------------------------------------------------------------------
         Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
         intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
         efermi_linear, delta_bins, calc_weighted_dos
        -------------------------------------------------------------------------------
         Modules Used: see below
        -------------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------------
         Necessary Conditions: One of linear, adaptive or fixed must be .true.
        -------------------------------------------------------------------------------
         Known Worries: (1) If more than one of linear, adaptive or fixed are set it
         uses the most complicated method.
         (2) It should be possible to pass optioinal arguments to sub programs as
         optional argumnets without checking whether they are there or not. g95 will
         allow this behaviour. gfotran will not.
        -------------------------------------------------------------------------------
         Written by : A J Morris December 2010
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_dos_utils_calculate_at_e(energy=energy, dos_at_e=dos_at_e, \
            weighted_dos_at_e=weighted_dos_at_e)
    
    @property
    def dos_adaptive(self):
        """
        Element dos_adaptive ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 47
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__dos_adaptive(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dos_adaptive = self._arrays[array_handle]
        else:
            dos_adaptive = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__dos_adaptive)
            self._arrays[array_handle] = dos_adaptive
        return dos_adaptive
    
    @dos_adaptive.setter
    def dos_adaptive(self, dos_adaptive):
        self.dos_adaptive[...] = dos_adaptive
    
    @property
    def dos_fixed(self):
        """
        Element dos_fixed ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 48
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__dos_fixed(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dos_fixed = self._arrays[array_handle]
        else:
            dos_fixed = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__dos_fixed)
            self._arrays[array_handle] = dos_fixed
        return dos_fixed
    
    @dos_fixed.setter
    def dos_fixed(self, dos_fixed):
        self.dos_fixed[...] = dos_fixed
    
    @property
    def dos_linear(self):
        """
        Element dos_linear ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__dos_linear(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dos_linear = self._arrays[array_handle]
        else:
            dos_linear = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__dos_linear)
            self._arrays[array_handle] = dos_linear
        return dos_linear
    
    @dos_linear.setter
    def dos_linear(self, dos_linear):
        self.dos_linear[...] = dos_linear
    
    @property
    def intdos_adaptive(self):
        """
        Element intdos_adaptive ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__intdos_adaptive(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            intdos_adaptive = self._arrays[array_handle]
        else:
            intdos_adaptive = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__intdos_adaptive)
            self._arrays[array_handle] = intdos_adaptive
        return intdos_adaptive
    
    @intdos_adaptive.setter
    def intdos_adaptive(self, intdos_adaptive):
        self.intdos_adaptive[...] = intdos_adaptive
    
    @property
    def intdos_fixed(self):
        """
        Element intdos_fixed ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__intdos_fixed(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            intdos_fixed = self._arrays[array_handle]
        else:
            intdos_fixed = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__intdos_fixed)
            self._arrays[array_handle] = intdos_fixed
        return intdos_fixed
    
    @intdos_fixed.setter
    def intdos_fixed(self, intdos_fixed):
        self.intdos_fixed[...] = intdos_fixed
    
    @property
    def intdos_linear(self):
        """
        Element intdos_linear ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__intdos_linear(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            intdos_linear = self._arrays[array_handle]
        else:
            intdos_linear = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__intdos_linear)
            self._arrays[array_handle] = intdos_linear
        return intdos_linear
    
    @intdos_linear.setter
    def intdos_linear(self, intdos_linear):
        self.intdos_linear[...] = intdos_linear
    
    @property
    def e(self):
        """
        Element e ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 55
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_dos_utils__array__e(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            e = self._arrays[array_handle]
        else:
            e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_dos_utils__array__e)
            self._arrays[array_handle] = e
        return e
    
    @e.setter
    def e(self, e):
        self.e[...] = e
    
    @property
    def efermi_fixed(self):
        """
        Element efermi_fixed ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 57
        
        """
        return _OptaPyDOS.f90wrap_od_dos_utils__get__efermi_fixed()
    
    @efermi_fixed.setter
    def efermi_fixed(self, efermi_fixed):
        _OptaPyDOS.f90wrap_od_dos_utils__set__efermi_fixed(efermi_fixed)
    
    @property
    def efermi_adaptive(self):
        """
        Element efermi_adaptive ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 58
        
        """
        return _OptaPyDOS.f90wrap_od_dos_utils__get__efermi_adaptive()
    
    @efermi_adaptive.setter
    def efermi_adaptive(self, efermi_adaptive):
        _OptaPyDOS.f90wrap_od_dos_utils__set__efermi_adaptive(efermi_adaptive)
    
    @property
    def efermi_linear(self):
        """
        Element efermi_linear ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/dos_utils.f90 line 59
        
        """
        return _OptaPyDOS.f90wrap_od_dos_utils__get__efermi_linear()
    
    @efermi_linear.setter
    def efermi_linear(self, efermi_linear):
        _OptaPyDOS.f90wrap_od_dos_utils__set__efermi_linear(efermi_linear)
    
    def __str__(self):
        ret = ['<od_dos_utils>{\n']
        ret.append('    dos_adaptive : ')
        ret.append(repr(self.dos_adaptive))
        ret.append(',\n    dos_fixed : ')
        ret.append(repr(self.dos_fixed))
        ret.append(',\n    dos_linear : ')
        ret.append(repr(self.dos_linear))
        ret.append(',\n    intdos_adaptive : ')
        ret.append(repr(self.intdos_adaptive))
        ret.append(',\n    intdos_fixed : ')
        ret.append(repr(self.intdos_fixed))
        ret.append(',\n    intdos_linear : ')
        ret.append(repr(self.intdos_linear))
        ret.append(',\n    e : ')
        ret.append(repr(self.e))
        ret.append(',\n    efermi_fixed : ')
        ret.append(repr(self.efermi_fixed))
        ret.append(',\n    efermi_adaptive : ')
        ret.append(repr(self.efermi_adaptive))
        ret.append(',\n    efermi_linear : ')
        ret.append(repr(self.efermi_linear))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_dos_utils = Od_Dos_Utils()

class Od_Electronic(f90wrap.runtime.FortranModule):
    """
    Module od_electronic
    
    
    Defined at ../src/electronic.f90 lines 31-1078
    
    """
    class Matrix_Weights_Array_Boundaries(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=matrix_weights_array_boundaries)
        
        
        Defined at ../src/electronic.f90 lines 58-64
        
        """
        def __init__(self, handle=None):
            """
            self = Matrix_Weights_Array_Boundaries()
            
            
            Defined at ../src/electronic.f90 lines 58-64
            
            
            Returns
            -------
            this : Matrix_Weights_Array_Boundaries
            	Object to be constructed
            
            
            Automatically generated constructor for matrix_weights_array_boundaries
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _OptaPyDOS.f90wrap_matrix_weights_array_boundaries_initialise()
        
        def __del__(self):
            """
            Destructor for class Matrix_Weights_Array_Boundaries
            
            
            Defined at ../src/electronic.f90 lines 58-64
            
            Parameters
            ----------
            this : Matrix_Weights_Array_Boundaries
            	Object to be destructed
            
            
            Automatically generated destructor for matrix_weights_array_boundaries
            """
            if self._alloc:
                _OptaPyDOS.f90wrap_matrix_weights_array_boundaries_finalise(this=self._handle)
        
        @property
        def norbitals(self):
            """
            Element norbitals ftype=integer  pytype=int
            
            
            Defined at ../src/electronic.f90 line 59
            
            """
            return \
                _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__get__norbitals(self._handle)
        
        @norbitals.setter
        def norbitals(self, norbitals):
            \
                _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__set__norbitals(self._handle, \
                norbitals)
        
        @property
        def nbands(self):
            """
            Element nbands ftype=integer  pytype=int
            
            
            Defined at ../src/electronic.f90 line 60
            
            """
            return \
                _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__get__nbands(self._handle)
        
        @nbands.setter
        def nbands(self, nbands):
            _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__set__nbands(self._handle, \
                nbands)
        
        @property
        def nkpoints(self):
            """
            Element nkpoints ftype=integer  pytype=int
            
            
            Defined at ../src/electronic.f90 line 61
            
            """
            return \
                _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__get__nkpoints(self._handle)
        
        @nkpoints.setter
        def nkpoints(self, nkpoints):
            _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__set__nkpoints(self._handle, \
                nkpoints)
        
        @property
        def nspins(self):
            """
            Element nspins ftype=integer  pytype=int
            
            
            Defined at ../src/electronic.f90 line 62
            
            """
            return \
                _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__get__nspins(self._handle)
        
        @nspins.setter
        def nspins(self, nspins):
            _OptaPyDOS.f90wrap_matrix_weights_array_boundaries__set__nspins(self._handle, \
                nspins)
        
        def __str__(self):
            ret = ['<matrix_weights_array_boundaries>{\n']
            ret.append('    norbitals : ')
            ret.append(repr(self.norbitals))
            ret.append(',\n    nbands : ')
            ret.append(repr(self.nbands))
            ret.append(',\n    nkpoints : ')
            ret.append(repr(self.nkpoints))
            ret.append(',\n    nspins : ')
            ret.append(repr(self.nspins))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    class Orbitals(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=orbitals)
        
        
        Defined at ../src/electronic.f90 lines 66-74
        
        """
        def __init__(self, handle=None):
            """
            self = Orbitals()
            
            
            Defined at ../src/electronic.f90 lines 66-74
            
            
            Returns
            -------
            this : Orbitals
            	Object to be constructed
            
            
            Automatically generated constructor for orbitals
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _OptaPyDOS.f90wrap_orbitals_initialise()
        
        def __del__(self):
            """
            Destructor for class Orbitals
            
            
            Defined at ../src/electronic.f90 lines 66-74
            
            Parameters
            ----------
            this : Orbitals
            	Object to be destructed
            
            
            Automatically generated destructor for orbitals
            """
            if self._alloc:
                _OptaPyDOS.f90wrap_orbitals_finalise(this=self._handle)
        
        @property
        def ion_no(self):
            """
            Element ion_no ftype=integer pytype=int
            
            
            Defined at ../src/electronic.f90 line 67
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _OptaPyDOS.f90wrap_orbitals__array__ion_no(self._handle)
            if array_handle in self._arrays:
                ion_no = self._arrays[array_handle]
            else:
                ion_no = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _OptaPyDOS.f90wrap_orbitals__array__ion_no)
                self._arrays[array_handle] = ion_no
            return ion_no
        
        @ion_no.setter
        def ion_no(self, ion_no):
            self.ion_no[...] = ion_no
        
        @property
        def species_no(self):
            """
            Element species_no ftype=integer pytype=int
            
            
            Defined at ../src/electronic.f90 line 68
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _OptaPyDOS.f90wrap_orbitals__array__species_no(self._handle)
            if array_handle in self._arrays:
                species_no = self._arrays[array_handle]
            else:
                species_no = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _OptaPyDOS.f90wrap_orbitals__array__species_no)
                self._arrays[array_handle] = species_no
            return species_no
        
        @species_no.setter
        def species_no(self, species_no):
            self.species_no[...] = species_no
        
        @property
        def rank_in_species(self):
            """
            Element rank_in_species ftype=integer pytype=int
            
            
            Defined at ../src/electronic.f90 line 69
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _OptaPyDOS.f90wrap_orbitals__array__rank_in_species(self._handle)
            if array_handle in self._arrays:
                rank_in_species = self._arrays[array_handle]
            else:
                rank_in_species = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _OptaPyDOS.f90wrap_orbitals__array__rank_in_species)
                self._arrays[array_handle] = rank_in_species
            return rank_in_species
        
        @rank_in_species.setter
        def rank_in_species(self, rank_in_species):
            self.rank_in_species[...] = rank_in_species
        
        @property
        def am_channel(self):
            """
            Element am_channel ftype=integer pytype=int
            
            
            Defined at ../src/electronic.f90 line 70
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _OptaPyDOS.f90wrap_orbitals__array__am_channel(self._handle)
            if array_handle in self._arrays:
                am_channel = self._arrays[array_handle]
            else:
                am_channel = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _OptaPyDOS.f90wrap_orbitals__array__am_channel)
                self._arrays[array_handle] = am_channel
            return am_channel
        
        @am_channel.setter
        def am_channel(self, am_channel):
            self.am_channel[...] = am_channel
        
        @property
        def shell(self):
            """
            Element shell ftype=integer pytype=int
            
            
            Defined at ../src/electronic.f90 line 71
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _OptaPyDOS.f90wrap_orbitals__array__shell(self._handle)
            if array_handle in self._arrays:
                shell = self._arrays[array_handle]
            else:
                shell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _OptaPyDOS.f90wrap_orbitals__array__shell)
                self._arrays[array_handle] = shell
            return shell
        
        @shell.setter
        def shell(self, shell):
            self.shell[...] = shell
        
        @property
        def am_channel_name(self):
            """
            Element am_channel_name ftype=character(len=10) pytype=str
            
            
            Defined at ../src/electronic.f90 line 72
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _OptaPyDOS.f90wrap_orbitals__array__am_channel_name(self._handle)
            if array_handle in self._arrays:
                am_channel_name = self._arrays[array_handle]
            else:
                am_channel_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _OptaPyDOS.f90wrap_orbitals__array__am_channel_name)
                self._arrays[array_handle] = am_channel_name
            return am_channel_name
        
        @am_channel_name.setter
        def am_channel_name(self, am_channel_name):
            self.am_channel_name[...] = am_channel_name
        
        def __str__(self):
            ret = ['<orbitals>{\n']
            ret.append('    ion_no : ')
            ret.append(repr(self.ion_no))
            ret.append(',\n    species_no : ')
            ret.append(repr(self.species_no))
            ret.append(',\n    rank_in_species : ')
            ret.append(repr(self.rank_in_species))
            ret.append(',\n    am_channel : ')
            ret.append(repr(self.am_channel))
            ret.append(',\n    shell : ')
            ret.append(repr(self.shell))
            ret.append(',\n    am_channel_name : ')
            ret.append(repr(self.am_channel_name))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    @staticmethod
    def elec_report_parameters():
        """
        elec_report_parameters()
        
        
        Defined at ../src/electronic.f90 lines 108-151
        
        
        =========================================================================
         Report the electronic properties in the calculation
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables nbands,num_electrons,nkpoints,kpoint_grid_dim
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by A J Morris Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_elec_report_parameters()
    
    @staticmethod
    def elec_read_band_gradient():
        """
        elec_read_band_gradient()
        
        
        Defined at ../src/electronic.f90 lines 154-282
        
        
        =========================================================================
         Read the .cst_ome file in paralell if appropriate. These are the
         gradients of the bands at each kpoint.
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables: band_gradient,nspins,nbands
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by A J Morris Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_elec_read_band_gradient()
    
    @staticmethod
    def elec_read_optical_mat():
        """
        elec_read_optical_mat()
        
        
        Defined at ../src/electronic.f90 lines 285-429
        
        
        =========================================================================
         Read the .cst_ome file in paralell if appropriate. These are the
         gradients of the bands at each kpoint.
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables: band_gradient,nspins,nbands
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by A J Morris Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_elec_read_optical_mat()
    
    @staticmethod
    def elec_read_band_energy():
        """
        elec_read_band_energy()
        
        
        Defined at ../src/electronic.f90 lines 432-608
        
        
        =========================================================================
         Read the .bands file in the kpoint list, kpoint weights and band energies
         also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables: band_energy, efermi_castep, num_electrons
         spin_polarised, electrons_per_state, nspins,nbands
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by A J Morris Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_elec_read_band_energy()
    
    @staticmethod
    def elec_read_elnes_mat():
        """
        elec_read_elnes_mat()
        
        
        Defined at ../src/electronic.f90 lines 611-852
        
        
        =========================================================================
         Read the .bands file in the kpoint list, kpoint weights and band energies
         also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables: band_energy, efermi_castep, num_electrons
         spin_polarised, electrons_per_state, nspins,nbands
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by A J Morris Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_elec_read_elnes_mat()
    
    @staticmethod
    def elec_pdos_read():
        """
        elec_pdos_read()
        
        
        Defined at ../src/electronic.f90 lines 855-1015
        
        
        =========================================================================
         Read in the full pdos_weights. Write out any variables that we find on the
         way. These will be checked for consistency in the dos module. We can't do it
         yet as we haven't read the bands file.
        -------------------------------------------------------------------------
         Arguments: None
        -------------------------------------------------------------------------
         Parent module variables: pw, pdos_weights, pdos_orbital
        -------------------------------------------------------------------------
         Modules used:  See below
        -------------------------------------------------------------------------
         Key Internal Variables: None
        -------------------------------------------------------------------------
         Necessary conditions: None
        -------------------------------------------------------------------------
         Known Worries: None
        -------------------------------------------------------------------------
         Written by A J Morris Dec 2010
        =========================================================================
        """
        _OptaPyDOS.f90wrap_elec_pdos_read()
    
    @staticmethod
    def elec_dealloc_pdos():
        """
        elec_dealloc_pdos()
        
        
        Defined at ../src/electronic.f90 lines 1017-1034
        
        
        """
        _OptaPyDOS.f90wrap_elec_dealloc_pdos()
    
    @staticmethod
    def elec_dealloc_elnes():
        """
        elec_dealloc_elnes()
        
        
        Defined at ../src/electronic.f90 lines 1036-1052
        
        
        """
        _OptaPyDOS.f90wrap_elec_dealloc_elnes()
    
    @staticmethod
    def elec_dealloc_band_gradient():
        """
        elec_dealloc_band_gradient()
        
        
        Defined at ../src/electronic.f90 lines 1054-1064
        
        
        """
        _OptaPyDOS.f90wrap_elec_dealloc_band_gradient()
    
    @staticmethod
    def elec_dealloc_optical():
        """
        elec_dealloc_optical()
        
        
        Defined at ../src/electronic.f90 lines 1066-1078
        
        
        """
        _OptaPyDOS.f90wrap_elec_dealloc_optical()
    
    @property
    def band_energy(self):
        """
        Element band_energy ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__band_energy(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            band_energy = self._arrays[array_handle]
        else:
            band_energy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__band_energy)
            self._arrays[array_handle] = band_energy
        return band_energy
    
    @band_energy.setter
    def band_energy(self, band_energy):
        self.band_energy[...] = band_energy
    
    @property
    def band_gradient(self):
        """
        Element band_gradient ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__band_gradient(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            band_gradient = self._arrays[array_handle]
        else:
            band_gradient = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__band_gradient)
            self._arrays[array_handle] = band_gradient
        return band_gradient
    
    @band_gradient.setter
    def band_gradient(self, band_gradient):
        self.band_gradient[...] = band_gradient
    
    @property
    def optical_mat(self):
        """
        Element optical_mat ftype=complex(kind=dp) pytype=complex
        
        
        Defined at ../src/electronic.f90 line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__optical_mat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            optical_mat = self._arrays[array_handle]
        else:
            optical_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__optical_mat)
            self._arrays[array_handle] = optical_mat
        return optical_mat
    
    @optical_mat.setter
    def optical_mat(self, optical_mat):
        self.optical_mat[...] = optical_mat
    
    @property
    def elnes_mat(self):
        """
        Element elnes_mat ftype=complex(kind=dp) pytype=complex
        
        
        Defined at ../src/electronic.f90 line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__elnes_mat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            elnes_mat = self._arrays[array_handle]
        else:
            elnes_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__elnes_mat)
            self._arrays[array_handle] = elnes_mat
        return elnes_mat
    
    @elnes_mat.setter
    def elnes_mat(self, elnes_mat):
        self.elnes_mat[...] = elnes_mat
    
    @property
    def efermi(self):
        """
        Element efermi ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 43
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__efermi()
    
    @efermi.setter
    def efermi(self, efermi):
        _OptaPyDOS.f90wrap_od_electronic__set__efermi(efermi)
    
    @property
    def efermi_set(self):
        """
        Element efermi_set ftype=logical pytype=bool
        
        
        Defined at ../src/electronic.f90 line 44
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__efermi_set()
    
    @efermi_set.setter
    def efermi_set(self, efermi_set):
        _OptaPyDOS.f90wrap_od_electronic__set__efermi_set(efermi_set)
    
    @property
    def unshifted_efermi(self):
        """
        Element unshifted_efermi ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 45
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__unshifted_efermi()
    
    @unshifted_efermi.setter
    def unshifted_efermi(self, unshifted_efermi):
        _OptaPyDOS.f90wrap_od_electronic__set__unshifted_efermi(unshifted_efermi)
    
    @property
    def efermi_castep(self):
        """
        Element efermi_castep ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 46
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__efermi_castep()
    
    @efermi_castep.setter
    def efermi_castep(self, efermi_castep):
        _OptaPyDOS.f90wrap_od_electronic__set__efermi_castep(efermi_castep)
    
    @property
    def num_electrons(self):
        """
        Element num_electrons ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 48
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__num_electrons(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            num_electrons = self._arrays[array_handle]
        else:
            num_electrons = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__num_electrons)
            self._arrays[array_handle] = num_electrons
        return num_electrons
    
    @num_electrons.setter
    def num_electrons(self, num_electrons):
        self.num_electrons[...] = num_electrons
    
    @property
    def nbands(self):
        """
        Element nbands ftype=integer pytype=int
        
        
        Defined at ../src/electronic.f90 line 51
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__nbands()
    
    @nbands.setter
    def nbands(self, nbands):
        _OptaPyDOS.f90wrap_od_electronic__set__nbands(nbands)
    
    @property
    def nspins(self):
        """
        Element nspins ftype=integer pytype=int
        
        
        Defined at ../src/electronic.f90 line 51
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__nspins()
    
    @nspins.setter
    def nspins(self, nspins):
        _OptaPyDOS.f90wrap_od_electronic__set__nspins(nspins)
    
    @property
    def electrons_per_state(self):
        """
        Element electrons_per_state ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 52
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__electrons_per_state()
    
    @electrons_per_state.setter
    def electrons_per_state(self, electrons_per_state):
        _OptaPyDOS.f90wrap_od_electronic__set__electrons_per_state(electrons_per_state)
    
    @property
    def spin_polarised(self):
        """
        Element spin_polarised ftype=logical pytype=bool
        
        
        Defined at ../src/electronic.f90 line 54
        
        """
        return _OptaPyDOS.f90wrap_od_electronic__get__spin_polarised()
    
    @spin_polarised.setter
    def spin_polarised(self, spin_polarised):
        _OptaPyDOS.f90wrap_od_electronic__set__spin_polarised(spin_polarised)
    
    @property
    def pdos_weights(self):
        """
        Element pdos_weights ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 77
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__pdos_weights(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pdos_weights = self._arrays[array_handle]
        else:
            pdos_weights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__pdos_weights)
            self._arrays[array_handle] = pdos_weights
        return pdos_weights
    
    @pdos_weights.setter
    def pdos_weights(self, pdos_weights):
        self.pdos_weights[...] = pdos_weights
    
    @property
    def all_kpoints(self):
        """
        Element all_kpoints ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/electronic.f90 line 82
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_electronic__array__all_kpoints(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            all_kpoints = self._arrays[array_handle]
        else:
            all_kpoints = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_electronic__array__all_kpoints)
            self._arrays[array_handle] = all_kpoints
        return all_kpoints
    
    @all_kpoints.setter
    def all_kpoints(self, all_kpoints):
        self.all_kpoints[...] = all_kpoints
    
    def __str__(self):
        ret = ['<od_electronic>{\n']
        ret.append('    band_energy : ')
        ret.append(repr(self.band_energy))
        ret.append(',\n    band_gradient : ')
        ret.append(repr(self.band_gradient))
        ret.append(',\n    optical_mat : ')
        ret.append(repr(self.optical_mat))
        ret.append(',\n    elnes_mat : ')
        ret.append(repr(self.elnes_mat))
        ret.append(',\n    efermi : ')
        ret.append(repr(self.efermi))
        ret.append(',\n    efermi_set : ')
        ret.append(repr(self.efermi_set))
        ret.append(',\n    unshifted_efermi : ')
        ret.append(repr(self.unshifted_efermi))
        ret.append(',\n    efermi_castep : ')
        ret.append(repr(self.efermi_castep))
        ret.append(',\n    num_electrons : ')
        ret.append(repr(self.num_electrons))
        ret.append(',\n    nbands : ')
        ret.append(repr(self.nbands))
        ret.append(',\n    nspins : ')
        ret.append(repr(self.nspins))
        ret.append(',\n    electrons_per_state : ')
        ret.append(repr(self.electrons_per_state))
        ret.append(',\n    spin_polarised : ')
        ret.append(repr(self.spin_polarised))
        ret.append(',\n    pdos_weights : ')
        ret.append(repr(self.pdos_weights))
        ret.append(',\n    all_kpoints : ')
        ret.append(repr(self.all_kpoints))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_electronic = Od_Electronic()

class Od_Io(f90wrap.runtime.FortranModule):
    """
    Module od_io
    
    
    Defined at ../src/io.f90 lines 36-187
    
    """
    @staticmethod
    def io_get_seedname():
        """
        io_get_seedname()
        
        
        Defined at ../src/io.f90 lines 59-85
        
        
        ==================================================================
         Get the seedname from the commandline
         Note iargc and getarg are not standard
         Some platforms require them to be external or provide
         equivalent routines. Not a problem in f2003
        ===================================================================
        """
        _OptaPyDOS.f90wrap_io_get_seedname()
    
    @staticmethod
    def io_error(error_msg):
        """
        io_error(error_msg)
        
        
        Defined at ../src/io.f90 lines 88-105
        
        Parameters
        ----------
        error_msg : str
        
        ==================================================================
         Aborts giving error message
        ===================================================================
        """
        _OptaPyDOS.f90wrap_io_error(error_msg=error_msg)
    
    @staticmethod
    def io_date():
        """
        cdate, ctime = io_date()
        
        
        Defined at ../src/io.f90 lines 108-130
        
        
        Returns
        -------
        cdate : str
        ctime : str
        
        ==================================================================
             Returns two strings containing the date and the time
             in human-readable format. Uses a standard f90 call.
        ===================================================================
        """
        cdate, ctime = _OptaPyDOS.f90wrap_io_date()
        return cdate, ctime
    
    @staticmethod
    def io_time():
        """
        io_time = io_time()
        
        
        Defined at ../src/io.f90 lines 133-161
        
        
        Returns
        -------
        io_time : float
        
        ==================================================================
         Returns elapsed CPU time in seconds since its first call
         uses standard f90 call
        ===================================================================
        """
        io_time = _OptaPyDOS.f90wrap_io_time()
        return io_time
    
    @staticmethod
    def io_file_unit():
        """
        io_file_unit = io_file_unit()
        
        
        Defined at ../src/io.f90 lines 164-187
        
        
        Returns
        -------
        io_file_unit : int
        
        ==================================================================
         Returns an unused unit number
         (so we can open a file on that unit
        ===================================================================
        """
        io_file_unit = _OptaPyDOS.f90wrap_io_file_unit()
        return io_file_unit
    
    @property
    def stdout(self):
        """
        Element stdout ftype=integer pytype=int
        
        
        Defined at ../src/io.f90 line 44
        
        """
        return _OptaPyDOS.f90wrap_od_io__get__stdout()
    
    @stdout.setter
    def stdout(self, stdout):
        _OptaPyDOS.f90wrap_od_io__set__stdout(stdout)
    
    @property
    def stderr(self):
        """
        Element stderr ftype=integer pytype=int
        
        
        Defined at ../src/io.f90 line 45
        
        """
        return _OptaPyDOS.f90wrap_od_io__get__stderr()
    
    @stderr.setter
    def stderr(self, stderr):
        _OptaPyDOS.f90wrap_od_io__set__stderr(stderr)
    
    @property
    def seedname(self):
        """
        Element seedname ftype=character(len=50) pytype=str
        
        
        Defined at ../src/io.f90 line 46
        
        """
        return _OptaPyDOS.f90wrap_od_io__get__seedname()
    
    @seedname.setter
    def seedname(self, seedname):
        _OptaPyDOS.f90wrap_od_io__set__seedname(seedname)
    
    @property
    def maxlen(self):
        """
        Element maxlen ftype=integer pytype=int
        
        
        Defined at ../src/io.f90 line 47
        
        """
        return _OptaPyDOS.f90wrap_od_io__get__maxlen()
    
    @property
    def filename_len(self):
        """
        Element filename_len ftype=integer pytype=int
        
        
        Defined at ../src/io.f90 line 48
        
        """
        return _OptaPyDOS.f90wrap_od_io__get__filename_len()
    
    def __str__(self):
        ret = ['<od_io>{\n']
        ret.append('    stdout : ')
        ret.append(repr(self.stdout))
        ret.append(',\n    stderr : ')
        ret.append(repr(self.stderr))
        ret.append(',\n    seedname : ')
        ret.append(repr(self.seedname))
        ret.append(',\n    maxlen : ')
        ret.append(repr(self.maxlen))
        ret.append(',\n    filename_len : ')
        ret.append(repr(self.filename_len))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_io = Od_Io()

class Od_Jdos(f90wrap.runtime.FortranModule):
    """
    Module od_jdos
    
    
    Defined at ../src/jdos.f90 lines 24-206
    
    """
    @staticmethod
    def jdos_calculate():
        """
        jdos_calculate()
        
        
        Defined at ../src/jdos.f90 lines 35-74
        
        
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_jdos_calculate()
    
    @staticmethod
    def write_jdos(e, dos, dos_name):
        """
        write_jdos(e, dos, dos_name)
        
        
        Defined at ../src/jdos.f90 lines 77-149
        
        Parameters
        ----------
        e : float array
        dos : float array
        dos_name : str
        
        ===============================================================================
         This routine receives an energy scale, a density of states and a file name
         and writes out the DOS to disk
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_write_jdos(e=e, dos=dos, dos_name=dos_name)
    
    @staticmethod
    def write_jdos_xmgrace(dos_name, e, dos):
        """
        write_jdos_xmgrace(dos_name, e, dos)
        
        
        Defined at ../src/jdos.f90 lines 154-206
        
        Parameters
        ----------
        dos_name : str
        e : float array
        dos : float array
        
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_write_jdos_xmgrace(dos_name=dos_name, e=e, dos=dos)
    
    _dt_array_initialisers = []
    

od_jdos = Od_Jdos()

class Od_Jdos_Utils(f90wrap.runtime.FortranModule):
    """
    Module od_jdos_utils
    
    
    Defined at ../src/jdos_utils.f90 lines 31-465
    
    """
    @staticmethod
    def jdos_utils_calculate(matrix_weights=None):
        """
        jdos_utils_calculate([matrix_weights])
        
        
        Defined at ../src/jdos_utils.f90 lines 64-184
        
        Parameters
        ----------
        matrix_weights : float array
        
        ===============================================================================
         Main routine in dos module, drives the calculation of Density of states for
         both task : dos and also if it is required elsewhere.
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_jdos_utils_calculate(matrix_weights=matrix_weights)
    
    @staticmethod
    def setup_energy_scale():
        """
        setup_energy_scale()
        
        
        Defined at ../src/jdos_utils.f90 lines 187-243
        
        
        ===============================================================================
         Sets up all broadening independent DOS concerns
         Calls the relevant dos calculator.
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_setup_energy_scale()
    
    @staticmethod
    def jdos_deallocate():
        """
        jdos_deallocate()
        
        
        Defined at ../src/jdos_utils.f90 lines 265-291
        
        
        ===============================================================================
        ===============================================================================
        """
        _OptaPyDOS.f90wrap_jdos_deallocate()
    
    @property
    def jdos_adaptive(self):
        """
        Element jdos_adaptive ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/jdos_utils.f90 line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_jdos_utils__array__jdos_adaptive(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jdos_adaptive = self._arrays[array_handle]
        else:
            jdos_adaptive = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_jdos_utils__array__jdos_adaptive)
            self._arrays[array_handle] = jdos_adaptive
        return jdos_adaptive
    
    @jdos_adaptive.setter
    def jdos_adaptive(self, jdos_adaptive):
        self.jdos_adaptive[...] = jdos_adaptive
    
    @property
    def jdos_fixed(self):
        """
        Element jdos_fixed ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/jdos_utils.f90 line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_jdos_utils__array__jdos_fixed(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jdos_fixed = self._arrays[array_handle]
        else:
            jdos_fixed = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_jdos_utils__array__jdos_fixed)
            self._arrays[array_handle] = jdos_fixed
        return jdos_fixed
    
    @jdos_fixed.setter
    def jdos_fixed(self, jdos_fixed):
        self.jdos_fixed[...] = jdos_fixed
    
    @property
    def jdos_linear(self):
        """
        Element jdos_linear ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/jdos_utils.f90 line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_jdos_utils__array__jdos_linear(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            jdos_linear = self._arrays[array_handle]
        else:
            jdos_linear = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_jdos_utils__array__jdos_linear)
            self._arrays[array_handle] = jdos_linear
        return jdos_linear
    
    @jdos_linear.setter
    def jdos_linear(self, jdos_linear):
        self.jdos_linear[...] = jdos_linear
    
    @property
    def jdos_nbins(self):
        """
        Element jdos_nbins ftype=integer pytype=int
        
        
        Defined at ../src/jdos_utils.f90 line 44
        
        """
        return _OptaPyDOS.f90wrap_od_jdos_utils__get__jdos_nbins()
    
    @jdos_nbins.setter
    def jdos_nbins(self, jdos_nbins):
        _OptaPyDOS.f90wrap_od_jdos_utils__set__jdos_nbins(jdos_nbins)
    
    @property
    def e(self):
        """
        Element e ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/jdos_utils.f90 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_jdos_utils__array__e(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            e = self._arrays[array_handle]
        else:
            e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_jdos_utils__array__e)
            self._arrays[array_handle] = e
        return e
    
    @e.setter
    def e(self, e):
        self.e[...] = e
    
    @property
    def delta_bins(self):
        """
        Element delta_bins ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/jdos_utils.f90 line 56
        
        """
        return _OptaPyDOS.f90wrap_od_jdos_utils__get__delta_bins()
    
    @delta_bins.setter
    def delta_bins(self, delta_bins):
        _OptaPyDOS.f90wrap_od_jdos_utils__set__delta_bins(delta_bins)
    
    @property
    def calc_weighted_jdos(self):
        """
        Element calc_weighted_jdos ftype=logical pytype=bool
        
        
        Defined at ../src/jdos_utils.f90 line 57
        
        """
        return _OptaPyDOS.f90wrap_od_jdos_utils__get__calc_weighted_jdos()
    
    @calc_weighted_jdos.setter
    def calc_weighted_jdos(self, calc_weighted_jdos):
        _OptaPyDOS.f90wrap_od_jdos_utils__set__calc_weighted_jdos(calc_weighted_jdos)
    
    @property
    def vb_max(self):
        """
        Element vb_max ftype=integer pytype=int
        
        
        Defined at ../src/jdos_utils.f90 line 58
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_jdos_utils__array__vb_max(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            vb_max = self._arrays[array_handle]
        else:
            vb_max = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_jdos_utils__array__vb_max)
            self._arrays[array_handle] = vb_max
        return vb_max
    
    @vb_max.setter
    def vb_max(self, vb_max):
        self.vb_max[...] = vb_max
    
    def __str__(self):
        ret = ['<od_jdos_utils>{\n']
        ret.append('    jdos_adaptive : ')
        ret.append(repr(self.jdos_adaptive))
        ret.append(',\n    jdos_fixed : ')
        ret.append(repr(self.jdos_fixed))
        ret.append(',\n    jdos_linear : ')
        ret.append(repr(self.jdos_linear))
        ret.append(',\n    jdos_nbins : ')
        ret.append(repr(self.jdos_nbins))
        ret.append(',\n    e : ')
        ret.append(repr(self.e))
        ret.append(',\n    delta_bins : ')
        ret.append(repr(self.delta_bins))
        ret.append(',\n    calc_weighted_jdos : ')
        ret.append(repr(self.calc_weighted_jdos))
        ret.append(',\n    vb_max : ')
        ret.append(repr(self.vb_max))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_jdos_utils = Od_Jdos_Utils()

class Od_Optics(f90wrap.runtime.FortranModule):
    """
    Module od_optics
    
    
    Defined at ../src/optics.f90 lines 24-1421
    
    """
    @staticmethod
    def optics_calculate():
        """
        optics_calculate()
        
        
        Defined at ../src/optics.f90 lines 78-155
        
        
        """
        _OptaPyDOS.f90wrap_optics_calculate()
    
    @property
    def matrix_weights(self):
        """
        Element matrix_weights ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__matrix_weights(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            matrix_weights = self._arrays[array_handle]
        else:
            matrix_weights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__matrix_weights)
            self._arrays[array_handle] = matrix_weights
        return matrix_weights
    
    @matrix_weights.setter
    def matrix_weights(self, matrix_weights):
        self.matrix_weights[...] = matrix_weights
    
    @property
    def dos_matrix_weights(self):
        """
        Element dos_matrix_weights ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 47
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__dos_matrix_weights(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dos_matrix_weights = self._arrays[array_handle]
        else:
            dos_matrix_weights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__dos_matrix_weights)
            self._arrays[array_handle] = dos_matrix_weights
        return dos_matrix_weights
    
    @dos_matrix_weights.setter
    def dos_matrix_weights(self, dos_matrix_weights):
        self.dos_matrix_weights[...] = dos_matrix_weights
    
    @property
    def weighted_jdos(self):
        """
        Element weighted_jdos ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 48
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__weighted_jdos(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            weighted_jdos = self._arrays[array_handle]
        else:
            weighted_jdos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__weighted_jdos)
            self._arrays[array_handle] = weighted_jdos
        return weighted_jdos
    
    @weighted_jdos.setter
    def weighted_jdos(self, weighted_jdos):
        self.weighted_jdos[...] = weighted_jdos
    
    @property
    def weighted_dos_at_e(self):
        """
        Element weighted_dos_at_e ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__weighted_dos_at_e(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            weighted_dos_at_e = self._arrays[array_handle]
        else:
            weighted_dos_at_e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__weighted_dos_at_e)
            self._arrays[array_handle] = weighted_dos_at_e
        return weighted_dos_at_e
    
    @weighted_dos_at_e.setter
    def weighted_dos_at_e(self, weighted_dos_at_e):
        self.weighted_dos_at_e[...] = weighted_dos_at_e
    
    @property
    def dos_at_e(self):
        """
        Element dos_at_e ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 50
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__dos_at_e(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dos_at_e = self._arrays[array_handle]
        else:
            dos_at_e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__dos_at_e)
            self._arrays[array_handle] = dos_at_e
        return dos_at_e
    
    @dos_at_e.setter
    def dos_at_e(self, dos_at_e):
        self.dos_at_e[...] = dos_at_e
    
    @property
    def epsilon(self):
        """
        Element epsilon ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__epsilon(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            epsilon = self._arrays[array_handle]
        else:
            epsilon = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__epsilon)
            self._arrays[array_handle] = epsilon
        return epsilon
    
    @epsilon.setter
    def epsilon(self, epsilon):
        self.epsilon[...] = epsilon
    
    @property
    def conduct(self):
        """
        Element conduct ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__conduct(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            conduct = self._arrays[array_handle]
        else:
            conduct = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__conduct)
            self._arrays[array_handle] = conduct
        return conduct
    
    @conduct.setter
    def conduct(self, conduct):
        self.conduct[...] = conduct
    
    @property
    def refract(self):
        """
        Element refract ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__refract(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            refract = self._arrays[array_handle]
        else:
            refract = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__refract)
            self._arrays[array_handle] = refract
        return refract
    
    @refract.setter
    def refract(self, refract):
        self.refract[...] = refract
    
    @property
    def loss_fn(self):
        """
        Element loss_fn ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 55
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__loss_fn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            loss_fn = self._arrays[array_handle]
        else:
            loss_fn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__loss_fn)
            self._arrays[array_handle] = loss_fn
        return loss_fn
    
    @loss_fn.setter
    def loss_fn(self, loss_fn):
        self.loss_fn[...] = loss_fn
    
    @property
    def absorp(self):
        """
        Element absorp ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__absorp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            absorp = self._arrays[array_handle]
        else:
            absorp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__absorp)
            self._arrays[array_handle] = absorp
        return absorp
    
    @absorp.setter
    def absorp(self, absorp):
        self.absorp[...] = absorp
    
    @property
    def reflect(self):
        """
        Element reflect ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 57
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__reflect(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            reflect = self._arrays[array_handle]
        else:
            reflect = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__reflect)
            self._arrays[array_handle] = reflect
        return reflect
    
    @reflect.setter
    def reflect(self, reflect):
        self.reflect[...] = reflect
    
    @property
    def intra(self):
        """
        Element intra ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/optics.f90 line 59
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_optics__array__intra(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            intra = self._arrays[array_handle]
        else:
            intra = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_optics__array__intra)
            self._arrays[array_handle] = intra
        return intra
    
    @intra.setter
    def intra(self, intra):
        self.intra[...] = intra
    
    def __str__(self):
        ret = ['<od_optics>{\n']
        ret.append('    matrix_weights : ')
        ret.append(repr(self.matrix_weights))
        ret.append(',\n    dos_matrix_weights : ')
        ret.append(repr(self.dos_matrix_weights))
        ret.append(',\n    weighted_jdos : ')
        ret.append(repr(self.weighted_jdos))
        ret.append(',\n    weighted_dos_at_e : ')
        ret.append(repr(self.weighted_dos_at_e))
        ret.append(',\n    dos_at_e : ')
        ret.append(repr(self.dos_at_e))
        ret.append(',\n    epsilon : ')
        ret.append(repr(self.epsilon))
        ret.append(',\n    conduct : ')
        ret.append(repr(self.conduct))
        ret.append(',\n    refract : ')
        ret.append(repr(self.refract))
        ret.append(',\n    loss_fn : ')
        ret.append(repr(self.loss_fn))
        ret.append(',\n    absorp : ')
        ret.append(repr(self.absorp))
        ret.append(',\n    reflect : ')
        ret.append(repr(self.reflect))
        ret.append(',\n    intra : ')
        ret.append(repr(self.intra))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_optics = Od_Optics()

class Od_Parameters(f90wrap.runtime.FortranModule):
    """
    Module od_parameters
    
    
    Defined at ../src/parameters.f90 lines 35-1576
    
    """
    @staticmethod
    def param_read():
        """
        param_read()
        
        
        Defined at ../src/parameters.f90 lines 143-458
        
        
        ==================================================================
         Read parameters and calculate derived values
        ===================================================================
        """
        _OptaPyDOS.f90wrap_param_read()
    
    @staticmethod
    def param_write_header():
        """
        param_write_header()
        
        
        Defined at ../src/parameters.f90 lines 497-535
        
        
        """
        _OptaPyDOS.f90wrap_param_write_header()
    
    @staticmethod
    def param_write_atomic_coord():
        """
        param_write_atomic_coord()
        
        
        Defined at ../src/parameters.f90 lines 537-586
        
        
        ==================================================================
         write atomic coodes to stdout
        ===================================================================
        """
        _OptaPyDOS.f90wrap_param_write_atomic_coord()
    
    @staticmethod
    def param_write():
        """
        param_write()
        
        
        Defined at ../src/parameters.f90 lines 589-804
        
        
        ==================================================================
         write parameters to stdout
        ===================================================================
        """
        _OptaPyDOS.f90wrap_param_write()
    
    @staticmethod
    def param_dealloc():
        """
        param_dealloc()
        
        
        Defined at ../src/parameters.f90 lines 807-824
        
        
        ==================================================================
         release memory from allocated parameters
        ===================================================================
        """
        _OptaPyDOS.f90wrap_param_dealloc()
    
    @staticmethod
    def param_dist():
        """
        param_dist()
        
        
        Defined at ../src/parameters.f90 lines 1501-1576
        
        
        -----------------------------------------------------
         Send the parameters from the root node to all others
        -----------------------------------------------------
        """
        _OptaPyDOS.f90wrap_param_dist()
    
    @property
    def output_format(self):
        """
        Element output_format ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 45
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__output_format()
    
    @output_format.setter
    def output_format(self, output_format):
        _OptaPyDOS.f90wrap_od_parameters__set__output_format(output_format)
    
    @property
    def devel_flag(self):
        """
        Element devel_flag ftype=character(len=100) pytype=str
        
        
        Defined at ../src/parameters.f90 line 46
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__devel_flag()
    
    @devel_flag.setter
    def devel_flag(self, devel_flag):
        _OptaPyDOS.f90wrap_od_parameters__set__devel_flag(devel_flag)
    
    @property
    def iprint(self):
        """
        Element iprint ftype=integer pytype=int
        
        
        Defined at ../src/parameters.f90 line 47
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__iprint()
    
    @iprint.setter
    def iprint(self, iprint):
        _OptaPyDOS.f90wrap_od_parameters__set__iprint(iprint)
    
    @property
    def energy_unit(self):
        """
        Element energy_unit ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 48
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__energy_unit()
    
    @energy_unit.setter
    def energy_unit(self, energy_unit):
        _OptaPyDOS.f90wrap_od_parameters__set__energy_unit(energy_unit)
    
    @property
    def length_unit(self):
        """
        Element length_unit ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 49
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__length_unit()
    
    @length_unit.setter
    def length_unit(self, length_unit):
        _OptaPyDOS.f90wrap_od_parameters__set__length_unit(length_unit)
    
    @property
    def legacy_file_format(self):
        """
        Element legacy_file_format ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 50
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__legacy_file_format()
    
    @legacy_file_format.setter
    def legacy_file_format(self, legacy_file_format):
        _OptaPyDOS.f90wrap_od_parameters__set__legacy_file_format(legacy_file_format)
    
    @property
    def kpoint_mp_grid(self):
        """
        Element kpoint_mp_grid ftype=integer pytype=int
        
        
        Defined at ../src/parameters.f90 line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_parameters__array__kpoint_mp_grid(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kpoint_mp_grid = self._arrays[array_handle]
        else:
            kpoint_mp_grid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_parameters__array__kpoint_mp_grid)
            self._arrays[array_handle] = kpoint_mp_grid
        return kpoint_mp_grid
    
    @kpoint_mp_grid.setter
    def kpoint_mp_grid(self, kpoint_mp_grid):
        self.kpoint_mp_grid[...] = kpoint_mp_grid
    
    @property
    def dos(self):
        """
        Element dos ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 55
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__dos()
    
    @dos.setter
    def dos(self, dos):
        _OptaPyDOS.f90wrap_od_parameters__set__dos(dos)
    
    @property
    def compare_dos(self):
        """
        Element compare_dos ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 56
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__compare_dos()
    
    @compare_dos.setter
    def compare_dos(self, compare_dos):
        _OptaPyDOS.f90wrap_od_parameters__set__compare_dos(compare_dos)
    
    @property
    def pdos(self):
        """
        Element pdos ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 57
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__pdos()
    
    @pdos.setter
    def pdos(self, pdos):
        _OptaPyDOS.f90wrap_od_parameters__set__pdos(pdos)
    
    @property
    def jdos(self):
        """
        Element jdos ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 58
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__jdos()
    
    @jdos.setter
    def jdos(self, jdos):
        _OptaPyDOS.f90wrap_od_parameters__set__jdos(jdos)
    
    @property
    def compare_jdos(self):
        """
        Element compare_jdos ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 59
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__compare_jdos()
    
    @compare_jdos.setter
    def compare_jdos(self, compare_jdos):
        _OptaPyDOS.f90wrap_od_parameters__set__compare_jdos(compare_jdos)
    
    @property
    def optics(self):
        """
        Element optics ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 60
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__optics()
    
    @optics.setter
    def optics(self, optics):
        _OptaPyDOS.f90wrap_od_parameters__set__optics(optics)
    
    @property
    def core(self):
        """
        Element core ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 61
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__core()
    
    @core.setter
    def core(self, core):
        _OptaPyDOS.f90wrap_od_parameters__set__core(core)
    
    @property
    def fixed(self):
        """
        Element fixed ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 64
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__fixed()
    
    @fixed.setter
    def fixed(self, fixed):
        _OptaPyDOS.f90wrap_od_parameters__set__fixed(fixed)
    
    @property
    def adaptive(self):
        """
        Element adaptive ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 65
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__adaptive()
    
    @adaptive.setter
    def adaptive(self, adaptive):
        _OptaPyDOS.f90wrap_od_parameters__set__adaptive(adaptive)
    
    @property
    def linear(self):
        """
        Element linear ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 66
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__linear()
    
    @linear.setter
    def linear(self, linear):
        _OptaPyDOS.f90wrap_od_parameters__set__linear(linear)
    
    @property
    def quad(self):
        """
        Element quad ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 67
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__quad()
    
    @quad.setter
    def quad(self, quad):
        _OptaPyDOS.f90wrap_od_parameters__set__quad(quad)
    
    @property
    def compute_band_energy(self):
        """
        Element compute_band_energy ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 70
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__compute_band_energy()
    
    @compute_band_energy.setter
    def compute_band_energy(self, compute_band_energy):
        _OptaPyDOS.f90wrap_od_parameters__set__compute_band_energy(compute_band_energy)
    
    @property
    def adaptive_smearing(self):
        """
        Element adaptive_smearing ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 71
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__adaptive_smearing()
    
    @adaptive_smearing.setter
    def adaptive_smearing(self, adaptive_smearing):
        _OptaPyDOS.f90wrap_od_parameters__set__adaptive_smearing(adaptive_smearing)
    
    @property
    def fixed_smearing(self):
        """
        Element fixed_smearing ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 72
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__fixed_smearing()
    
    @fixed_smearing.setter
    def fixed_smearing(self, fixed_smearing):
        _OptaPyDOS.f90wrap_od_parameters__set__fixed_smearing(fixed_smearing)
    
    @property
    def linear_smearing(self):
        """
        Element linear_smearing ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 73
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__linear_smearing()
    
    @linear_smearing.setter
    def linear_smearing(self, linear_smearing):
        _OptaPyDOS.f90wrap_od_parameters__set__linear_smearing(linear_smearing)
    
    @property
    def dos_per_volume(self):
        """
        Element dos_per_volume ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 74
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__dos_per_volume()
    
    @dos_per_volume.setter
    def dos_per_volume(self, dos_per_volume):
        _OptaPyDOS.f90wrap_od_parameters__set__dos_per_volume(dos_per_volume)
    
    @property
    def efermi_user(self):
        """
        Element efermi_user ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 75
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__efermi_user()
    
    @efermi_user.setter
    def efermi_user(self, efermi_user):
        _OptaPyDOS.f90wrap_od_parameters__set__efermi_user(efermi_user)
    
    @property
    def efermi_choice(self):
        """
        Element efermi_choice ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 76
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__efermi_choice()
    
    @efermi_choice.setter
    def efermi_choice(self, efermi_choice):
        _OptaPyDOS.f90wrap_od_parameters__set__efermi_choice(efermi_choice)
    
    @property
    def finite_bin_correction(self):
        """
        Element finite_bin_correction ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 77
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__finite_bin_correction()
    
    @finite_bin_correction.setter
    def finite_bin_correction(self, finite_bin_correction):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__finite_bin_correction(finite_bin_correction)
    
    @property
    def hybrid_linear(self):
        """
        Element hybrid_linear ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 78
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__hybrid_linear()
    
    @hybrid_linear.setter
    def hybrid_linear(self, hybrid_linear):
        _OptaPyDOS.f90wrap_od_parameters__set__hybrid_linear(hybrid_linear)
    
    @property
    def hybrid_linear_grad_tol(self):
        """
        Element hybrid_linear_grad_tol ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 79
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__hybrid_linear_grad_tol()
    
    @hybrid_linear_grad_tol.setter
    def hybrid_linear_grad_tol(self, hybrid_linear_grad_tol):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__hybrid_linear_grad_tol(hybrid_linear_grad_tol)
    
    @property
    def numerical_intdos(self):
        """
        Element numerical_intdos ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 80
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__numerical_intdos()
    
    @numerical_intdos.setter
    def numerical_intdos(self, numerical_intdos):
        _OptaPyDOS.f90wrap_od_parameters__set__numerical_intdos(numerical_intdos)
    
    @property
    def compute_band_gap(self):
        """
        Element compute_band_gap ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 81
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__compute_band_gap()
    
    @compute_band_gap.setter
    def compute_band_gap(self, compute_band_gap):
        _OptaPyDOS.f90wrap_od_parameters__set__compute_band_gap(compute_band_gap)
    
    @property
    def set_efermi_zero(self):
        """
        Element set_efermi_zero ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 83
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__set_efermi_zero()
    
    @set_efermi_zero.setter
    def set_efermi_zero(self, set_efermi_zero):
        _OptaPyDOS.f90wrap_od_parameters__set__set_efermi_zero(set_efermi_zero)
    
    @property
    def dos_min_energy(self):
        """
        Element dos_min_energy ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 84
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__dos_min_energy()
    
    @dos_min_energy.setter
    def dos_min_energy(self, dos_min_energy):
        _OptaPyDOS.f90wrap_od_parameters__set__dos_min_energy(dos_min_energy)
    
    @property
    def dos_max_energy(self):
        """
        Element dos_max_energy ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 85
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__dos_max_energy()
    
    @dos_max_energy.setter
    def dos_max_energy(self, dos_max_energy):
        _OptaPyDOS.f90wrap_od_parameters__set__dos_max_energy(dos_max_energy)
    
    @property
    def dos_spacing(self):
        """
        Element dos_spacing ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 86
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__dos_spacing()
    
    @dos_spacing.setter
    def dos_spacing(self, dos_spacing):
        _OptaPyDOS.f90wrap_od_parameters__set__dos_spacing(dos_spacing)
    
    @property
    def dos_nbins(self):
        """
        Element dos_nbins ftype=integer pytype=int
        
        
        Defined at ../src/parameters.f90 line 87
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__dos_nbins()
    
    @dos_nbins.setter
    def dos_nbins(self, dos_nbins):
        _OptaPyDOS.f90wrap_od_parameters__set__dos_nbins(dos_nbins)
    
    @property
    def pdos_string(self):
        """
        Element pdos_string ftype=character(len=maxlen) pytype=str
        
        
        Defined at ../src/parameters.f90 line 90
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__pdos_string()
    
    @pdos_string.setter
    def pdos_string(self, pdos_string):
        _OptaPyDOS.f90wrap_od_parameters__set__pdos_string(pdos_string)
    
    @property
    def jdos_max_energy(self):
        """
        Element jdos_max_energy ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 94
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__jdos_max_energy()
    
    @jdos_max_energy.setter
    def jdos_max_energy(self, jdos_max_energy):
        _OptaPyDOS.f90wrap_od_parameters__set__jdos_max_energy(jdos_max_energy)
    
    @property
    def jdos_spacing(self):
        """
        Element jdos_spacing ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 95
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__jdos_spacing()
    
    @jdos_spacing.setter
    def jdos_spacing(self, jdos_spacing):
        _OptaPyDOS.f90wrap_od_parameters__set__jdos_spacing(jdos_spacing)
    
    @property
    def scissor_op(self):
        """
        Element scissor_op ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 96
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__scissor_op()
    
    @scissor_op.setter
    def scissor_op(self, scissor_op):
        _OptaPyDOS.f90wrap_od_parameters__set__scissor_op(scissor_op)
    
    @property
    def exclude_bands(self):
        """
        Element exclude_bands ftype=integer pytype=int
        
        
        Defined at ../src/parameters.f90 line 97
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_parameters__array__exclude_bands(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            exclude_bands = self._arrays[array_handle]
        else:
            exclude_bands = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_parameters__array__exclude_bands)
            self._arrays[array_handle] = exclude_bands
        return exclude_bands
    
    @exclude_bands.setter
    def exclude_bands(self, exclude_bands):
        self.exclude_bands[...] = exclude_bands
    
    @property
    def num_exclude_bands(self):
        """
        Element num_exclude_bands ftype=integer pytype=int
        
        
        Defined at ../src/parameters.f90 line 99
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__num_exclude_bands()
    
    @num_exclude_bands.setter
    def num_exclude_bands(self, num_exclude_bands):
        _OptaPyDOS.f90wrap_od_parameters__set__num_exclude_bands(num_exclude_bands)
    
    @property
    def optics_geom(self):
        """
        Element optics_geom ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 102
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__optics_geom()
    
    @optics_geom.setter
    def optics_geom(self, optics_geom):
        _OptaPyDOS.f90wrap_od_parameters__set__optics_geom(optics_geom)
    
    @property
    def optics_qdir(self):
        """
        Element optics_qdir ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 103
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_parameters__array__optics_qdir(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            optics_qdir = self._arrays[array_handle]
        else:
            optics_qdir = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_parameters__array__optics_qdir)
            self._arrays[array_handle] = optics_qdir
        return optics_qdir
    
    @optics_qdir.setter
    def optics_qdir(self, optics_qdir):
        self.optics_qdir[...] = optics_qdir
    
    @property
    def optics_intraband(self):
        """
        Element optics_intraband ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 104
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__optics_intraband()
    
    @optics_intraband.setter
    def optics_intraband(self, optics_intraband):
        _OptaPyDOS.f90wrap_od_parameters__set__optics_intraband(optics_intraband)
    
    @property
    def optics_drude_broadening(self):
        """
        Element optics_drude_broadening ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 105
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__optics_drude_broadening()
    
    @optics_drude_broadening.setter
    def optics_drude_broadening(self, optics_drude_broadening):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__optics_drude_broadening(optics_drude_broadening)
    
    @property
    def optics_lossfn_gaussian(self):
        """
        Element optics_lossfn_gaussian ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 106
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__optics_lossfn_gaussian()
    
    @optics_lossfn_gaussian.setter
    def optics_lossfn_gaussian(self, optics_lossfn_gaussian):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__optics_lossfn_gaussian(optics_lossfn_gaussian)
    
    @property
    def optics_lossfn_broadening(self):
        """
        Element optics_lossfn_broadening ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 107
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__optics_lossfn_broadening()
    
    @optics_lossfn_broadening.setter
    def optics_lossfn_broadening(self, optics_lossfn_broadening):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__optics_lossfn_broadening(optics_lossfn_broadening)
    
    @property
    def core_geom(self):
        """
        Element core_geom ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 111
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__core_geom()
    
    @core_geom.setter
    def core_geom(self, core_geom):
        _OptaPyDOS.f90wrap_od_parameters__set__core_geom(core_geom)
    
    @property
    def core_qdir(self):
        """
        Element core_qdir ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 112
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_parameters__array__core_qdir(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            core_qdir = self._arrays[array_handle]
        else:
            core_qdir = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_parameters__array__core_qdir)
            self._arrays[array_handle] = core_qdir
        return core_qdir
    
    @core_qdir.setter
    def core_qdir(self, core_qdir):
        self.core_qdir[...] = core_qdir
    
    @property
    def core_lai_broadening(self):
        """
        Element core_lai_broadening ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 113
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__core_lai_broadening()
    
    @core_lai_broadening.setter
    def core_lai_broadening(self, core_lai_broadening):
        _OptaPyDOS.f90wrap_od_parameters__set__core_lai_broadening(core_lai_broadening)
    
    @property
    def core_type(self):
        """
        Element core_type ftype=character(len=20) pytype=str
        
        
        Defined at ../src/parameters.f90 line 114
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__core_type()
    
    @core_type.setter
    def core_type(self, core_type):
        _OptaPyDOS.f90wrap_od_parameters__set__core_type(core_type)
    
    @property
    def lai_gaussian_width(self):
        """
        Element lai_gaussian_width ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 115
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lai_gaussian_width()
    
    @lai_gaussian_width.setter
    def lai_gaussian_width(self, lai_gaussian_width):
        _OptaPyDOS.f90wrap_od_parameters__set__lai_gaussian_width(lai_gaussian_width)
    
    @property
    def lai_gaussian(self):
        """
        Element lai_gaussian ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 116
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lai_gaussian()
    
    @lai_gaussian.setter
    def lai_gaussian(self, lai_gaussian):
        _OptaPyDOS.f90wrap_od_parameters__set__lai_gaussian(lai_gaussian)
    
    @property
    def lai_lorentzian_width(self):
        """
        Element lai_lorentzian_width ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 117
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lai_lorentzian_width()
    
    @lai_lorentzian_width.setter
    def lai_lorentzian_width(self, lai_lorentzian_width):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__lai_lorentzian_width(lai_lorentzian_width)
    
    @property
    def lai_lorentzian_scale(self):
        """
        Element lai_lorentzian_scale ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 118
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lai_lorentzian_scale()
    
    @lai_lorentzian_scale.setter
    def lai_lorentzian_scale(self, lai_lorentzian_scale):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__lai_lorentzian_scale(lai_lorentzian_scale)
    
    @property
    def lai_lorentzian_offset(self):
        """
        Element lai_lorentzian_offset ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 119
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lai_lorentzian_offset()
    
    @lai_lorentzian_offset.setter
    def lai_lorentzian_offset(self, lai_lorentzian_offset):
        \
            _OptaPyDOS.f90wrap_od_parameters__set__lai_lorentzian_offset(lai_lorentzian_offset)
    
    @property
    def lai_lorentzian(self):
        """
        Element lai_lorentzian ftype=logical pytype=bool
        
        
        Defined at ../src/parameters.f90 line 120
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lai_lorentzian()
    
    @lai_lorentzian.setter
    def lai_lorentzian(self, lai_lorentzian):
        _OptaPyDOS.f90wrap_od_parameters__set__lai_lorentzian(lai_lorentzian)
    
    @property
    def lenconfac(self):
        """
        Element lenconfac ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/parameters.f90 line 122
        
        """
        return _OptaPyDOS.f90wrap_od_parameters__get__lenconfac()
    
    @lenconfac.setter
    def lenconfac(self, lenconfac):
        _OptaPyDOS.f90wrap_od_parameters__set__lenconfac(lenconfac)
    
    def __str__(self):
        ret = ['<od_parameters>{\n']
        ret.append('    output_format : ')
        ret.append(repr(self.output_format))
        ret.append(',\n    devel_flag : ')
        ret.append(repr(self.devel_flag))
        ret.append(',\n    iprint : ')
        ret.append(repr(self.iprint))
        ret.append(',\n    energy_unit : ')
        ret.append(repr(self.energy_unit))
        ret.append(',\n    length_unit : ')
        ret.append(repr(self.length_unit))
        ret.append(',\n    legacy_file_format : ')
        ret.append(repr(self.legacy_file_format))
        ret.append(',\n    kpoint_mp_grid : ')
        ret.append(repr(self.kpoint_mp_grid))
        ret.append(',\n    dos : ')
        ret.append(repr(self.dos))
        ret.append(',\n    compare_dos : ')
        ret.append(repr(self.compare_dos))
        ret.append(',\n    pdos : ')
        ret.append(repr(self.pdos))
        ret.append(',\n    jdos : ')
        ret.append(repr(self.jdos))
        ret.append(',\n    compare_jdos : ')
        ret.append(repr(self.compare_jdos))
        ret.append(',\n    optics : ')
        ret.append(repr(self.optics))
        ret.append(',\n    core : ')
        ret.append(repr(self.core))
        ret.append(',\n    fixed : ')
        ret.append(repr(self.fixed))
        ret.append(',\n    adaptive : ')
        ret.append(repr(self.adaptive))
        ret.append(',\n    linear : ')
        ret.append(repr(self.linear))
        ret.append(',\n    quad : ')
        ret.append(repr(self.quad))
        ret.append(',\n    compute_band_energy : ')
        ret.append(repr(self.compute_band_energy))
        ret.append(',\n    adaptive_smearing : ')
        ret.append(repr(self.adaptive_smearing))
        ret.append(',\n    fixed_smearing : ')
        ret.append(repr(self.fixed_smearing))
        ret.append(',\n    linear_smearing : ')
        ret.append(repr(self.linear_smearing))
        ret.append(',\n    dos_per_volume : ')
        ret.append(repr(self.dos_per_volume))
        ret.append(',\n    efermi_user : ')
        ret.append(repr(self.efermi_user))
        ret.append(',\n    efermi_choice : ')
        ret.append(repr(self.efermi_choice))
        ret.append(',\n    finite_bin_correction : ')
        ret.append(repr(self.finite_bin_correction))
        ret.append(',\n    hybrid_linear : ')
        ret.append(repr(self.hybrid_linear))
        ret.append(',\n    hybrid_linear_grad_tol : ')
        ret.append(repr(self.hybrid_linear_grad_tol))
        ret.append(',\n    numerical_intdos : ')
        ret.append(repr(self.numerical_intdos))
        ret.append(',\n    compute_band_gap : ')
        ret.append(repr(self.compute_band_gap))
        ret.append(',\n    set_efermi_zero : ')
        ret.append(repr(self.set_efermi_zero))
        ret.append(',\n    dos_min_energy : ')
        ret.append(repr(self.dos_min_energy))
        ret.append(',\n    dos_max_energy : ')
        ret.append(repr(self.dos_max_energy))
        ret.append(',\n    dos_spacing : ')
        ret.append(repr(self.dos_spacing))
        ret.append(',\n    dos_nbins : ')
        ret.append(repr(self.dos_nbins))
        ret.append(',\n    pdos_string : ')
        ret.append(repr(self.pdos_string))
        ret.append(',\n    jdos_max_energy : ')
        ret.append(repr(self.jdos_max_energy))
        ret.append(',\n    jdos_spacing : ')
        ret.append(repr(self.jdos_spacing))
        ret.append(',\n    scissor_op : ')
        ret.append(repr(self.scissor_op))
        ret.append(',\n    exclude_bands : ')
        ret.append(repr(self.exclude_bands))
        ret.append(',\n    num_exclude_bands : ')
        ret.append(repr(self.num_exclude_bands))
        ret.append(',\n    optics_geom : ')
        ret.append(repr(self.optics_geom))
        ret.append(',\n    optics_qdir : ')
        ret.append(repr(self.optics_qdir))
        ret.append(',\n    optics_intraband : ')
        ret.append(repr(self.optics_intraband))
        ret.append(',\n    optics_drude_broadening : ')
        ret.append(repr(self.optics_drude_broadening))
        ret.append(',\n    optics_lossfn_gaussian : ')
        ret.append(repr(self.optics_lossfn_gaussian))
        ret.append(',\n    optics_lossfn_broadening : ')
        ret.append(repr(self.optics_lossfn_broadening))
        ret.append(',\n    core_geom : ')
        ret.append(repr(self.core_geom))
        ret.append(',\n    core_qdir : ')
        ret.append(repr(self.core_qdir))
        ret.append(',\n    core_lai_broadening : ')
        ret.append(repr(self.core_lai_broadening))
        ret.append(',\n    core_type : ')
        ret.append(repr(self.core_type))
        ret.append(',\n    lai_gaussian_width : ')
        ret.append(repr(self.lai_gaussian_width))
        ret.append(',\n    lai_gaussian : ')
        ret.append(repr(self.lai_gaussian))
        ret.append(',\n    lai_lorentzian_width : ')
        ret.append(repr(self.lai_lorentzian_width))
        ret.append(',\n    lai_lorentzian_scale : ')
        ret.append(repr(self.lai_lorentzian_scale))
        ret.append(',\n    lai_lorentzian_offset : ')
        ret.append(repr(self.lai_lorentzian_offset))
        ret.append(',\n    lai_lorentzian : ')
        ret.append(repr(self.lai_lorentzian))
        ret.append(',\n    lenconfac : ')
        ret.append(repr(self.lenconfac))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_parameters = Od_Parameters()

class Od_Pdos(f90wrap.runtime.FortranModule):
    """
    Module od_pdos
    
    
    Defined at ../src/pdos.F90 lines 23-926
    
    """
    @staticmethod
    def pdos_calculate():
        """
        pdos_calculate()
        
        
        Defined at ../src/pdos.F90 lines 66-110
        
        
        """
        _OptaPyDOS.f90wrap_pdos_calculate()
    
    @property
    def dos_partial(self):
        """
        Element dos_partial ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/pdos.F90 line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_pdos__array__dos_partial(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dos_partial = self._arrays[array_handle]
        else:
            dos_partial = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_pdos__array__dos_partial)
            self._arrays[array_handle] = dos_partial
        return dos_partial
    
    @dos_partial.setter
    def dos_partial(self, dos_partial):
        self.dos_partial[...] = dos_partial
    
    @property
    def matrix_weights(self):
        """
        Element matrix_weights ftype=real(kind=dp) pytype=float
        
        
        Defined at ../src/pdos.F90 line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _OptaPyDOS.f90wrap_od_pdos__array__matrix_weights(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            matrix_weights = self._arrays[array_handle]
        else:
            matrix_weights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _OptaPyDOS.f90wrap_od_pdos__array__matrix_weights)
            self._arrays[array_handle] = matrix_weights
        return matrix_weights
    
    @matrix_weights.setter
    def matrix_weights(self, matrix_weights):
        self.matrix_weights[...] = matrix_weights
    
    def __str__(self):
        ret = ['<od_pdos>{\n']
        ret.append('    dos_partial : ')
        ret.append(repr(self.dos_partial))
        ret.append(',\n    matrix_weights : ')
        ret.append(repr(self.matrix_weights))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

od_pdos = Od_Pdos()

class Xmgrace_Utils(f90wrap.runtime.FortranModule):
    """
    Module xmgrace_utils
    
    
    Defined at ../src/xmgrace_utils.f90 lines 23-507
    
    """
    @staticmethod
    def xmgu_setup(unit):
        """
        xmgu_setup(unit)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 45-131
        
        Parameters
        ----------
        unit : int
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_setup(unit=unit)
    
    @staticmethod
    def xmgu_legend(unit):
        """
        xmgu_legend(unit)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 137-160
        
        Parameters
        ----------
        unit : int
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_legend(unit=unit)
    
    @staticmethod
    def xmgu_title(unit, min_x, max_x, min_y, max_y, title):
        """
        xmgu_title(unit, min_x, max_x, min_y, max_y, title)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 167-207
        
        Parameters
        ----------
        unit : int
        min_x : float
        max_x : float
        min_y : float
        max_y : float
        title : str
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_title(unit=unit, min_x=min_x, max_x=max_x, min_y=min_y, \
            max_y=max_y, title=title)
    
    @staticmethod
    def xmgu_subtitle(unit, subtitle):
        """
        xmgu_subtitle(unit, subtitle)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 213-231
        
        Parameters
        ----------
        unit : int
        subtitle : str
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_subtitle(unit=unit, subtitle=subtitle)
    
    @staticmethod
    def xmgu_axis(unit, axis, label):
        """
        xmgu_axis(unit, axis, label)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 237-342
        
        Parameters
        ----------
        unit : int
        axis : str
        label : str
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_axis(unit=unit, axis=axis, label=label)
    
    @staticmethod
    def xmgu_data(unit, field, x_data, y_data):
        """
        xmgu_data(unit, field, x_data, y_data)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 347-376
        
        Parameters
        ----------
        unit : int
        field : int
        x_data : float array
        y_data : float array
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_data(unit=unit, field=field, x_data=x_data, \
            y_data=y_data)
    
    @staticmethod
    def xmgu_data_header(unit, field, colour, legend):
        """
        xmgu_data_header(unit, field, colour, legend)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 379-469
        
        Parameters
        ----------
        unit : int
        field : int
        colour : int
        legend : str
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_data_header(unit=unit, field=field, colour=colour, \
            legend=legend)
    
    @staticmethod
    def xmgu_vertical_line(unit, x_coord, y_max, y_min):
        """
        xmgu_vertical_line(unit, x_coord, y_max, y_min)
        
        
        Defined at ../src/xmgrace_utils.f90 lines 472-507
        
        Parameters
        ----------
        unit : int
        x_coord : float
        y_max : float
        y_min : float
        
        ==================================================================
        """
        _OptaPyDOS.f90wrap_xmgu_vertical_line(unit=unit, x_coord=x_coord, y_max=y_max, \
            y_min=y_min)
    
    _dt_array_initialisers = []
    

xmgrace_utils = Xmgrace_Utils()

def my_dcopy(n, dx, incx, dy, incy):
    """
    my_dcopy(n, dx, incx, dy, incy)
    
    
    Defined at ../src/comms.F90 lines 667-730
    
    Parameters
    ----------
    n : int
    dx : float array
    incx : int
    dy : float array
    incy : int
    
    """
    _OptaPyDOS.f90wrap_my_dcopy(n=n, dx=dx, incx=incx, dy=dy, incy=incy)

def my_zcopy(n, zx, incx, zy, incy):
    """
    my_zcopy(n, zx, incx, zy, incy)
    
    
    Defined at ../src/comms.F90 lines 732-777
    
    Parameters
    ----------
    n : int
    zx : complex array
    incx : int
    zy : complex array
    incy : int
    
    """
    _OptaPyDOS.f90wrap_my_zcopy(n=n, zx=zx, incx=incx, zy=zy, incy=incy)

def my_icopy(n, zx, incx, zy, incy):
    """
    my_icopy(n, zx, incx, zy, incy)
    
    
    Defined at ../src/comms.F90 lines 779-821
    
    Parameters
    ----------
    n : int
    zx : int array
    incx : int
    zy : int array
    incy : int
    
    """
    _OptaPyDOS.f90wrap_my_icopy(n=n, zx=zx, incx=incx, zy=zy, incy=incy)

