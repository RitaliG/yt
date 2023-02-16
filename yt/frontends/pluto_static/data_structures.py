import os
import weakref
from abc import ABC, abstractmethod
import numpy as np

from yt.data_objects.index_subobjects.stretched_grid import StretchedGrid
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import PlutoStaticFieldInfo
from .definitions import pluto_def_constants

class PlutoStaticGrid(StretchedGrid):
    _id_offset = 0

    def __init__(self, id, cell_widths, filename, index, level, dims):
        super().__init__(id=id, filename=filename, index=index, cell_widths=cell_widths)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dims


class PlutoStaticHierarchy(GridIndex, ABC):
    grid = PlutoStaticGrid

    def __init__(self, ds, dataset_type="pluto_static"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        with h5py.File(self.index_filename, mode="r") as h5f:
            self.field_list = [("pluto_static", "%s"%k) for k in h5f[list(h5f.keys())[0]+'/vars/']]
        
    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grid_particle_count[0][0] = 0
        self.grid_levels[0][0] = 0
        self.max_level = 0
        
    def _parse_grid_data(self, gridtxt):
        start = 10
        nx1 = int(gridtxt[start][:-1])
        start = start + (nx1+1)
        nx2 = int(gridtxt[start][:-1])
        start = start + (nx2+1)
        nx3 = int(gridtxt[start][:-1])
        
        start = 11
        cell_width1 = np.array([float(gridtxt[i].split()[-1])-float(gridtxt[i].split()[-2]) for i in range(start,start+nx1)])
        start = start + (nx1+1)
        cell_width2 = np.array([float(gridtxt[i].split()[-1])-float(gridtxt[i].split()[-2]) for i in range(start,start+nx2)])
        start = start + (nx2+1)
        cell_width3 = np.array([float(gridtxt[i].split()[-1])-float(gridtxt[i].split()[-2]) for i in range(start,start+nx3)])
        self._cell_widths = (cell_width1,cell_width2,cell_width3)

    def _populate_grid_objects(self):
        assert self.num_grids == 1
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            grid_file = os.path.join(os.path.dirname(self.index_filename), 'grid.out')        
            with open(grid_file, 'r') as gridtxt:
                txt = gridtxt.readlines()
                self._parse_grid_data(txt)
            g = self.grid(id=i, index=self, filename=self.index_filename, cell_widths=self._cell_widths,
                           level=self.grid_levels.flat[i], dims=self.grid_dimensions[i], )
            
            g._prepare_grid()
            g._setup_dx()
            self.grids[i] = g


class PlutoStaticDataset(Dataset, ABC):
    _index_class = PlutoStaticHierarchy
    _field_info_class = PlutoStaticFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="pluto_static",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
    ):
        self.fluid_types += ("pluto_static",)
        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the defaults, but if they are listed
        # in the HDF5 attributes for a file, which is loaded first, then those are
        # used instead.
        #
        length_unit =  None
        mass_unit = None
        time_unit = None
        velocity_unit = None
        density_unit  = None
        magnetic_unit = None
        
        if (os.path.exists(os.path.join(os.path.dirname(self.parameter_filename), 'definitions.h'))):
            def_file = os.path.join(os.path.dirname(self.parameter_filename), 'definitions.h')
            with open(def_file, 'r') as deftxt:
                for line in deftxt:
                    if 'UNIT_LENGTH' in line:
                        length_unit = line.split()[-1]
                    if 'UNIT_DENSITY' in line:
                        density_unit = line.split()[-1]
                    if 'UNIT_VELOCITY' in line:
                        velocity_unit = line.split()[-1]
            
            constant_names = list(pluto_def_constants.keys())
            try:
                length_unit = float(length_unit)
            except ValueError:
                for name in constant_names:
                    length_unit = length_unit.replace(name, 'pluto_def_constants[\"%s\"]'%name)
                length_unit = eval(length_unit)
            
            try:
                density_unit = float(density_unit)
            except ValueError:
                for name in constant_names:
                    density_unit = density_unit.replace(name, 'pluto_def_constants[\"%s\"]'%name)
                density_unit = eval(density_unit)
                
            try:
                velocity_unit = float(velocity_unit)
            except ValueError:
                for name in constant_names:
                    velocity_unit = velocity_unit.replace(name, 'pluto_def_constants[\"%s\"]'%name)
                velocity_unit = eval(velocity_unit)
                
            mass_unit = density_unit*length_unit**3
        
        if not length_unit:
            self.length_unit = self.quan(1.0, "AU")
        else:
            self.length_unit = self.quan(length_unit, "cm")
        if not mass_unit:
            mp = pluto_def_constants['CONST_mp']
            self.mass_unit = self.quan(self.quan(mp, "g/cm**3")*self.length_unit**3/self.quan(1.0,"Msun"),"Msun")
        else:
            self.mass_unit = self.quan(mass_unit, "g")
        if not velocity_unit:
            self.velocity_unit = self.quan(1.0, "km/s")
            self.time_unit = self.length_unit/self.velocity_unit
        else:
            self.velocity_unit = self.quan(velocity_unit, "cm/s")
            self.time_unit = self.length_unit/self.velocity_unit
        if not magnetic_unit:
            magnetic_unit = self.quan(1.0, "gauss")
        else:
            self.magnetic_unit = self.quan(magnetic_unit, "gauss")
            
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self): 
        data_file = self.parameter_filename
        grid_file = os.path.join(os.path.dirname(self.parameter_filename), 'grid.out')
        xmf_file  = self.parameter_filename[:-2]+'xmf'
        out_file  = os.path.join(os.path.dirname(self.parameter_filename), self.parameter_filename[-6:]+'.out')        
        with open(grid_file, 'r') as gridtxt:
            txt = gridtxt.readlines()
            self.domain_left_edge  = np.array([ float(txt[i].replace(',', '').replace('[','').replace(']','').split()[2]) for i in [6,7,8] ])
            self.domain_right_edge = np.array([ float(txt[i].replace(',', '').replace('[','').replace(']','').split()[3]) for i in [6,7,8] ])
            self.dimensionality    = int(txt[4][-2])
            self.domain_dimensions = np.array([ int(txt[i].replace(',', '').replace('[','').replace(']','').split()[4]) for i in [6,7,8] ] [:self.dimensionality])
            
            self.geometry = Geometry((txt[5].split()[-1]).lower())
            '''
            if (txt[5].split()[-1]=='CARTESIAN'):
                self.geometry = Geometry('CARTESIAN'.lower())#.CARTESIAN
            elif(txt[5].split()[-1]=='CYLINDRICAL'):
                self.geometry = Geometry.CYLINDRICAL
            elif(txt[5].split()[-1]=='POLAR'):
                self.geometry = Geometry.POLAR
            elif(txt[5].split()[-1]=='SPHERICAL'):
                self.geometry = Geometry.SPHERICAL
            '''
        with open(out_file, 'r') as outttxt:   
            txt = outttxt.readlines()
            entry = int(os.path.basename(self.parameter_filename).replace('.flt.h5','').replace('data.','').replace('.dbl.h5',''))
            self.current_time =  float(txt[entry].split()[1])
            self.ntracers = 0
            search = 'tr1'
            while(search in txt[entry].split()):
                self.ntracers += 1
                search = 'tr%d'%(self.ntracers+1)
            
        self.gamma = 5.0 / 3.0
        self.mu = 0.61          
        self.refine_by = 1 # no mesh refinement
        self._periodicity = (True, True, True)
        # PlutoStatic cannot yet be run as a cosmological simulation
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        test1 = ("dbl.h5" in filename) or ("flt.h5" in filename)
        directory = os.path.dirname(filename)
        test2 = os.path.exists(filename[:-2]+'xmf')
        test3 = os.path.exists(os.path.join(directory,'grid.out'))
        test4 = os.path.exists(os.path.join(directory, filename[-6:]+'.out'))
        try:
            fileh = h5py.File(filename, mode="r")
        except (ImportError, OSError):
            return False
        else:
            entries = list(fileh.keys())
            return test1 and test2 and test3 and test4 and "cell_coords" in entries and "node_coords" in entries
        finally:
            fileh.close()
