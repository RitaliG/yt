import numpy as np
import os

from yt.utilities.io_handler import BaseIOHandler, BaseParticleIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py


class PlutoStaticIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "pluto_static"

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError

    def _read_fluid_selection(self, chunks, selector, fields, size):
        data = {}
        entry = int(os.path.basename(self.ds.parameter_filename).replace('.flt.h5','').replace('data.','').replace('.dbl.h5',''))
        
        for field in fields:
            data[field] = np.empty(size, dtype="float64")

        with h5py.File(self.ds.parameter_filename, "r") as fh:
            ind = 0
            for chunk in chunks:
                for grid in chunk.objs:
                    nd = 0
                    for field in fields:
                        ftype, fname = field
                        position = '/Timestep_%d/vars/%s'%(entry,fname)
                        values = np.transpose(fh[position][:].astype("=f8"), (2, 1, 0))
                        nd = grid.select(selector, values, data[field], ind)
                    ind += nd
        return data

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError
