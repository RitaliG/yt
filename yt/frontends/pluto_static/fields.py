import numpy as np
from unyt import Zsun

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import kboltz, mh

# Copied from Athena frontend
pres_units = "code_pressure"
erg_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / code_length**2 / code_time"
vel_units = "code_length/code_time"

'''
def velocity_field(comp):
    def _velocity(field, data):
        return data["PlutoStatic", f"momentum_{comp}"] / data["PlutoStatic", "density"]

    return _velocity
'''

class PlutoStaticFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("rho", (rho_units, ["density", "rho"], None)),
        ("vx1", (vel_units, ["vel_x"], None)),
        ("vx2", (vel_units, ["vel_y"], None)),
        ("vx3", (vel_units, ["vel_z"], None)),
        ("prs", (pres_units, ["prs", "pres", "pressure"], None)),
        # ("scalar0", (rho_units, [], None)),
    )

    known_particle_fields = ()
    
    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system
        
        # Add tracer fields
        for i in range(1, self.ds.ntracers+1):
            if ("pluto_static", "tr%d"%i) in self.field_list:
                self.add_output_field(
                    ("pluto_static", "tr%d"%i), sampling_type="cell",
                    units="",
                )
                self.alias(
                    ("gas", "Tracer %d"%i), 
                    ("pluto_static", "tr%d"%i), units="",
                )
                self.alias(
                    ("gas", "tracer %d"%i),  
                    ("pluto_static", "tr%d"%i), units="",
                )
        
        if ("pluto_static", "Temp") in self.field_list:
            self.add_output_field(
                    ("pluto_static", "Temp"), sampling_type="cell",
                    units=unit_system["temperature"],
                )
            self.alias(
                ("gas", "Temperature"), 
                ("pluto_static", "Temp"), units=unit_system["temperature"],
            )
            self.alias(
                ("gas", "temperature"), 
                ("pluto_static", "Temp"), units=unit_system["temperature"],
            )
        elif ("pluto_static", "temp") in self.field_list:
            self.add_output_field(
                    ("pluto_static", "temp"), sampling_type="cell",
                    units=unit_system["temperature"],
                )
            self.alias(
                ("gas", "Temperature"), 
                ("pluto_static", "temp"), units=unit_system["temperature"],
            )
            self.alias(
                ("gas", "temperature"), 
                ("pluto_static", "temp"), units=unit_system["temperature"],
            )
        elif ("pluto_static", "Temperature") in self.field_list:
            self.add_output_field(
                    ("pluto_static", "Temperature"), sampling_type="cell",
                    units=unit_system["temperature"],
                )
            self.alias(
                ("gas", "Temperature"), 
                ("pluto_static", "Temperature"), units=unit_system["temperature"],
            )
            self.alias(
                ("gas", "temperature"), 
                ("pluto_static", "Temperature"), units=unit_system["temperature"],
            )
        elif ("pluto_static", "temperature") in self.field_list:
            self.add_output_field(
                    ("pluto_static", "temperature"), sampling_type="cell",
                    units=unit_system["temperature"],
                )
            self.alias(
                ("gas", "Temperature"), 
                ("pluto_static", "temperature"), units=unit_system["temperature"],
            )
            self.alias(
                ("gas", "temperature"), 
                ("pluto_static", "temperature"), units=unit_system["temperature"],
            )
        else:
            def _temperature(field, data):
                return data.ds.mu * (mh / kboltz) * (data[("gas", "pressure")] / data[("gas", "density")])
                
            self.add_field(
                    ("gas", "Temperature"), sampling_type="cell",
                    function = _temperature,
                    units=unit_system["temperature"],
                )
            self.alias(
                ("gas", "temperature"), 
                ("gas", "Temperature"), units=unit_system["temperature"],
            )
            
        if ("pluto_static", "mach") in self.field_list:
            self.add_output_field(
                    ("pluto_static", "mach"), sampling_type="cell",
                    units="",
                )
            self.alias(
                ("gas", "Mach"), 
                ("pluto_static", "mach"), units="",
            )
            self.alias(
                ("gas", "mach"), 
                ("pluto_static", "mach"), units="",
            )
        else:
            def _mach(field, data):
                return np.sqrt(data[("pluto_static", "vx1")]**2 + data[("pluto_static", "vx2")]**2 + data[("pluto_static", "vx3")]**2)\
                    /np.sqrt(data.ds.gamma * (data[("gas", "pressure")] / data[("gas", "density")]))
                    
            self.add_field(
                    ("gas", "Mach"), sampling_type="cell",
                    function = _mach,
                    units="",
                )
            self.alias(
                ("gas", "mach"), 
                ("gas", "Mach"), units="",
            )
        
        if ("pluto_static", "ndens") in self.field_list:
            self.add_output_field(
                    ("pluto_static", "ndens"), sampling_type="cell",
                    units=unit_system["length"]**-3,
                )
            self.alias(
                ("gas", "ndens"), 
                ("pluto_static", "ndens"), units=unit_system["length"]**-3,
            )
            self.alias(
                ("gas", "Number Density"), 
                ("pluto_static", "ndens"), units=unit_system["length"]**-3,
            )
            self.alias(
                ("gas", "number density"), 
                ("pluto_static", "ndens"), units=unit_system["length"]**-3,
            )
            self.alias(
                ("gas", "Number density"), 
                ("pluto_static", "ndens"), units=unit_system["length"]**-3,
            )
        else:
           def _ndens(field, data):
                return data["gas", "density"]/(data.ds.mu * mh)
                    
           self.add_field(
                    ("gas", "Number Density"), sampling_type="cell",
                    function = _ndens,
                    units=unit_system["length"]**-3,
           )
           self.alias(
                ("gas", "number density"), 
                ("gas", "Number Density"), units=unit_system["length"]**-3,
           )
           self.alias(
                ("gas", "Number density"), 
                ("gas", "Number Density"), units=unit_system["length"]**-3,
           ) 
        
        def _velMag(field, data):
            print('PLUTO ', data)
            return np.sqrt(data[("pluto_static", "vx1")]**2 + data[("pluto_static", "vx2")]**2 + data[("pluto_static", "vx3")]**2)
            
        self.add_field(
                ("gas", "Speed"), sampling_type="cell",
                function = _velMag,
                units=unit_system["velocity"],
            )       
        self.alias(
            ("gas", "speed"), 
            ("gas", "Speed"), units=unit_system["velocity"],
        )
        self.alias(
            ("gas", "Velocity Magnitude"), 
            ("gas", "Speed"), units=unit_system["velocity"],
        )
        
            # Using color field to define metallicity field, where a color of 1
            # indicates solar metallicity
        '''
            def _metallicity(field, data):
                # Ensuring that there are no negative metallicities
                return np.clip(data[("PlutoStatic", "color")], 0, np.inf) * Zsun

            self.add_field(
                ("PlutoStatic", "metallicity"),
                sampling_type="cell",
                function=_metallicity,
                units="Zsun",
            )

            self.alias(
                ("gas", "metallicity"),
                ("PlutoStatic", "metallicity"),
                units="Zsun",
            )
        '''
    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
