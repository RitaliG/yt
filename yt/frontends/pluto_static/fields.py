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

    # In PlutoStatic, conservative variables are written out.

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system
        
        # Add tracer fields
        for i in range(1, self.ds.ntracers+1):
            if ("pluto_static", "tr%d"%i) in self.field_list:
                self.add_output_field(
                    ("pluto_static", "Tracer %d"%i), sampling_type="cell",
                    units="",
                )
                self.alias(
                    ("pluto_static", "Tracer %d"%i),
                    ("gas", "Tracer %d"%i),
                )
                self.alias(
                    ("pluto_static", "Tracer %d"%i),
                    ("gas", "tracer %d"%i),
                )
            
        '''
        # Add velocity fields
        for comp in "xyz":
            self.add_field(
                ("gas", f"velocity_{comp}"),
                sampling_type="cell",
                function=velocity_field(comp),
                units=unit_system["velocity"],
            )

        # Add pressure field
        if ("PlutoStatic", "GasEnergy") in self.field_list:
            self.add_output_field(
                ("PlutoStatic", "GasEnergy"), sampling_type="cell", units=pres_units
            )
            self.alias(
                ("gas", "thermal_energy"),
                ("PlutoStatic", "GasEnergy"),
                units=unit_system["pressure"],
            )

            def _pressure(field, data):
                return (data.ds.gamma - 1.0) * data["PlutoStatic", "GasEnergy"]

        else:

            def _pressure(field, data):
                return (data.ds.gamma - 1.0) * (
                    data["PlutoStatic", "Energy"] - data["gas", "kinetic_energy_density"]
                )

        self.add_field(
            ("gas", "pressure"),
            sampling_type="cell",
            function=_pressure,
            units=unit_system["pressure"],
        )

        def _specific_total_energy(field, data):
            return data["PlutoStatic", "Energy"] / data["PlutoStatic", "density"]

        self.add_field(
            ("gas", "specific_total_energy"),
            sampling_type="cell",
            function=_specific_total_energy,
            units=unit_system["specific_energy"],
        )

        # Add temperature field
        def _temperature(field, data):
            return (
                data.ds.mu
                * data["gas", "pressure"]
                / data["gas", "density"]
                * mh
                / kboltz
            )

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=unit_system["temperature"],
        )

        # Add color field if present (scalar0 / density)
        if ("PlutoStatic", "scalar0") in self.field_list:
            self.add_output_field(
                ("PlutoStatic", "scalar0"),
                sampling_type="cell",
                units=rho_units,
            )

            def _color(field, data):
                return data["PlutoStatic", "scalar0"] / data["PlutoStatic", "density"]

            self.add_field(
                ("PlutoStatic", "color"),
                sampling_type="cell",
                function=_color,
                units="",
            )

            self.alias(
                ("gas", "color"),
                ("PlutoStatic", "color"),
                units="",
            )

            # Using color field to define metallicity field, where a color of 1
            # indicates solar metallicity

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
