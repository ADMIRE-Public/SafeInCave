# YAML case definitions

A simulation case can be described in a YAML file instead of a Python script.
The YAML file is converted into a standalone Python script that uses the
public SafeInCave API, and that script is what runs:

```bash
sic run case.yaml          # convert and run
sic y2p case.yaml          # write the equivalent case.py and stop
```

Python case scripts keep working exactly as before — `sic run main.py` is
unchanged, and the two workflows can live side by side in the same folder.

## The schema is the API

There is no separate list of allowed keywords to keep in sync with the code.
Every `type:` value is a real SafeInCave class name and every key under it is
a real constructor parameter, validated against `inspect.signature` when the
file is read. Adding a new constitutive model file makes it usable in YAML
immediately, under its own parameter names.

Mistakes are reported with the location, the legal alternatives, the closest
match and the real signature, so the vocabulary is discoverable from the
errors themselves:

```
material.non_elastic[0]: unknown type 'DislocationCrep'.
Did you mean: DislocationCreep?
Legal types for 'material.non_elastic': DislocationCreep, LinearDashpot, ...

equations.momentum: invalid parameters for LinearMomentum.
Unknown parameter: 'thta'. Did you mean 'theta'?
Signature: LinearMomentum(grid, theta, solver_name='gmres', ...)
```

Objects that the transpiler wires up itself — `grid`, `eq_mom`, `eq_heat`,
`t_control`, `outputs`, `caverns`, `simulation_logger` — are never written in
YAML; they are connected automatically from the sections that define them.

## Values

All values are plain numbers in **SI units**: stresses in Pa, times in
seconds, temperatures in K, angles as the underlying class expects. There is
no unit syntax and no expressions, so `102*GPa` must be written `1.02e11` and
a boundary condition ending at 24 hours uses `time_values: [0.0, 86400.0]`.

Material parameters are per-element fields internally. A scalar is broadcast
to all elements; a mapping assigns values per mesh region (the physical volume
names of the Gmsh file):

```yaml
material:
  density: 2000.0
  elastic:
    - type: Spring
      E: {Salt: 1.02e11, Overburden: 1.8e11}   # per region
      nu: 0.3                                  # same everywhere
```

An unknown region name is reported when the case runs, together with the
region names the mesh actually defines.

## Structure of a case file

```yaml
grid:
  type: GridHandlerGMSH
  geometry_name: geom
  grid_folder: ../../grids/cube

equations:                       # 'momentum', 'heat', or both
  momentum: {type: LinearMomentum, theta: 0.5}

material:
  density: 2000.0                # plus specific_heat_capacity,
                                 # thermal_conductivity, thermal_expansion
  elastic:                       # -> Material.add_to_elastic
    - {type: Spring, E: 1.02e11, nu: 0.3, name: spring}
  non_elastic:                   # -> Material.add_to_non_elastic, in order
    - {type: Viscoelastic, eta: 1.05e13, E: 1.0e10, nu: 0.32, name: kelvin}
    - {type: DislocationCreep, A: 1.9e-20, Q: 51600, n: 3.0, name: creep}
  thermoelastic:                 # -> Material.add_to_thermoelastic
    - {type: Thermoelastic, alpha: 4.4e-5}

body_force: [0.0, 0.0, -9.81]    # gravity vector, m/s^2
initial_temperature: 293         # K, uniform

stages:
  - name: loading
    time:
      type: TimeController
      dt: 0.5
      initial_time: 0.0
      final_time: 24
      time_unit: hour
    bcs:
      momentum:
        - {type: DirichletBC, boundary_name: WEST, component: 0,
           values: [0.0, 0.0], time_values: [0.0, 86400.0]}
        - {type: NeumannBC, boundary_name: TOP, direction: 2, density: 0.0,
           ref_pos: 0.0, g: 0.0,
           values: [4.1e6, 1.6e7], time_values: [0.0, 7200.0]}
    outputs:
      - equation: momentum
        output_format: xdmf      # or vtkhdf
        folder: output/case_0
        fields:
          u: "Displacement (m)"
          eps_tot: "Total strain (-)"
    simulator:
      type: Simulator_M
      compute_elastic_response: true
```

The order of `non_elastic` entries is preserved, because the mechanisms are
composed in the order they are added.

`values` and `time_values` are the piecewise-linear tables the boundary
condition classes already use: the value is interpolated between the listed
times.

### Stages

`stages` is a list, and each entry runs a full simulation with its own time
controller, boundary conditions, outputs and simulator, while the grid,
equations and material are built once and carried over. This is how a
geostatic equilibrium step is followed by an operational step, using
`GeostaticStep` to check/equilibrate/commit the in-situ stress before the
operational stage continues from it:

```yaml
stages:
  - name: equilibrium
    time: {type: TimeControllerParabolic, n_time_steps: 20, initial_time: 0.0,
           final_time: 100, time_unit: day}
    bcs:
      momentum: [...]           # same tractions/rollers the operational stage uses
    outputs:
      - equation: momentum
        folder: output/equilibrium
        fields: {sig: "Stress (Pa)"}
    simulator:
      type: GeostaticStep
      tolerance: 1.0e-8
  - name: operation
    time: {type: TimeControllerAdaptive, initial_time: 0.0, initial_dt: 0.01,
           final_time: 30, time_unit: day, dt_max: 1.0}
    ...
    simulator:
      type: Simulator_M
      compute_elastic_response: false   # continue from the committed geostatic state
```

`GeostaticStep` checks the initial stress set on the momentum equation
against the applied loads/BCs, solves to equilibrium (always with
`compute_elastic_response=false` internally, so two documented defects in
the elastic pre-solve path are never reached), then commits the equilibrated
stress as the new reference and zeros displacement, strain and inelastic
history — hence `operation` should run with `compute_elastic_response:
false` to continue from that committed state. Note there is currently no
YAML key to set an initial stress before stage 0 (`apply_initial_stress` is
Python-API only), so a fully YAML-only geostatic workflow needs that set up
separately, e.g. via `sic y2p` and a small edit to the generated script.

### Heat and thermomechanics

Declaring a `heat` equation enables heat boundary conditions (including
`RobinBC`) and the thermal material properties. `Simulator_TM` requires both
`momentum` and `heat` to be defined; if one is missing the case is rejected
before anything runs.

```yaml
equations:
  momentum: {type: LinearMomentum, theta: 0.5}
  heat: {type: HeatDiffusion}
...
    bcs:
      momentum: [...]
      heat:
        - {type: RobinBC, boundary_name: Cavern, values: [303, 303], h: 5.0,
           time_values: [0.0, 86400.0]}
    simulator: {type: Simulator_TM}
```

### Caverns and monitoring

A stage may add caverns and a monitoring point; both are wired into the
simulator automatically:

```yaml
    caverns:
      - {type: Cavern_T, cavern_name: Cavern, T_values: [303, 303],
         time_values: [0.0, 86400.0]}
    logging:
      type: SimulationLogging
      target_point: [0.0, 0.0, 1.0]
      variables_to_track: [sxx, syy, szz]
```

## Splitting a case across files

Any part of a case can live in its own file, referenced with `!include`:

```yaml
# case.yaml
grid:
  type: GridHandlerGMSH
  geometry_name: geom
  grid_folder: ../../grids/cube

equations:
  momentum: {type: LinearMomentum, theta: 0.5}

material: !include materials/salt.yaml

body_force: [0.0, 0.0, -9.81]
initial_temperature: 293

stages:
  - !include stages/equilibrium.yaml
  - !include stages/operation.yaml
```

The included file holds the content of the node it replaces, so
`materials/salt.yaml` is the material mapping itself — not a file with a
`material:` key inside it:

```yaml
# materials/salt.yaml
density: 2000.0
elastic:
  - {type: Spring, E: 1.02e11, nu: 0.3}
non_elastic:
  - {type: DislocationCreep, A: 1.9e-20, Q: 51600, n: 3.0}
```

Paths are relative to the file that contains the tag, so an included file can
itself include others relative to its own location. Each tag takes exactly one
path; to combine several files, tag each entry, as the `stages:` list above
does. A file that ends up including itself, directly or through a chain, is
reported as a circular include rather than looping.

This is a good way to share one material definition or one loading schedule
between several cases. The generated script lists the files it was built from
in its header, and `sic y2p` still produces a single self-contained
script.

## Limitations, and the way around them

A YAML case covers what the standard API can express declaratively. It cannot
express custom Python: equation subclasses that publish extra fields (as
`examples/mechanics/1_triaxial/main.py` does with `run_after_solve`),
spatially varying initial temperature, or values computed from the geometry
file.

For those cases, export the script and continue in Python:

```bash
sic y2p case.yaml -o main.py
```

The exported script is exactly what `sic run case.yaml` executes, so nothing
changes behaviourally when you switch. Custom classes installed through the
[extensions mechanism](../implementation/intro.md) are picked up by the
registry automatically and can be used as `type:` values like any built-in
class.

Two YAML details worth knowing: values such as `1.9e-20` (no decimal point)
are read as strings by YAML 1.1 and converted to numbers automatically, and
unquoted `yes`/`no`/`on`/`off` are booleans — quote them if you mean strings.
