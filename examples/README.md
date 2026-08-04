# Finch examples

The example sources stay separate from generated data. The run scripts stage
inputs under an ignored `output/` directory, so running a case does not add
BP5/BOV fields or CSV files beside the input files.

## Available examples

- `create_scan_path`: creates rotated hatch scan paths.
- `single_line/inputs_full_physics.json`: complete 3 mm line scan with
  temperature-dependent properties, independent mixed face conditions, a
  tabulated source with transient liquidus-depth feedback, solidification
  data collection, and scheduled solidus/liquidus melt-pool dimensions.

After installing Finch to the default `build/install` location, run:

```bash
examples/create_scan_path/run_example.sh
examples/single_line/run_example.sh full_physics 8
```

The second single-line argument is the number of MPI ranks. To use another
installation, set `FINCH_EXECUTABLE=/path/to/finch` or pass the executable as
argument three. `FINCH_MPI_ARGS` supplies launcher-specific options, for
example `FINCH_MPI_ARGS="--bind-to core --map-by core"`. Set
`FINCH_COMBINE_SOLIDIFICATION=1` to combine per-rank event files after a run.

All single-line inputs use `"ranks_per_dim": [0, 0, 0]`, allowing Finch to
select its geometry-aware Cartesian decomposition for the requested MPI size.

Finch does not require lock-test files. The default ADIOS2 backend does not
create them. They are temporary filesystem probes made by Open MPI's OMPIO
component only when the BOV fallback collectively opens its raw field; a
normal close removes them, and the example runner confines them to its output
directory. To prevent even those temporary probes, an Open MPI deployment can
select ROMIO with `FINCH_MPI_ARGS="--mca io romio321"`. Finch does not force
that MPI-specific choice because OMPIO versus ROMIO performance depends on the
target filesystem.


# Finch inputs

All inputs are passed to Finch in a JSON format.

## Temporal parameters (`time`)

- `start_time`: Simulation start time
  - units: `s`
- `end_time`:	Simulation end time
  - units: `s`
- `maximum_fourier_number`: Maximum explicit diffusion number
  - units: unitless
  - optional (defaults to 0.125; must not exceed 1/6)
- `total_monitor_steps`: Desired number of evenly distributed timing reports (zero disables progress reports)
  - units: unitless
  - optional (defaults to zero)

Finch calculates a conservative diffusion-stable timestep from the grid and
the complete material-property range. It then makes small timestep adjustments
to land exactly on true scan-path discontinuities, requested function times,
and `end_time`; the user can specify ordinary decimal times such as `0.00485`.
Collinear subdivisions with unchanged velocity and power do not split a
timestep interval. Scheduling is performed once on the host and adds no
per-step device reduction or MPI collective. Progress monitoring is based on
completed-step counts and does not introduce additional time stops. At each
progress report, Finch also integrates the discrete heat source with nodal
trapezoidal weights and reports commanded power, instantaneous absorptivity,
target absorbed power, integrated power, and relative error. Transient sources
additionally report measured and effective depths, principal and equivalent
D4σ diameters, profile azimuth, and absorption aspect ratio. This adds one
compact-support device reduction and one MPI reduction only at the requested
monitor reports; setting `total_monitor_steps` to zero disables both the
reports and this diagnostic.

## Spatial parameters (`space`)

- `initial_temperature`: Initial temperature at all grid points
  - units: `K`
- `cell_size`: Grid cell size
  - units: `m`
- `global_low_corner`: Bottom corner of the physical domain
  - units: `m`
- `global_high_corner`: Top corner of the physical domain
  - units: `m`
- `ranks_per_dim`: MPI ranks per dimension; `[0, 0, 0]` selects the
  geometry-aware automatic decomposition
  - units: unitless
  - optional (defaults to a geometry-aware Cartesian decomposition)

## Material properties (`properties`)

- `density`: Material density
  - units: `kg/m^3`
- `specific_heat`: Specific heat capacity
  - units: `J/kg/K`
  - accepts either a positive scalar or a temperature table
- `thermal_conductivity`: Thermal conductivity
  - units: `W/m/K`
  - accepts either a positive scalar or a temperature table
- `latent_heat`: Latent heat of fusion
  - units: `J/kg`
- `solidus`: Solidus temperature
  - units: `K`
- `liquidus`: Liquidus temperature
  - units: `K`

Temperature-dependent heat capacity and conductivity use the same table
format. Temperatures must be finite and strictly increasing, property values
must be finite and positive, and the two arrays must have the same length:

```json
"specific_heat": {
  "temperature": [300.0, 1000.0, 1600.0, 2500.0],
  "values": [750.0, 790.0, 830.0, 870.0]
}
```

Finch uses linear interpolation between input temperatures and clamps values
outside the supplied range. It resamples input curves once during
initialization into an internal uniform table, so device evaluation has
constant-time indexing and no binary search. Scalar properties, and tables
whose values are all identical, automatically use the dedicated
constant-property kernel.

## Laser source parameters (`source`)

- `type`: Heat-source model
  - options: `gaussian` (default) and `tabulated`
- `absorption`: Absorption model object
  - `type`: `constant` or `kelly`
- `two_sigma`: Gaussian beam radius (half D4_sigma beam diameter)
  - units: `m`
  - required for `gaussian`
- `scan_path_file`: File containing laser path information

The tabulated model replaces the Gaussian planar profile with bilinear
interpolation from a uniform two-dimensional table. Its source object is:

```json
"source": {
  "type": "tabulated",
  "absorption": {
    "type": "kelly",
    "geometry": "cone",
    "fresnel_absorptivity": 0.28,
    "conduction_absorptivity": 0.35,
    "transition_aspect_ratio": 1.0
  },
  "profile_file": "tabulated_profile.txt",
  "coordinate_frame": "scan_path",
  "minimum_depth": 60e-6,
  "axial_profile": {
    "exponent_slope": 0.0,
    "exponent_intercept": 1.0
  },
  "transient_depth": {
    "temperature": "liquidus"
  },
  "scan_path_file": "scan_path.txt"
}
```

The profile file starts with `nx ny`, followed by `x0 y0 dx dy`, then
`nx * ny` nonnegative values in row-major order (x varying fastest). Values
need not be normalized: Finch trapezoid-integrates the profile once and
normalizes the three-dimensional source to absorbed beam power. Exact zero
borders are cropped during initialization while retaining one interpolation
cell around positive values. The included `tabulated_profile.txt` samples a
circular Gaussian on a 10 µm grid; its bilinearly interpolated second moments
give a 100 µm D4σ diameter.

Absorption uses one canonical object. Constant absorption is

```json
"absorption": {
  "type": "constant",
  "coefficient": 0.35
}
```

The `kelly` model represents enhanced absorption from internal reflections in
a conical or cylindrical cavity. It requires tabulated-source liquidus-depth
feedback. Finch calculates the principal D4σ diameters and their azimuth from
the complete bilinear-profile second-moment tensor. The area-equivalent
diameter

`D4sigma_equivalent = sqrt(D4sigma_major * D4sigma_minor)`

is invariant to profile rotation and avoids exaggerating the aspect ratio for
an elliptical beam. Kelly absorption uses

`aspect_ratio = 2 * measured_liquidus_depth / D4sigma_equivalent`.

This area-preserving equivalent diameter is an axisymmetric closure for an
elliptical profile; it does not make the Kelly cavity model fully elliptical.
For a strongly elongated keyhole, a directional optical model would be a
separate higher-fidelity absorption model.

The measured depth is zero when no liquidus isotherm exists, selecting
`conduction_absorptivity`. Source geometry remains finite by separately using
`max(minimum_depth, measured_liquidus_depth)`. `transition_aspect_ratio`
defaults to 1.0. The Kelly calculation is performed once on the host after the
existing depth reduction; the device source kernel receives only the final
normalization scalar.

`coordinate_frame` is `global` or `scan_path`; the latter rotates the table
with the current scan direction. The reported profile azimuth is measured in
that selected frame; the principal diameters and equivalent diameter are
unchanged by azimuthal rotation. The one-sided axial profile below the beam is
`exp(-3 * (depth / current_depth)^p)`, where
`p = 2^(clamp(exponent_slope * log2(2 * current_depth /
D4sigma_equivalent) +
exponent_intercept, 0, 9))`. The common exponents 1, 2, 4, and 8 avoid a
general device `pow` evaluation.

`transient_depth` is optional. When present, Finch measures the selected
solidus or liquidus depth in a compact region under the source from the
completed previous temperature field, performs one MPI maximum reduction,
and uses the result on the next explicit source update. The lateral search is
restricted to the source support and the vertical search reaches the physical
domain bottom. When no selected isotherm is present, detected depth is zero
and the source uses `minimum_depth`. Omitting `transient_depth` uses
`minimum_depth` without an every-step traversal or collective.

## Boundary conditions (`boundary`)

The optional `boundary` object assigns conditions independently to `x_min`,
`x_max`, `y_min`, `y_max`, `z_min`, and `z_max`. Omitted faces, or the entire
omitted section, default to `adiabatic`.

```json
"boundary": {
  "z_min": {
    "type": "dirichlet",
    "value": 300.0
  },
  "z_max": {
    "type": "convection_radiation",
    "h": 10.0,
    "emissivity": 0.4,
    "ambient_temperature": 300.0
  }
}
```

Supported face types are:

- `adiabatic`: zero heat flux; no additional entries.
- `dirichlet`: fixed ghost-face temperature `value` in K.
- `neumann`: outward temperature `gradient` in K/m.
- `convection_radiation`: mixed convection and radiation with `h` in
  W/m2/K, `emissivity` in [0,1], and `ambient_temperature` in K.

The mixed condition applies convection and radiation as

`q = h (T - ambient_temperature) + emissivity sigma
(T^4 - ambient_temperature^4)`,

using the local temperature-dependent conductivity and a radiation
linearization about the completed previous boundary temperature.

The explicit three-dimensional diffusion stencil requires
`0 < maximum_fourier_number <= 1/6`.
All domain extents must be positive and evenly divisible by `cell_size`.



## Runtime functions (`functions`)

The optional `functions` object configures calculations and output without
placing strings, virtual calls, or schedule interpretation in device kernels.
Finch parses these entries at initialization and builds a concrete host-side
execution schedule. The supported function types are
`solidification_data`, `melt_pool_dimensions`, and `field_output`:

```json
"functions": {
  "solidification": {
    "type": "solidification_data",
    "execute": {"control": "every_step"},
    "write": {"control": "end"},
    "format": "default",
    "directory": "solidification"
  },
  "melt_pool": {
    "type": "melt_pool_dimensions",
    "execute": {"control": "output_count", "count": 20},
    "write": {"control": "execute"},
    "isotherms": ["solidus", "liquidus"],
    "coordinate_frame": "scan_path",
    "directory": "melt_pool_dimensions"
  },
  "temperature": {
    "type": "field_output",
    "execute": {"control": "output_count", "count": 2},
    "fields": ["temperature", "volumetric_heat_source"],
    "format": "adios2"
  }
}
```

`field_output.format` is optional. It defaults to `adios2` when Finch was
built with `Finch_ENABLE_ADIOS2=ON`, otherwise it defaults to `bov`. The
`adios2` backend keeps one `fields.bp` BP5 dataset open for the run, records
each execution as a time step, and embeds Fides and VTK image-data metadata.
Open `fields.bp` in ParaView with the ADIOS2 BP4/5 Fides reader. The `bov`
fallback writes one raw field and XDMF sidecar per execution.

`field_output.fields` accepts `temperature` and
`volumetric_heat_source` (instantaneous absorbed power density in W/m^3).
Derived source values are materialized only at output executions and only over
the compact source support. Temperature is written in simulation precision;
the derived visualization field uses single precision to reduce memory and I/O
without materially affecting source-integration diagnostics. Requesting the
source field does not affect normal timesteps, but it adds one dense field to
each BP5 step. Derived fields require the `adios2` format; the BOV fallback
supports temperature only.

The solidification function records spatial position, melting and
solidification times, cooling rate, and temperature gradients every step, then
writes per-rank CSV files at the end. `format: "exaca"` omits the gradient
columns.

The melt-pool function traverses owned nodes once for all selected isotherms,
includes linearly interpolated grid-edge crossings, and uses one fixed-size MPI
reduction per execution. `coordinate_frame` can be `global`, producing x/y/z
extents, or `scan_path`, producing length/width/depth along the current scan
direction. Rank zero writes `dimensions.csv`. An absent isotherm has
`active = 0` and zero dimensions. The field-output function writes the
temperature field at its scheduled execution times. Disabled or omitted
functions have no traversal, allocation, or MPI-collective cost.


# Scan path creation inputs

- `min_point`: Lower corner of scan path region
  - units: `m`
- `max_point`: Upper corner of scan path region
  - units: `m`
- `hatch`: Hatch spacing
  - units: `m`
- `angle`: Scan angle
  - units: `degrees`
- `power`: Laser power
  - units: `W`
- `speed`: Laser scan speed
  - units: `m/s`
- `dwell_time`: Dwell time
  - units: `s`
- `bi_direction`: If true, reverse the scan direction for every line
  - boolean
  - optional (defaults to true)
