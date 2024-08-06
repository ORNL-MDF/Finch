# Finch examples

Examples included in Finch are scan path creation and various versions of a single line additive case.


# Finch inputs

All inputs are passed to Finch in a JSON format.

## Temporal parameters (`time`)

- `start_time`: Simulation start time
  - units: `s`
- `end_time`:	Simulation end time
  - units: `s`
- `Co`:	Courant number
  - units: unitless
- `total_output_steps`: Frequency to output simulation data (across the full grid)
  - units: unitless
- `total_monitor_steps`: Frequency to output timing information
  - units: unitless

## Spatial parameters (`space`)

- `initial_temperature`: Initial temperature at all grid points
  - units: `K`
- `cell_size`: Grid cell size
  - units: `m`
- `global_low_corner`: Bottom corner of the physical domain
  - units: `m`
- `global_high_corner`: Top corner of the physical domain
  - units: `m`
- `ranks_per_dim`: MPI ranks per dimension (replaced if incompatible with resource set)
  - units: unitless
  - optional (defaults to MPI-determined cartesian domain decomposition)

## Material properties (`properties`)

- `density`: Material density
  - units: `kg/m^3`
- `specific_heat`: Specific heat capacity
  - units: `J/g/K`
- `thermal_conductivity`: Thermal conductivity
  - units: `W/m/K`
- `latent_heat`: Latent heat of fusion
  - units: `J/kg`
- `solidus`: Solidus temperature
  - units: `K`
- `liquidus`: Liquidus temperature
  - units: `K`

## Laser source parameters (`source`)
- `absorption`: Laser absorption
  - units: unitless
- `two_sigma`: Laser beam radius (half D4_sigma beam diameter)
  - units: `m`
- `scan_path_file`: File containing laser path information



## Output sampling (`sampling`)
This entire section is optional.

- `type`: Type of sampling
  - options: `solidification_data` (outputs sampled solidification data with spatial position x, y, z; melting time tm; solidification start time (time at which the location goes below the liquidus temperature) ts; solidification rate R; and (optionally) temperature gradients Gx, Gy, Gz)
- `format`: Output format
  - options: `default` (output sampled solidification data) and `exaca` (output only sampled solidification data relevant to ExaCA microstructure prediction: does not output Gx, Gy, Gz)
- `directory_name`: Path to save output
  - optional (defaults to "solidification/", within the current directory)


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
  - units: `degrees`
- `bi_direction`: If true, reverse the scan direction for every line
  - boolean
  - optional (defaults to true)
