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
- `total_output_steps`: Desired number of evenly distributed field outputs (zero disables field output)
  - units: unitless
- `total_monitor_steps`: Desired number of evenly distributed timing reports (zero disables progress reports)
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
- `ranks_per_dim`: MPI ranks per dimension (replaced by a geometry-aware decomposition if incompatible with the communicator size)
  - units: unitless
  - optional (defaults to a geometry-aware Cartesian decomposition)

## Material properties (`properties`)

- `density`: Material density
  - units: `kg/m^3`
- `specific_heat`: Specific heat capacity
  - units: `J/kg/K`
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
  - must be between 0 and 1
- `two_sigma`: Laser beam radius (half D4_sigma beam diameter)
  - units: `m`
- `scan_path_file`: File containing laser path information

The explicit three-dimensional diffusion stencil requires `0 < Co <= 1/6`.
All domain extents must be positive and evenly divisible by `cell_size`.



## Output sampling (`sampling`)
This entire section is optional.

- `type`: Type of sampling
  - options: `solidification_data` (outputs sampled solidification data with spatial position x, y, z; melting time tm; solidification time ts; cooling rate R in K/s; and, optionally, temperature gradients Gx, Gy, Gz in K/m)
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
  - units: `s`
- `bi_direction`: If true, reverse the scan direction for every line
  - boolean
  - optional (defaults to true)
