# Integrating cpp_oti_lib into Finch — spike record and level of effort

Branch `oti-spike`, commit `88a5f26`, on top of upstream `9f79a05`.

This document records what was actually changed to make Finch's heat solve
differentiable with respect to its material and source parameters using
`cpp_oti_lib` forward-mode AD, so that the effort can be judged from evidence
rather than estimate. Friction points are recorded whether or not they turned
out to matter.

**Scope of the spike:** temperature-field sensitivities w.r.t. six parameters,
first order, CPU (OpenMP) build, validated against finite differences on 1 and
4 MPI ranks. Out of scope and explicitly deferred: solidification-event (G, R)
sensitivities, beam parameters, GPU build, second order.

---

## 1. Result

One OTI solve produces d(QoI)/dp for all six parameters. Validated against
central finite differences from 12 ordinary double solves:

| parameter | d(T_sum)/dp — OTI | central FD | rel. diff |
|---|---:|---:|---:|
| density [kg/m³] | −1.204057e+03 | −1.204056e+03 | 2.0e−07 |
| specific_heat [J/kg/K] | −1.199533e+04 | −1.199532e+04 | 6.1e−07 |
| thermal_conductivity [W/m/K] | 1.136934e+04 | 1.136936e+04 | 1.6e−06 |
| latent_heat [J/kg] | −1.696343e−01 | −1.698053e−01 | 1.0e−03 |
| absorption [−] | 2.915397e+07 | 2.915391e+07 | 1.9e−06 |
| two_sigma [m] | −2.446787e+10 | −2.446805e+10 | 7.1e−06 |

The second QoI (single-node probe temperature) agrees to 1e−7…1e−9 on all six.
Field values from the OTI run match the double run to all printed digits
(`T_sum` 2.1919490395e+07 both).

**Cost** (533 steps, 51×41×21 nodes, 16 OpenMP threads, after warm-up):

| | time | vs one double solve |
|---|---:|---:|
| 1 double solve | 5.89 s | 1.00× |
| 1 OTI solve → all 6 derivatives | 12.17 s | **2.07×** |
| 12 double solves → same 6 derivatives by FD | 92.43 s | 15.70× |

The OTI jet carries 7 coefficients per node, so 2.07× is markedly sub-linear in
coefficient count at this problem size.

**Multi-rank:** on 4 ranks (2×2×1 decomposition) `T_sum` and every derivative
match the 1-rank values to all printed digits. **No communication code was
written or modified.**

---

## 2. Changes to Finch

Five existing headers modified, one added. Excluding comments and blank lines:

| file | +lines | −lines | what changed |
|---|---:|---:|---|
| `src/Finch_Solver.hpp` | ~85 | 40 | `Scalar` deduced from view; `MaterialProperties` struct; kernel locals retyped; `Math::exp`; return types pinned |
| `src/Finch_Grid.hpp` | ~35 | 12 | `Scalar` template param; array value type; BOV output path |
| `src/Finch_SolidificationData.hpp` | ~20 | 12 | `Scalar` template param; value projection into the event records |
| `src/Finch_Run.hpp` | ~6 | 5 | `Scalar` template param threaded through `Layer` |
| `src/Finch_Boundary.hpp` | ~6 | 3 | `Scalar` template param; boundary values retyped |
| `src/Finch_Scalar.hpp` | 110 (new) | — | math dispatch + `ScalarValue` trait |
| **total in `src/`** | **152** | **72** | |

Plus, outside the library: `applications/Finch_OTI.hpp` (162 lines, the entire
Finch↔cpp_oti_lib coupling) and `applications/Sensitivity.cpp` (306 lines,
the driver and validation harness — a demo, not infrastructure).

### 2.1 The design decision that kept this small

`Solver` needed **no new template parameter**. The field scalar type is
recovered from the view it already had:

```cpp
using scalar_type = typename ViewType::non_const_value_type;
```

Changing `Cabana::Grid::Array<double,…>` to `Array<Scalar,…>` in `Grid` was
therefore enough to retype the entire solver. Only `Grid`, `Boundary`, `Layer`
and `SolidificationData` gained a `typename Scalar = double` parameter.

### 2.2 The one abstraction added

Finch called `Kokkos::exp` fully qualified, which cannot be extended for a
user-defined type. `Finch_Scalar.hpp` replaces those calls with
`Finch::Math::exp`, which forwards to `Kokkos::exp` for arithmetic types (same
device path, same generated code) and otherwise resolves unqualified so ADL
finds the type's own overload.

Consequence: **Finch core has no dependency on cpp_oti_lib, or on any AD
library.** The coupling is one trait specialization, in the application:

```cpp
template <int M, int N, class Coeff>
struct ScalarValue<oti::otinum<M, N, Coeff>> {
    KOKKOS_INLINE_FUNCTION static double value( const oti::otinum<M,N,Coeff>& x )
    { return static_cast<double>( x.real() ); }
};
```

That is the whole interface. Everything else — arithmetic, comparisons, `exp`,
`fmin`/`fmax` — is found by ADL in namespace `oti` with no glue at all.

### 2.3 What deliberately stayed `double`

- **Mesh geometry.** `UniformMesh<double>`, cell size, node coordinates. These
  are not field values. `weight()` computes `dist_to_beam` in double and only
  the final accumulation becomes `Scalar` via `A_inv_`.
- **`solidus` / `liquidus`.** They only ever appear in comparisons that select a
  branch, never in arithmetic.
- **Solidification event records.** `View<double**>`, the hand-off format to
  ExaCA. Values are projected out with `Math::value`. This is the main deferred
  item — see §5.

---

## 3. Friction encountered

The honest tally. Two compile errors total, **neither in Finch's physics and
neither in the library interop**:

1. **`Finch::Properties` name collision.** `Finch_Inputs.hpp` already defines
   `struct Properties`. Renamed mine to `MaterialProperties`. Cost: ~2 minutes.
   Entirely my own naming choice.
2. **`otinum::data()` returns a const reference**, so `&local.data().data()`
   would not bind for an `MPI_Allreduce` output buffer. Used `&local[0]`
   instead. Cost: ~2 minutes. In my driver's MPI glue, not in Finch.

Everything else compiled first time, including the whole solver, the Cabana
array allocation, `ArrayOp::assign`, `Kokkos::deep_copy`, host mirrors, and the
halo exchange.

### 3.1 Things that needed thought but not debugging

- **`source()` return type.** The host overload returned `0.0` from two branches
  and an expression from a third under `auto`. Deduction would have failed once
  one branch became OTI. Fixed pre-emptively by pinning the return type to
  `scalar_type` and returning `scalar_type(0)`.
- **BOV output.** `BovWriter` hard-requires both `MpiTraits<value_type>` and
  `BovFormat<value_type>`, so a compound scalar cannot be written. `Grid::output`
  now projects the field onto a temporary double array (~20 lines under
  `if constexpr`) and writes that. Derivative components are not in the BOV
  file; the application reports them separately.
- **`fmin`/`fmax` argument types.** My dispatch templates take both arguments as
  the same `T`, so `fmax(m, 0.0)` with `m` of scalar type would not deduce. Call
  sites use `Scalar(0)` / `Scalar(1)`.

### 3.2 Environment setup — not AD effort, but real time

Listed separately because it inflates any naive stopwatch reading and has
nothing to do with the integration:

- Local Kokkos is **5.1.99**; Cabana refuses ≥4.6.99 unless built with
  `Kokkos_ENABLE_IMPL_VIEW_LEGACY=ON`. Finch's CI pins **Kokkos 4.1.00 +
  Cabana 0.6.1**, so that pair was built from source.
- Cabana 0.6.1 configure fails with `Cabana HDF5 support requires parallel
  HDF5` when it finds a serial HDF5; needs `-DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON`.

Roughly 40 minutes of build time, one-time.

### 3.3 Incidental pre-existing bug (not fixed — upstream's call)

`SolidificationData::enabled_` is a `bool` with no initializer, and the default
constructor `SolidificationData() {}` leaves it indeterminate. `Layer` only
constructs the sampler when `inputs.sampling.enabled`, then calls
`solidification_data_.update()` unconditionally, which reads `enabled_` and
returns early. With sampling disabled this reads an uninitialized bool. It
appears benign in practice but is undefined behaviour. The sensitivity driver
uses its own time loop and so does not depend on it.

---

## 4. The accuracy investigation

Worth recording because it is the part that would have been easy to get wrong,
and because the first result looked like a failure.

**First run (FD relative step 1e−6):** `thermal_conductivity`, `latent_heat` and
`two_sigma` agreed to ~1e−8. `density`, `specific_heat` and `absorption`
disagreed by **28%**.

The pattern was informative rather than alarming: `density` and `specific_heat`
enter the solver only through the product `rho*cp`, and their disagreements were
28.0% and 28.1% — so OTI was at least self-consistent.

Two hypotheses were considered and one was ruled out by inspection: a
near-zero-denominator artifact in the relative-difference metric. It does not
fit — the two values were −1.204e3 and −8.665e2, both order 10³, and the
normalized sensitivity `p·dQ/dp` was −9.0e6 against a QoI of 2.2e7, about 41% of
the total.

**Step-size sweep** (`FINCH_FD_STEP`), `d(T_sum)/d(density)`, OTI = −1.204057e+03:

| FD rel. step | central FD | rel. diff |
|---|---:|---:|
| 1e−8 | −1.204056e+03 | 2.0e−07 |
| 1e−7 | −1.204056e+03 | 5.6e−08 |
| 1e−6 | −8.665150e+02 | 2.8e−01 |
| 1e−5 | −1.564970e+03 | 2.3e−01 |
| 1e−4 | −1.171179e+03 | 2.7e−02 |
| 1e−3 | −1.208358e+03 | 3.6e−03 |

Accuracy *improving* as the step shrinks, with no roundoff floor, is backwards
for ordinary finite differences. That is the signature of a **discontinuity
staircase**: the finite difference is straddling branch flips, not resolving a
smooth slope.

**Confirmation.** Re-running with `latent_heat = 0` removes the mushy-zone
branch (`rho_Lf_by_dT` becomes 0, so both sides of the conditional are equal).
At the same 1e−6 step that previously failed, `density`, `specific_heat`,
`absorption` and `two_sigma` all agree to ~1e−8.

**Conclusion.** OTI was right throughout. The discrepancy was a finite-difference
artifact of the latent-heat branch

```cpp
( x >= solidus_ && x <= liquidus_ ) ? dt_/( rho_cp_ + rho_Lf_by_dT_ )
                                    : dt_/( rho_cp_ )
```

where the switch is large: `rho_cp` = 5.63e6 against `rho_Lf_by_dT` = 7.14e6, so
a single node crossing the threshold changes its effective timestep by ~2.3×.
With thousands of nodes near the liquidus over 533 steps, a perturbation of
1e−6 relative flips enough of them to swamp the smooth signal.

This is the non-smoothness flagged before the spike started, and it is worth
being precise about what it means:

- **The OTI derivative is the correct derivative of the smooth branch actually
  taken.** For the temperature-field QoIs used here it matched FD wherever FD
  was trustworthy.
- **It does not account for the branch boundary itself moving.** For a QoI that
  depends on *which* nodes are mushy, there is a second contribution the jet
  does not see. That matters much more for solidification-event quantities
  (§5) than for the temperature field.
- **Practical consequence: FD is the unreliable reference here, not AD.** Anyone
  validating this work must do a step-size sweep; a single-step FD check at a
  plausible-looking 1e−6 would have "disproved" a correct derivative.

Two further observations from the `latent_heat = 0` run:

- `d(T_sum)/d(thermal_conductivity)` becomes 1.9e−10 — structurally zero, because
  with adiabatic walls and no latent heat, conduction only redistributes heat and
  cannot change the temperature sum. OTI returns this exactly; the FD returns its
  noise floor, 5.4e−03. The reporting metric was given an absolute floor so such
  rows read `negligible` rather than a misleading 100% disagreement.
- `d(T_sum)/d(latent_heat)` at Lf = 0 is 1.96e−02 from OTI, while the FD is
  `nan` — the step `h = rel·|0|` is zero. AD gives a derivative at a parameter
  value where finite differencing structurally cannot.

---

## 5. What is not done

In rough order of value:

1. **Solidification-event sensitivities (G, R, cooling rate).** The highest-value
   item — these feed ExaCA. Events are currently `View<double**>` and values are
   projected out. Making them `Scalar`-valued has a useful side effect: the
   crossing interpolation `m = (temp − liquidus)/(temp − temp0)` becomes a jet,
   so the solidification *time* carries its own derivative. This is also where
   the moving-branch caveat above bites hardest, since event capture is itself a
   threshold test. **Estimate: 1–2 weeks including validation.**
2. **Beam parameters** (power, scan speed, spot position). Requires
   `MovingBeam`/`Segment` to become scalar-typed. **3–5 days.**
3. **GPU build.** cpp_oti_lib's Kokkos/CUDA path is already exercised by its own
   heat-equation study, and nothing added here is host-only, but it is untested
   in this combination. **2–4 days, mostly validation.**
4. **Second order** (`otinum<6,2>` → 28 coefficients) for the parameter Hessian.
   A one-line change to the driver's typedef; the cost question is the open one.
   **Days for the change, longer to characterize performance.**
5. **Sensitivities of `solidus`/`liquidus`** would require handling the branch
   predicate itself and is a genuinely different problem.
6. **Upstreaming.** Every change is `Scalar = double`-defaulted and the double
   path is bit-identical, so existing users see nothing. That is the argument to
   make to ORNL-MDF.

---

## 6. Revised level of effort

The pre-spike estimate was 2–4 days for a demo and 4–8 weeks for a validated
version. The demo took well under that. Revising:

| phase | estimate | confidence |
|---|---|---|
| ~~Spike: field sensitivities, CPU, validated~~ | **done** | — |
| Solidification-event sensitivities | 1–2 weeks | medium |
| Beam parameters | 3–5 days | high |
| GPU validation | 2–4 days | medium |
| Performance characterization, higher order | 1 week | low |
| Upstreamable PR (tests, docs, review cycles) | 1–2 weeks | medium |

**Total to a validated, upstreamable capability: 4–6 weeks**, revised down from
4–8. The reduction comes from the halo requiring no work at all and from the
solver templating being smaller than expected.

The characterization "plug and play" is fair for the **mechanical** integration:
152 substantive lines, two trivial compile errors, no communication code, and a
bit-identical double path. It is not fair for the **numerical** part — the 28%
discrepancy was a real investigation, and the branch-derivative question in §4
is a genuine open issue that will need a defensible answer before these
sensitivities are used for optimization or UQ.

---

## 7. Reproducing

```sh
# Dependencies (one-time, ~40 min)
cmake -S kokkos -B kokkos/build -DCMAKE_INSTALL_PREFIX=<prefix>/kokkos \
  -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON
cmake --build kokkos/build --parallel --target install

cmake -S cabana -B cabana/build -DCMAKE_INSTALL_PREFIX=<prefix>/cabana \
  -DCMAKE_PREFIX_PATH=<prefix>/kokkos -DCMAKE_BUILD_TYPE=Release \
  -DCabana_ENABLE_GRID=ON -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON
cmake --build cabana/build --parallel --target install

# Finch. cpp_oti_lib is auto-detected at ../cpp_oti_lib/include,
# or pass -DCPP_OTI_LIB_DIR=/path/to/cpp_oti_lib/include
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=<prefix>/cabana
cmake --build build --parallel

cd examples/single_line
FINCH_FD_STEP=1e-8 mpirun -np 1 ../../build/applications/finch_sensitivity \
  -i inputs_small.json
```

`FINCH_FD_STEP` controls the finite-difference relative step. Use it: the
default of 1e−6 sits inside the bad zone documented in §4, and is left there
deliberately so the failure is reproducible.
