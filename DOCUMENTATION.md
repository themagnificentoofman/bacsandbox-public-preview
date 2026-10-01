# BacSandbox: Formal Documentation

## 1. Introduction

BacSandbox is a micro-scale bacterial environment simulator built with Python and Pygame. It models independent bacterial species in a 2D petri-dish-like environment where cells move, sense organic nutrient gradients, consume resources, compete, reproduce, age, and die under configurable environmental pressures.

The simulation is real-time: biological rates are expressed per second and multiplied by elapsed `dt`. The `FPS` constant controls visual smoothness, not the amount of biology performed per update.

Food nodes represent an abstract organic nutrient pool measured in pg glucose-equivalent. This lets respiration, fermentation, uptake, and biomass synthesis share one energy/carbon accounting unit without implying that all environmental matter is literal glucose.

## 2. Setup And Usage

The project standard is CPython 3.14, with development and CI pinned to **3.14.7**
in `.python-version`. Package metadata requires `>=3.14,<3.15`. Use the standard
GIL-enabled interpreter; free-threaded `3.14t` is not part of this validation target.
Later 3.14 patch upgrades should update `.python-version` and pass the same QA gates.

Follow the [README setup commands](README.md#setup) to create and activate `.venv`,
then choose a requirements-only installation or an editable native installation.
`requirements.txt` is the single runtime dependency list used by both routes:

| Dependency | Range | Role |
| --- | --- | --- |
| pygame-ce | `~=2.5.8` | Pygame UI; imported as `pygame` |
| NumPy | `>=2.3.4,<3` | Numerical buffers and reference kernels |
| pandas | `>=2.3.3,<4` | CSV analysis |
| Matplotlib | `>=3.10.7,<4` | Plotting |

`pyproject.toml` separately declares isolated native build requirements:
`setuptools>=80.9`, `wheel`, `pybind11>=3.0.1,<4`, and `numpy>=2.3.4,<3`.
The optional extension still requires a C++17 compiler. OpenMP is selected by
`BACSANDBOX_OPENMP=auto|on|off`; `off` disables OpenMP but still compiles C++.
Both install routes can now be used in the same environment without installing
competing Pygame distributions. Do not install `pygame` alongside `pygame-ce`.

The app records population and event history in `simulation_stats.csv` for
spreadsheet or external analysis workflows. For isolated UI runs after an editable
install, run from a separate directory containing copies of `species_data.csv` and
`env_config.json` to preserve tracked inputs and outputs.

### 2.1 Python 3.14 Migration And Verification

An existing Python 3.12 virtual environment cannot be upgraded in place. Before
switching, save its `python -m pip freeze` output outside the checkout and preserve
it as a backup. Deactivate it, move the old `.venv` out of the checkout, then create
a fresh `.venv` using the README's Python 3.14 commands. A moved virtual environment
is only a backup: restore it to its original path before using it for rollback.
Rebuild `bacteria_c_kernel` for the new interpreter; a CPython 3.12 binary is not a
CPython 3.14 build. Keep the backup until the new environment passes verification.

Run from the repository root in the new environment:

```bash
python --version
python -m pip check
python -c "import sys, pygame; print(sys.executable); print(pygame.version.ver); assert pygame.IS_CE"
python -m unittest discover
```

The version should match `.python-version` and the executable should be inside
`.venv`. For a native installation, also run:

```bash
python -c "import bacteria_c_kernel as k; print('OpenMP:', k.openmp_enabled()); assert hasattr(k, 'MetabolicKernel')"
```

`OpenMP: True` is required for an explicitly requested `on` build; `False` is
expected for `off`. Native tests skip intentionally in a requirements-only
installation. Scale tests and timing benchmarks have separate opt-in flags; a
successful suite with skips does not mean those checks ran. Tests that inject
synthetic FBA failures print expected rollback tracebacks; use the final unittest
summary to distinguish these from failures.

VS Code defaults to `.venv` and unittest discovery. Reselect **Python: Select
Interpreter** after migration if the workspace already has an interpreter choice.
The C++ editor configuration uses C++17 and leaves compiler discovery to the host;
any additional Python/pybind11 header paths must come from the active installation,
not another developer's absolute paths. Editor settings do not control the build.
Record the actual Python version, installed package versions, OpenMP state, and Git
commit with each research run. Passing this migration's tests does not establish
bit-for-bit equivalence with previous Python 3.12 stochastic outputs.

## 3. Core Units

* Distance: 1 world unit equals approximately $1\thinspace\mu\mathrm{m}$.
* Time: all simulation rates use elapsed seconds.
* Nutrient amount: 1 food unit equals $1\thinspace\mathrm{pg}$ glucose-equivalent organic nutrient.
* Cell size: `base_size` is a structural cell size in micrometers. The rendered `size` only varies from 90% to 100% with fullness, so stored food does not inflate storage capacity.
* ECM concentration: quorum-sensing molecules, bacteriocins, acetate, and phage are stored on a coarse spatial grid and updated every simulation tick.
* Screen scale: `zoom` is pixels per micrometer in the Grid tab.

The Grid scale indicator uses the same transform as the renderer:

$$
x_{\text{screen}} = x_{\text{world}} z
$$

Therefore a bar of length $p$ pixels represents:

$$
d_{\mu\mathrm{m}} = \frac{p}{z}
$$

The current UI displays this as `100 px = ...` next to the minimap.

## 4. Mathematical And Biological Models

### 4.1 Nutrient Storage And Cell Size

Nutrient capacity scales with structural volume, not temporary fullness<sup>[33](#ref-33)</sup>:

$$
C_{\text{food}} = k_{\text{cap}} \cdot s_{\text{base}}^3
$$

The rendered cell radius uses fullness only as a small starvation cue<sup>[33](#ref-33)</sup>:

$$
s_{\text{rendered}} = s_{\text{base}} \left(0.9 + 0.1 \cdot \frac{F}{C_{\text{food}}}\right)
$$

This avoids the unrealistic behavior where eating more food directly increases the maximum storage cap.

### 4.2 Metabolism And Uptake

Basal metabolism follows an allometric Kleiber-style relationship<sup>[17](#ref-17)</sup>. Because cell mass scales with volume, a mass exponent of 3/4 becomes a size exponent of 2.25:

$$
B = \beta \cdot s_{\text{base}}^{2.25}
$$

Nutrient uptake is surface-area limited<sup>[33](#ref-33)</sup>:

$$
U_{\max} = u \cdot s_{\text{base}}^2
$$

During each elapsed-time update:

$$
F_{t+\Delta t} = F_t - B \cdot \Delta t \cdot A(t)
$$

where $A(t)$ is the senescence maintenance factor.

### 4.3 Food Sensing And Chemotaxis

Bacteria do not teleport or directly seek a target. They sense local nutrient fields and bias their swimming direction and run duration. Each food node contributes a softened Euclidean concentration proxy<sup>[32](#ref-32)</sup>:

$$
C_i = \frac{A_i}{\pi(d_i^2 + a^2)}
$$

where $A_i$ is node amount, $d_i$ is Euclidean distance, and $a = 5.0\thinspace\mu\mathrm{m}$ prevents a zero-distance singularity.

Chemoreceptor response is log-compressed to approximate saturation<sup>[15](#ref-15)</sup>:

$$
S_i = g_{\text{eff}} \log\left(1 + \frac{C_i}{C_{\text{ref}}}\right)
$$

where $C_{\text{ref}} = 10^{-3}\thinspace\mathrm{pg}/\mu\mathrm{m}^2$ and $g_{\text{eff}}$ increases as the bacterium becomes hungry.

The local gradient direction is a weighted sum of node direction vectors<sup>[32](#ref-32)</sup>:

$$
\vec{G} = \sum_i S_i \frac{\vec{x}_i - \vec{x}_{\text{cell}}}{\lVert \vec{x}_i - \vec{x}_{\text{cell}} \rVert}
$$

### 4.4 Run-And-Tumble Motion

BacSandbox uses an E. coli-inspired run-and-tumble model<sup>[3](#ref-3)</sup>. Runs are sampled around 1 second and clamped to biological bounds:

$$
t_{\text{run}} \sim \mathrm{clamp}(\mathcal{N}(1.0, 0.2^2), 0.25, 3.0)
$$

Tumbles last about 0.1 seconds:

$$
t_{\text{tumble}} = 0.1\thinspace\mathrm{s}
$$

When the nutrient score improves, the next run is extended:

$$
t'_{\text{run}} = 2.5 \cdot t_{\text{run}}
$$

During a run, heading gently rotates toward the local gradient while movement remains forward-only:

$$
\Delta \theta_{\text{run}} = \mathrm{clamp}(\theta_G - \theta, -\omega g_{\text{eff}} \Delta t, \omega g_{\text{eff}} \Delta t)
$$

During a tumble, the cell performs one event-scale reorientation rather than a tiny frame-scaled wobble. The unbiased magnitude is centered near the classic observed E. coli reorientation angle of 68 degrees<sup>[3](#ref-3)</sup>:

$$
\lvert \Delta \theta_{\text{tumble}} \rvert \sim \mathcal{N}(1.19, 0.45^2)\thinspace\mathrm{rad}
$$

Translation during tumbles is multiplied by 0.05, representing motor activity being spent on reorientation rather than propulsion. Speed is never negative; apparent reversal is a new heading after a tumble, not backward swimming.

### 4.5 Rotational Diffusion

Brownian rotational noise follows the low-Reynolds scaling for small particles in fluid<sup>[4](#ref-4)</sup>:

$$
D_{\text{rot}} \propto \frac{1}{s_{\text{base}}^3}
$$

The code uses a reference value of about $0.16\thinspace\mathrm{rad}^2/\mathrm{s}$ for a $1\thinspace\mu\mathrm{m}$ sphere in water.

Because $D_{\text{rot}}$ is an angular variance rate, the sampled angular standard deviation is:

$$
\sigma_{\theta} = \sqrt{2D_{\text{rot}}\Delta t}
$$

This keeps Brownian heading noise tied to elapsed seconds instead of tick count.

### 4.6 Boundary Reflection

The petri dish is circular. If a bacterium reaches the dish boundary, it is moved slightly inside the rim and its velocity heading reflects off the outward normal only when it was moving outward<sup>[38](#ref-38)</sup>:

$$
\vec{v}' = \vec{v} - 2(\vec{v} \cdot \hat{n})\hat{n}
$$

This is applied only when:

$$
\vec{v} \cdot \hat{n} \gt 0
$$

### 4.7 Reproduction

Population-level reproduction probability follows Verhulst logistic pressure through carrying capacity<sup>[29](#ref-29)</sup>:

$$
\frac{dN}{dt} = rN\left(1 - \frac{N}{K}\right)
$$

At the individual level, binary fission requires stored nutrient above the reproduction threshold and a synthesis cost inspired by Stouthamer-style ATP/biomass costs<sup>[28](#ref-28)</sup>:

$$
F_{\text{cost}} = 40.0 \cdot s_{\text{base}}^3
$$

The parent divides only if both parent and child retain enough nutrient to survive the next update.

The implementation exposes this conversion as `FISSION_COST_PG_PER_UM3 = 40.0` so tests and model code share the same biological assumption.

### 4.8 Death, Predation, And Necromass Recycling

Dead bacteria do not become pure glucose. Death releases a lossy glucose-equivalent nutrient value from cytosolic leakage plus labile dry necromass<sup>[34](#ref-34)</sup>:

$$
F_{\text{recovered}} = \max(0, F_{\text{stored}}) f_{\text{cytosol}} + V_{\text{cell}} m_{\text{wet}} f_{\text{dry}} f_{\text{labile}}
$$

The implemented constants are:

$$
f_{\text{cytosol}} = 0.75 \quad m_{\text{wet}} = 1.0\thinspace\mathrm{pg}/\mu\mathrm{m}^3 \quad f_{\text{dry}} = 0.30 \quad f_{\text{labile}} = 0.35
$$

Most structural biomass is therefore unavailable, refractory, respired, or lost to the environment. Predation uses the same digestible nutrient logic instead of consuming a fixed fraction of the prey's storage capacity.

### 4.9 ECM Diffusion Grid

The extracellular matrix tracks quorum-sensing molecules, bacteriocins, cross-feeding acetate (Section 4.15), and phage on a spatial hash grid; every simulation tick applies an explicit Fickian 2D stencil<sup>[9](#ref-9)</sup> with a stability clamp<sup>[5](#ref-5)</sup> to every active grid cell:

$$
\alpha = \min\left(\alpha_{\max}, \frac{D\Delta t}{\Delta x^2}\right) \quad \alpha_{\max}=0.24
$$

Decay is applied as:

$$
m = \mathrm{clamp}(1 - \lambda\Delta t, 0, 1)
$$

The finite-difference update is:

$$
C'_{i,j}=m\left((1-4\alpha)C_{i,j}+\alpha(C_{i-1,j}+C_{i+1,j}+C_{i,j-1}+C_{i,j+1})\right)
$$

The simulation uses a $20\thinspace\mu\mathrm{m}$ ECM cell size. Quorum molecules use a higher diffusion coefficient than bacteriocins because the latter are larger antimicrobial peptides/proteins.

Concentrations below $10^{-3}$ are pruned from the sparse grid. This threshold is numerical housekeeping: those values are below the simulator's biological resolution and would otherwise accumulate as dead sparse-grid entries.

### 4.10 Gas Exchange, Fermentation, And pH

Aerobic and anaerobic environmental effects are based on each cell's actual size- and age-scaled nutrient consumption during the current update. Each bacterium stores `last_metabolic_consumption_pg` immediately after its metabolism step:

$$
G_{\Delta t}=\sum_i \min(F_i, B_i\Delta t A_i)
$$

where $G_{\Delta t}$ is consumed glucose-equivalent nutrient, $F_i$ is the cell's available stored nutrient before metabolism, $B_i$ is its size-scaled basal metabolism, and $A_i$ is the age maintenance factor.

For aerobic species with positive oxygen sensitivity $s_{\mathrm{O_2}}$, respiration follows<sup>[31](#ref-31)</sup>:

$$
\mathrm{C}_6\mathrm{H}_{12}\mathrm{O}_6 + 6\mathrm{O}_2 \rightarrow 6\mathrm{CO}_2 + 6\mathrm{H}_2\mathrm{O}
$$

The implementation converts glucose-equivalent pg into dissolved gas fractions with molar mass ratios and model capacities<sup>[37](#ref-37)</sup>:

$$
\Delta \mathrm{O}_2 = \frac{G_{\Delta t}\left(\frac{192.0}{180.16}\right)s_{\mathrm{O_2}}}{4.0\times10^7}
$$

$$
\Delta \mathrm{CO}_2 = \frac{G_{\Delta t}\left(\frac{264.0}{180.16}\right)s_{\mathrm{O_2}}}{1.5\times10^9}
$$

Carbonic acid contributes a small pH decrement:

$$
\Delta \text{acid}_{\mathrm{CO_2}}=0.02\Delta \mathrm{CO}_2
$$

For anaerobic species with negative oxygen sensitivity, the simulator uses a simplified lactic-acid pathway plus trace heterolactic carbon dioxide<sup>[31](#ref-31)</sup>:

$$
\Delta \mathrm{CO}_{2,\text{anaerobic}} = \frac{0.1G_{\Delta t}\lvert s_{\mathrm{O_2}} \rvert}{1.5\times10^9}
$$

$$
\Delta \text{acid}_{\text{lactic}} = \frac{2.5G_{\Delta t}\lvert s_{\mathrm{O_2}} \rvert}{4.0\times10^7}
$$

The update sequence is: subtract oxygen use, add carbon dioxide production, subtract acid from pH, equilibrate gases toward an open dish<sup>[9](#ref-9)</sup>, then buffer pH toward neutral<sup>[14](#ref-14)</sup>. Gas and buffer restoration are<sup>[37](#ref-37)</sup>:

$$
\mathrm{O}_{2,t+\Delta t}=\mathrm{O}_{2,t}+(0.21-\mathrm{O}_{2,t})k_{\text{gas}}\Delta t
$$

$$
\mathrm{CO}_{2,t+\Delta t}=\mathrm{CO}_{2,t}+(0-\mathrm{CO}_{2,t})k_{\text{gas}}\Delta t
$$

$$
\mathrm{pH}_{t+\Delta t}=\mathrm{pH}_t+(7.0-\mathrm{pH}_t)k_{\mathrm{buffer}}\Delta t
$$

### 4.11 Senescence

Maintenance cost rises with age<sup>[36](#ref-36)</sup>:

$$
A(t) = 1 + \frac{t}{180 \cdot L}
$$

Intrinsic death risk follows a Gompertz-style exponential penalty<sup>[12](#ref-12), [13](#ref-13)</sup>:

$$
\mu(t) = r_{\text{aging}} \Delta t \cdot e^{\alpha t}
$$

### 4.12 Competition, Bacteriocins, And T6SS

Diffusible bacteriocin production consumes stored nutrient and emits antimicrobial material into the ECM<sup>[35](#ref-35)</sup>:

$$
F_{\text{cost}} = p_{\text{bac}} k_{\text{cost}} \Delta t
$$

$$
C_{\text{bac,emit}} = 10F_{\text{cost}}
$$

Damage from environmental bacteriocin is reduced by immunity but is intentionally not allowed to reach perfect invulnerability:

$$
D_{\text{bac}} = C_{\text{local}} k_{\text{damage}} (1 - \min(0.99, 1.5I))\Delta t
$$

Type VI secretion system contact damage scales with attacker surface area<sup>[35](#ref-35)</sup>:

$$
D_{\text{T6SS}} = a_{\text{T6SS}} s^2 k_{\text{T6SS}}\Delta t(1-I_{\text{defense}})
$$

where:

$$
I_{\text{defense}} = \min\left(0.95, \frac{I_{\text{victim}}s_{\text{victim}}}{s_{\text{attacker}} + 0.01}\right)
$$

Predation requires the attacker to be at least twice the victim's linear size, which corresponds to roughly an eightfold volume advantage.

### 4.13 Named Model Constants

`bacteria/constants.py` is the canonical location for numerical model assumptions. UI-only layout values remain in UI code when they are not biological or mathematical model parameters.

Simulation timing and scale:

* `FPS = 144`: render/update target. Biological rates still use elapsed `dt`.
* `WORLD_UNIT_NAME = 'um'`: one world unit is one micrometer.
* `SCALE_BAR_REFERENCE_PIXELS = 100`: Grid scale indicator reference length.

Default species seeds:

* `DEFAULT_FOOD_CAPACITY_PG = 10.0`
* `DEFAULT_METABOLIC_RATE_PG_S = 0.1`
* `DEFAULT_CONSUMPTION_SPEED_PG_S = 1.0`
* `DEFAULT_BACTERIA_SPEED_UM_S = 40.0`
* `DEFAULT_GREEDINESS = 0.5`
* `DEFAULT_T6SS_ACTIVITY = 0.5`
* `DEFAULT_BACTERIOCIN_PRODUCTION = 0.5`
* `DEFAULT_IMMUNITY_FACTOR = 0.5`
* `DEFAULT_REPRODUCTION_THRESHOLD = 0.8`
* `DEFAULT_SHARE_THRESHOLD = 0.6`
* `DEFAULT_PREDATION_THRESHOLD = 0.9`
* `DEFAULT_BASE_SIZE_UM = 1.0`
* `DEFAULT_SENSOR_RADIUS_UM = 200.0`

Numerical tolerances:

* `NUMERIC_DISTANCE_EPSILON_UM = 1.0e-3`
* `CHEMOTAXIS_GRADIENT_EPSILON = 1.0e-3`
* `MIN_LIFESPAN_VARIANCE_DENOMINATOR = 1.0e-3`

Movement, rendering, and senescence<sup>[12](#ref-12)</sup>:

* `BACTERIA_RUN_MEAN_SEC = 1.0`
* `BACTERIA_RUN_STD_SEC = 0.2`
* `BACTERIA_RUN_MIN_SEC = 0.25`
* `BACTERIA_RUN_MAX_SEC = 3.0`
* `BACTERIA_TUMBLE_DURATION_SEC = 0.1`
* `BACTERIA_TUMBLE_MEAN_TURN_RAD = 1.19`
* `BACTERIA_TUMBLE_TURN_STD_RAD = 0.45`
* `BACTERIA_TUMBLE_TRANSLATION_FACTOR = 0.05`
* `BACTERIA_CHEMOTAXIS_RUN_EXTENSION_FACTOR = 2.5`
* `BACTERIA_ROTATIONAL_DIFFUSION_COEFFICIENT_RAD2_S = 0.16`
* `BROWNIAN_ROTATIONAL_VARIANCE_FACTOR = 2.0`
* `BACTERIA_GRAZING_SPEED_FACTOR = 0.5`
* `BACTERIA_BOUNDARY_OFFSET_UM = 0.1`
* `BACTERIA_BOUNDARY_REFLECTION_NOISE_RAD = 0.2`
* `BACTERIA_STARVED_RENDER_SCALE = 0.9`
* `BACTERIA_FULL_RENDER_SCALE_RANGE = 0.1`
* `LIFESPAN_VARIANCE_MIN = 0.4`
* `LIFESPAN_VARIANCE_MEAN = 1.0`
* `LIFESPAN_VARIANCE_STD = 0.2`
* `SENESCENCE_BASE_LIFESPAN_SEC = 180.0`
* `GOMPERTZ_MORTALITY_ALPHA = 0.04`
* `FISSION_COST_PG_PER_UM3 = 40.0`

Chemotaxis and nutrient geometry:

* `CHEMOTAXIS_CORE_RADIUS_UM = 5.0`
* `CHEMOTAXIS_REFERENCE_PG_PER_UM2 = 1.0e-3`
* `CHEMOTAXIS_STEER_RATE_RAD_S = 1.5`
* `MIN_FOOD_NODE_AMOUNT_PG = 1.0e-3`
* `FOOD_NODE_MIN_RADIUS_UM = 0.5`
* `FOOD_NODE_MIN_VOLUME_PG = 0.1`
* `FOOD_EATING_DISTANCE_MARGIN_UM = 1.0`
* `FOOD_EATING_DISTANCE_GRACE_UM = 3.0`

Necromass and death recycling:

* `CELL_WET_MASS_PG_PER_UM3 = 1.0`
* `CELL_DRY_MASS_FRACTION = 0.30`
* `LABILE_NECROMASS_FRACTION = 0.35`
* `CYTOSOL_RECOVERY_FRACTION = 0.75`

Environmental fitness, gas, and pH:

* `PH_TOLERANCE_GAUSSIAN_SD = 1.5`
* `AEROBE_O2_SATURATION_LEVEL = 0.10`
* `ANAEROBE_O2_TOXICITY_THRESHOLD = 0.02`
* `ENV_TOXIC_DEATH_RATE_PER_SEC = 0.05`
* `ENV_100_O2_PG = 4.0e7`
* `ENV_100_CO2_PG = 1.5e9`
* `GLUCOSE_MOLAR_MASS_G = 180.16`
* `O2_SIX_MOLES_MASS_G = 192.0`
* `CO2_SIX_MOLES_MASS_G = 264.0`
* `CARBONIC_ACID_PH_EFFECT = 0.02`
* `HETEROLACTIC_CO2_TRACE_FRACTION = 0.1`
* `LACTIC_ACID_PH_EFFECT = 2.5`
* `GAS_DIFFUSION_RATE_SEC = 0.288`
* `BUFFER_RESTORATION_RATE_SEC = 0.144`
* `ATMOSPHERIC_O2_FRACTION = 0.21`
* `ATMOSPHERIC_CO2_FRACTION = 0.0`
* `NEUTRAL_PH = 7.0`

Competition and contact:

* `BACTERIOCIN_PRODUCTION_COST_RATE = 0.5`
* `BACTERIOCIN_ECM_SECRETION_MULTIPLIER = 10.0`
* `BACTERIOCIN_MAX_IMMUNITY_MITIGATION = 0.99`
* `BACTERIOCIN_IMMUNITY_SCALING = 1.5`
* `BACTERIOCIN_DAMAGE_COEFFICIENT = 0.1`
* `T6SS_DAMAGE_SCALING_FACTOR = 5.0`
* `T6SS_MAX_IMMUNITY_DEFENSE = 0.95`
* `BACTERIA_CONTACT_MARGIN_UM = 2.0`
* `T6SS_IMMUNITY_SIZE_EPSILON_UM = 0.01`
* `PHAGOCYTOSIS_SIZE_RATIO_THRESHOLD = 2.0`
* `CONSUMED_PREY_SENTINEL_FOOD_LEVEL = -9999.0`

ECM and sparse-grid performance:

* `ECM_ALPHA_MAX = 0.24`
* `ECM_QS_DIFFUSION_UM2_S = 400.0`
* `ECM_QS_DECAY_S = 0.00385`
* `ECM_BAC_DIFFUSION_UM2_S = 60.0`
* `ECM_BAC_DECAY_S = 0.00578`
* `ECM_MIN_CONCENTRATION = 1.0e-3`
* `ECM_CELL_SIZE_UM = 20`

Spatial bounds and lookup:

* `FOOD_SPATIAL_HASH_CELL_SIZE_UM = 50`
* `BACTERIA_INTERACTION_HASH_CELL_SIZE_UM = 20`
* `BACTERIA_INTERACTION_AXIS_CUTOFF_UM = 15`
* `PETRI_RADIUS_PADDING_UM = 100`

Initial inoculum spawning (Section 4.23):

* `INITIAL_SPAWN_MODE = 'clustered'` — `'clustered'` | `'uniform_disc'` | `'central'`
* `INITIAL_SPAWN_LUMEN_BIAS = 0.9`
* `INITIAL_SPAWN_LUMEN_RADIUS_FRACTION` $=$ `MUCUS_LAYER_RADIUS_FRACTION` $(= 0.75)$
* `INITIAL_SPAWN_CLUSTER_COUNT = 3`
* `INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION = 0.06`
* `INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION = 0.15`
* `INITIAL_SPECIES_ARRANGEMENT = 'independent'` — `'independent'` | `'mixed'` | `'segregated'`

Temperature kinetics (Ratkowsky model)<sup>[23](#ref-23)</sup>:

* `RATKOWSKY_B_SQRT_PER_C = 0.0410`
* `RATKOWSKY_C_PER_C = 0.161`
* `DEFAULT_TEMP_MIN_C = 4.0`
* `DEFAULT_TEMP_MAX_C = 47.0`
* `DEFAULT_ENV_TEMPERATURE_C = 37.0`
* `ENV_TEMPERATURE_MIN_C = 0.0`
* `ENV_TEMPERATURE_MAX_C = 60.0`
* `RATKOWSKY_NORMALIZATION_SAMPLES = 200`

Cross-feeding (acetate overflow, see Section 4.15)<sup>[30](#ref-30)</sup>:

* `ACETATE_MOLAR_MASS_G = 60.05`
* `ACETATE_PER_GLUCOSE_MAX_MASS_RATIO` $= (2 \times 60.05) / 180.16 \approx 0.667$
* `GLUCOSE_AEROBIC_ATP_PER_MOL = 32.0`
* `ACETATE_AEROBIC_ATP_PER_MOL = 10.0`
* `ACETATE_ENERGY_YIELD_GLUCOSE_EQUIV` $\approx 0.938$
* `O2_PER_ACETATE_MASS_RATIO = 64.0 / 60.05`
* `CO2_PER_ACETATE_MASS_RATIO = 88.02 / 60.05`
* `ECM_ACETATE_DIFFUSION_UM2_S = 1000.0`
* `ECM_ACETATE_DECAY_S = 0.001`
* `ACETATE_ECM_UNITS_PER_PG = 100.0`
* `DEFAULT_ACETATE_SECRETION = 0.0`
* `DEFAULT_ACETATE_UPTAKE = 0.0`

Lag phase (Baranyi & Roberts dynamic state, see Section 4.16)<sup>[2](#ref-2)</sup>:

* `DEFAULT_BARANYI_Q0 = 0.05` (fresh-inoculum physiological state; $h_0 = \ln(1 + 1/q_0) \approx 3.04$)
* `BARANYI_Q_STATE_MAX` $= 1.0 \times 10^{6}$ (numerical cap where $\alpha \gt 0.999999$)
* `BARANYI_Q0_MIN` $= 1.0 \times 10^{-4}$ (floor; guarantees a finite lag)
* `BARANYI_Q0_MAX` $= 1.0 \times 10^{3}$ (ceiling; effectively disables the lag)

Phenotypic switching (persister bet-hedging, see Section 4.17)<sup>[1](#ref-1), [20](#ref-20)</sup>:

* `DEFAULT_PERSISTER_SWITCH = 0.0` (per-species switching propensity in $[0,1]$; $0$ disables persistence)
* `DEFAULT_PERSISTER_ENTER_THRESH = 0.15` (starvation fullness threshold $S_1$)
* `DEFAULT_PERSISTER_EXIT_THRESH = 0.50` (recovery fullness threshold $S_2$)
* `PERSISTER_SPONTANEOUS_SWITCH_RATE_PER_SEC = 0.002` (baseline bet-hedging entry rate $k_0^{+}$)
* `PERSISTER_STRESS_SWITCH_RATE_PER_SEC = 0.2` (starvation-induced entry rate $k_s^{+}$)
* `PERSISTER_WAKE_RATE_PER_SEC = 0.2` (nutrient-triggered resuscitation rate $k_r^{-}$)
* `PERSISTER_SPONTANEOUS_WAKE_RATE_PER_SEC = 0.005` (baseline reversion rate $k_0^{-}$)
* `PERSISTER_METABOLIC_FACTOR = 0.1` (dormant maintenance burn as a fraction of basal metabolism)
* `PERSISTER_RESISTANCE_FACTOR = 0.1` (incoming bacteriocin/T6SS/senescence hazard multiplier, $\approx 10\times$ tolerance)
* `PERSISTER_THRESHOLD_MIN_GAP = 0.05` (minimum $S_2 - S_1$ separation enforced by the constructor)

These rates are sim-accelerated (the simulation runs in real seconds, not hours) so the bet-hedging dynamics are observable within a session while preserving the relative ordering of the spontaneous, stress-induced, and resuscitation pathways.

Phage predation and horizontal gene transfer (see Section 4.18)<sup>[8](#ref-8), [19](#ref-19), [22](#ref-22), [26](#ref-26)</sup>:

* `ECM_PHAGE_DIFFUSION_UM2_S = 500.0` (free-virion diffusion coefficient in the ECM phage layer)
* `ECM_PHAGE_DECAY_S = 0.02` (virion inactivation / loss-of-infectivity rate per second)
* `PHAGE_ECM_UNITS_PER_VIRION = 1.0` (ECM grid units per virion; kept $O(1)$ so a single adsorption is meaningful)
* `DEFAULT_PHAGE_ADSORPTION_RATE = 0.02` (mass-action adsorption rate constant $k$ per virion-unit per second)
* `DEFAULT_PHAGE_LATENT_PERIOD_SEC = 3.0` (eclipse + latent period $\tau$ before lysis, sim-accelerated)
* `DEFAULT_PHAGE_BURST_SIZE = 20.0` (progeny virions $\beta$ released per lysed cell)
* `DEFAULT_PHAGE_INFLOW = 0.0` (global seed-packet introduction rate in packets/sec; $0$ = no outbreak)
* `PHAGE_SEED_PACKET_VIRIONS = 10.0` (virions deposited per introduced seed packet)
* `PHAGE_MIN_ADSORBABLE` $= 1.0 \times 10^{-3}$ (local phage below this ECM prune floor is negligible for adsorption)
* `PHAGE_RESISTANCE_MUTATION_RATE = 0.001` (spontaneous susceptible $\to$ resistant mutation probability per division)
* `DEFAULT_PHAGE_RESISTANCE = 0.0` (per-species initial resistant fraction in $[0,1]$)
* `DEFAULT_HGT_COMPETENCE = 0.0` (per-species HGT donate/receive competence in $[0,1]$; $0$ disables transfer)
* `HGT_TRANSFER_RATE_PER_SEC = 0.1` (contact-dependent horizontal transfer rate constant)
* `HGT_METABOLIC_BURDEN = 0.2` (fractional rise in basal metabolism from carrying the transferred plasmid)

These rates are likewise sim-accelerated: natural latent periods span tens of minutes, but here the lytic cycle resolves in seconds so an epidemic crash and the resistant rebound are visible within a session.

### 4.14 Temperature Kinetics (Ratkowsky Square-Root Model)

Temperature is the single most influential environmental driver in predictive microbiology.
BacSandbox scales growth, metabolism, nutrient uptake, and motility by a temperature multiplier
derived from the four-parameter Ratkowsky square-root model<sup>[23](#ref-23)</sup>, which is biologically validated across the entire biokinetic temperature range:

$$
\sqrt{\mu(T)} = b\thinspace(T - T_{\min})\thinspace\bigl(1 - e^{\thinspace c\thinspace(T - T_{\max})}\bigr)
$$

where $T$ is the environment temperature in degrees Celsius, $T_{\min}$ and $T_{\max}$ are the
species' theoretical minimum and maximum growth temperatures (its *cardinal* temperatures),
$b$ is the low-temperature square-root slope, and $c$ governs the steepness of high-temperature
thermal inactivation. Outside the open interval $(T_{\min}, T_{\max})$ growth is zero.

The raw rate $\mu(T) = \bigl(\sqrt{\mu(T)}\bigr)^2$ is converted into a dimensionless multiplier in
$[0, 1]$ by normalizing against each species' own optimum $\sqrt{\mu_{\text{peak}}}$, located once at
construction by uniform sampling across $[T_{\min}, T_{\max}]$:

$$
f_{T}(T) = \left( \frac{\sqrt{\mu(T)}}{\sqrt{\mu_{\text{peak}}}} \right)^{2} \qquad f_{T} \in [0, 1]
$$

Because the multiplier is normalized to the optimum, the slope parameter $b$ cancels and only the
curve *shape* (set by $T_{\min}$, $T_{\max}$, and $c$) matters. The factor reaches $1.0$ at the
species optimum and falls to $0.0$ at its cardinal limits. It is applied as follows:

* **Reproduction:** $f_{T}$ multiplies the intrinsic Verhulst growth rate $r$ alongside the
  carrying-capacity term and environmental fitness (see Section 4.7).
* **Metabolism:** basal nutrient demand is scaled by $f_{T}$, so chilled cultures deplete their
  nutrient reserves far more slowly, reproducing the long shelf life of refrigerated cells.
* **Nutrient uptake:** membrane transporter throughput (`consumption_speed`) is scaled by $f_{T}$.
* **Motility:** flagellar swimming speed is scaled by $f_{T}$, since motor output is enzyme-driven.

Near and beyond the cardinal limits $f_{T} \to 0$, driving cells toward thermal dormancy rather than
immediate death. Different cardinal temperatures yield distinct thermal niches: a psychrophile
(low $T_{\min}$, $T_{\max}$) peaks at a colder optimum than a mesophile, enabling temperature-based
niche separation between species.


### 4.15 Metabolic Cross-Feeding (Acetate Overflow)

Real microbial communities are knit together by *cross-feeding*: the metabolic waste of one
organism is the food of another. The most studied example in *Escherichia coli* is the **acetate
switch**<sup>[30](#ref-30)</sup>: during fast aerobic growth on
glucose the cell's central metabolism outpaces its respiratory capacity and excretes *overflow*
acetate; when glucose runs low (or in neighbouring taxa) that acetate is re-imported and oxidized
as a secondary carbon source. BacSandbox models acetate as a third diffusible ECM layer (alongside
quorum-sensing inducers and bacteriocins) coupling acetate **producers** to acetate **consumers**.

**Secretion.** A producer with secretion trait $s \in [0, 1]$ converts a fraction of the glucose it
took up on the previous tick, $U_{\text{glc}}$ (pg), into excreted acetate, bounded by the
stoichiometric ceiling of one glucose splitting into two acetate molecules by mass:

$$
A_{\text{sec}} = \min\negthinspace\Bigl(\thinspace s \thinspace U_{\text{glc}} \thinspace Y_{a/g}, F \thinspace\Bigr) \qquad
Y_{a/g} = \frac{2\thinspace M_{\text{ace}}}{M_{\text{glc}}} = \frac{2 \times 60.05}{180.16} \approx 0.667
$$

where $F$ is the cell's current internal nutrient store (acetate cannot be excreted from carbon the
cell does not hold) and $M_{\text{ace}}, M_{\text{glc}}$ are the molar masses of acetate and glucose.
The secreted mass is removed from the cell's store and deposited into the local acetate ECM cell.

**Uptake.** A consumer with uptake trait $u \in [0, 1]$ oxidizes locally available acetate at a rate
limited by its transporter budget, by local availability $A_{\text{avail}}$, and by remaining storage
room. Because acetate oxidation is aerobic and enzyme-driven, uptake is gated by both the species'
Ratkowsky temperature factor $f_T$ (Section 4.14) and a dissolved-oxygen gate:

$$
g_{\mathrm{O_2}} = \min\negthinspace\left(1, \frac{[\mathrm{O}_2]}{[\mathrm{O}_2]_{\text{sat}}}\right) \qquad
B = c_{\text{eat}} \thinspace u \thinspace \Delta t \thinspace f_T \thinspace g_{\mathrm{O_2}}
$$

$$
A_{\text{con}} = \max\negthinspace\Bigl(0, \min\bigl(A_{\text{avail}}, B, R \bigr)\Bigr) \qquad
R = \frac{F_{\max} - F}{Y_E}
$$

where $c_{\text{eat}}$ is the cell's size-scaled `consumption_speed`, $\Delta t$ the real elapsed
time, and $R$ the storage room expressed in acetate mass. Consumed acetate is removed from the ECM
and converted to stored nutrient at the **glucose-energy-equivalent yield** $Y_E$, which reflects the
lower ATP return of acetate relative to glucose:

$$
Y_E = \frac{\text{ATP}_{\text{ace}} / M_{\text{ace}}}{\text{ATP}_{\text{glc}} / M_{\text{glc}}}
    = \frac{10 / 60.05}{32 / 180.16} \approx 0.938
$$

so $\Delta F = A_{\text{con}} \thinspace Y_E$. The ATP yields ($\sim 32\thinspace\mathrm{ATP}/\mathrm{glucose}$<sup>[25](#ref-25)</sup>; $\sim 10\thinspace\mathrm{ATP}/\mathrm{acetate}$
via acetyl-CoA → TCA cycle + electron transport chain) make acetate a slightly less
energy-dense carbon source per unit mass, as observed experimentally.

**Gas coupling.** Acetate oxidation follows $\mathrm{C}_2\mathrm{H}_4\mathrm{O}_2 + 2\mathrm{O}_2 \rightarrow 2\mathrm{CO}_2 + 2\mathrm{H}_2\mathrm{O}$<sup>[31](#ref-31)</sup>.
Each picogram of acetate consumed therefore draws down dissolved $\mathrm{O_2}$ and releases
$\mathrm{CO_2}$ (with a corresponding carbonic-acid pH effect) independently of the glucose-driven
respiration in Section 4.10:

$$
\Delta[\mathrm{O}_2] = -\thinspace A_{\text{con}}\thinspace\frac{64.0}{M_{\text{ace}}}\Big/ C_{\mathrm{O_2}} \qquad
\Delta[\mathrm{CO}_2] = +\thinspace A_{\text{con}}\thinspace\frac{88.02}{M_{\text{ace}}}\Big/ C_{\mathrm{CO_2}}
$$

Here, $64.0\thinspace\mathrm{g}$ and $88.02\thinspace\mathrm{g}$ are the corresponding oxygen and carbon dioxide masses per mole of acetate.

The reaction consumes $2\thinspace\mathrm{mol}\thinspace\mathrm{O}_2$ and produces $2\thinspace\mathrm{mol}\thinspace\mathrm{CO}_2$ per mole of acetate.

The modeled dissolved-gas saturation capacities are $C_{\mathrm{O_2}}$ and $C_{\mathrm{CO_2}}$.

**Numerical note.** Picogram acetate masses are scaled into ECM grid units by a constant factor
(`ACETATE_ECM_UNITS_PER_PG`) so traces survive the sparse-grid pruning floor; the identical factor is
inverted on uptake so secreted and consumed mass are conserved. The acetate layer diffuses on the
same explicit Fickian stencil as the other ECM layers (Section 4.9). With the default traits
($s = u = 0$) cross-feeding is inert and pre-M2 behaviour is reproduced exactly.


### 4.16 Lag Phase (Baranyi & Roberts Dynamic Physiological State)

Newly inoculated cells do not divide immediately: they must first synthesize the ribosomes, enzymes,
and membrane components needed for growth. This delay is the **lag phase**. BacSandbox models it with
the dynamic approach of Baranyi & Roberts<sup>[2](#ref-2)</sup>, which augments each cell with a dimensionless
*physiological state* variable $q$ that captures the "work still to be done" before division.

The realized growth rate is the maximum (temperature-scaled) rate multiplied by an **adjustment
function** that depends only on $q$:

$$
\alpha(q) = \frac{q}{1 + q}
$$

When $q$ is small (fresh inoculum) $\alpha \approx 0$ and reproduction is suppressed; as the cell
adapts $q \to \infty$ and $\alpha \to 1$, restoring the full rate. The physiological state itself
grows in proportion to the specific growth rate,

$$
\frac{dq}{dt} = \mu_{\max}\thinspace \phi(T)\thinspace q
$$

where $\mu_{\max}$ is the species' base reproduction rate (`rep`) and $\phi(T)$ is the Ratkowsky
temperature factor (Section 4.14). Because the right-hand side is linear in $q$, the equation is
integrated **analytically** over each tick rather than with forward Euler, so the result is exact and
frame-rate independent:

$$
q_{t+\Delta t} = q_t \thinspace e^{\mu_{\max}\thinspace\phi(T)\thinspace\Delta t}.
$$

This temperature coupling reproduces the observed fact that cold cultures lag far longer than
cultures at their thermal optimum, since $\phi(T) \to 0$ slows adaptation.

**Lag duration.** Integrating the adjustment function gives the classic Baranyi lag time

$$
\lambda = \frac{h_0}{\mu_{\max}} \qquad h_0 = \ln\negthinspace\left(1 + \frac{1}{q_0}\right)
$$

where $q_0$ is the initial physiological state (`baranyi_q0`) and $h_0$ is the dimensionless "work to
be done". A smaller $q_0$ means a longer lag. The reproduction probability per tick (Section 4.14)
therefore becomes

$$
P_{rep} = (\text{rep} \cdot \Delta t)\thinspace (1 - N/K)\thinspace f_{env}\thinspace \phi(T)\thinspace \alpha(q).
$$

**Inheritance and reset.** Only the **initial inoculum** starts at $q_0$; daughter cells produced by
fission inherit the parent's already-elevated $q$, so an actively growing population does not re-enter
lag on every division. A fresh inoculum (or a culture seeded into a new patch via reset) therefore
shows the characteristic flat-then-exponential growth curve.

**Numerical bounds.** $q_0$ is clamped to $[\thinspace 10^{-4}, 10^{3}\thinspace]$ so the lag is always finite and
non-infinite, and the running state $q$ is capped at `BARANYI_Q_STATE_MAX` $= 10^6$ (where
$\alpha \gt 0.999999$) to prevent floating-point overflow on long runs. Setting `baranyi_q0` to its
ceiling makes $\alpha \approx 1$ from the outset, reproducing the pre-M3 (no-lag) growth behaviour.


### 4.17 Phenotypic Switching (Persister Bet-Hedging)

Even in clonal, genetically identical populations a small fraction of cells stochastically enter a
**persister** state: a dormant, slow-growing, stress-tolerant phenotype that survives antimicrobial
exposure lethal to actively dividing cells. When conditions improve the persisters resuscitate back
to the proliferative state. This reversible bistable switch is the canonical explanation for
antibiotic *tolerance* (as distinct from inherited resistance)<sup>[1](#ref-1), [20](#ref-20)</sup>. BacSandbox models it as a two-state
Markov process whose transition rates depend on the cell's internal **fullness**
$f = \min(1, \max(0, \text{food}/\text{capacity}))$, used as a self-consistent proxy for sustained
local nutrient availability.

**Entry (proliferative $\rightarrow$ persister).** Below the starvation threshold $S_1$ the entry
rate combines a constant bet-hedging term with a stress term that grows as reserves empty:

$$
k^{+}(f) = k_0^{+} + \mathbb{1}[f \lt S_1]\thinspace k_s^{+}\left(1 - \frac{f}{\max(S_1, \varepsilon)}\right)
$$

**Exit (persister $\rightarrow$ proliferative).** Above the recovery threshold $S_2$ the resuscitation
rate combines a constant reversion term with a recovery term that grows as reserves refill:

$$
k^{-}(f) = k_0^{-} + \mathbb{1}[f \gt S_2]\thinspace k_r^{-}\thinspace\frac{f - S_2}{\max(1 - S_2, \varepsilon)}
$$

**Per-tick probability.** Each transition is realized with the exact exponential survival
complement over the real elapsed time $\Delta t$, so the dynamics are frame-rate independent:

$$
P_{\text{switch}} = 1 - e^{-\thinspace k\thinspace\Delta t}.
$$

The entry rate is additionally scaled by the per-species propensity $s \in [0,1]$
(`persister_switch`), i.e. $k_{\text{eff}}^{+} = s\thinspace k^{+}(f)$, so $s = 0$ disables persistence
entirely and reproduces the pre-M4 behaviour exactly. The two fullness thresholds are constrained by
the constructor so that $S_2 \ge S_1 +$ `PERSISTER_THRESHOLD_MIN_GAP`, preventing the
enter/exit bands from overlapping and causing rapid chatter.

**Persister phenotype.** While dormant, a cell:

* **Halts division** entirely (fission is skipped in `BacteriaSpecies.step`).
* **Throttles maintenance metabolism** to a fraction `PERSISTER_METABOLIC_FACTOR` $= 0.1$ of its
  basal demand, conserving stored nutrient and extending survival during the crash.
* **Gains stress tolerance**: incoming bacteriocin exposure, T6SS contact damage, and the combined
  environmental/senescence death hazard are all multiplied by `PERSISTER_RESISTANCE_FACTOR` $= 0.1$
  ($\approx 10\times$ tolerance), modelling the dramatically reduced kill rate of dormant cells.

Persisters still move and feed normally, which is what allows their fullness to recover and trigger
resuscitation. Because dividing cells are always proliferative, daughter cells are born proliferative
and the switch must be re-entered stochastically — no inheritance of the persister state is required.

### 4.18 Phage Predation and Horizontal Gene Transfer (SIR Lytic Cycle)

Lytic bacteriophages are the most abundant biological entities on Earth and the dominant top-down
control on bacterial populations. BacSandbox models them not as discrete agents but as a diffusible
**free-virion field** carried in the ECM (`phage` layer), coupled to a per-cell susceptible /
infected / resistant (SIR) state machine. This reuses the existing reaction–diffusion infrastructure
and keeps the model lightweight while reproducing the classic predator–prey boom–bust oscillation<sup>[19](#ref-19)</sup>.

**Adsorption (susceptible $\rightarrow$ infected).** A susceptible (non-resistant, uninfected) cell
adsorbs phage by mass action against the local free-virion concentration $[P]$ sampled at its
position. Over the real elapsed time $\Delta t$ the infection probability is the exponential survival
complement, so the dynamics are frame-rate independent:

$$
P_{\text{adsorb}} = 1 - e^{-\thinspace k\thinspace[P]\thinspace\Delta t}
$$

where $k$ is `DEFAULT_PHAGE_ADSORPTION_RATE`. A successful adsorption removes one virion from the
field, sets the cell to the infected state, and starts its latent timer at $\tau$
(`DEFAULT_PHAGE_LATENT_PERIOD_SEC`).

**Latency and lysis (infected $\rightarrow$ lysed).** Infected cells cannot divide — the virus
commandeers the replication machinery. The latent timer counts down by $\Delta t$ each tick; when it
reaches zero the host **lyses**: it dies, deposits labile necromass as a food node (cell debris, with
no cytosol recovery since the cytoplasm was hijacked), and releases a burst of $\beta$
(`DEFAULT_PHAGE_BURST_SIZE`) progeny virions into the local phage field. This one-step growth
structure (a fixed latent period followed by a discrete burst) follows Ellis & Delbrück<sup>[8](#ref-8)</sup>.

**Free-virion field.** Released virions diffuse (`ECM_PHAGE_DIFFUSION_UM2_S`) and lose infectivity
(`ECM_PHAGE_DECAY_S`) via the same Fickian stencil used for the other ECM layers (Section 4.9),
spreading the epidemic spatially and decaying it where hosts are absent.

**Resistance, mutation, and rebound.** Each cell is born resistant with probability equal to its
species' `phage_resistance` trait (a heritable, chromosomal fraction). Resistant cells never adsorb
phage, so they survive the crash. On every division a susceptible lineage can spontaneously mutate to
resistant with probability `PHAGE_RESISTANCE_MUTATION_RATE` (a Luria–Delbrück-style fluctuation)<sup>[21](#ref-21)</sup>,
seeding the post-crash recovery. Daughters inherit the parent's resistance state and carriage cost.
The net result is the textbook sequence: phage outbreak $\rightarrow$ susceptible population crash
$\rightarrow$ resistant-only survivors $\rightarrow$ rebound.

**Horizontal gene transfer (HGT).** Resistance also spreads laterally<sup>[22](#ref-22)</sup>. When a resistant donor and a
susceptible recipient are in physical contact (the same contact test used for sharing/T6SS), and both
species are HGT-competent (`hgt_competence` $\gt 0$), the trait transfers by a conjugation-like process
with probability

$$
P_{\text{HGT}} = 1 - e^{-\thinspace k_H\thinspace c_{\text{donor}}\thinspace c_{\text{recip}}\thinspace\Delta t}
$$

where $k_H$ is `HGT_TRANSFER_RATE_PER_SEC` and $c$ are the two cells' competences. The recipient
becomes resistant and inherits a **plasmid carriage cost**: its basal maintenance metabolism is
multiplied by $1 +$ `HGT_METABOLIC_BURDEN`, modelling the well-documented fitness penalty of
maintaining a resistance plasmid<sup>[26](#ref-26)</sup>. This cost is inherited by daughters, so
in the absence of phage pressure the burdened lineage is slightly outcompeted — capturing the
cost–benefit tension that governs the spread of mobile resistance elements.

**Outbreak control.** The global `phage_inflow` environment parameter seeds free virions via a
Poisson process ($\lambda\thinspace\Delta t$ packets of `PHAGE_SEED_PACKET_VIRIONS` at random petri
locations), so a user can trigger or halt an epidemic live from the Globals tab. A value of $0$
(the default) means no outbreak and exactly reproduces the pre-M5 behaviour.


### 4.19 Diet-Driven Substrate Metabolism (Production pg Path)

The production Pygame simulation extends its abstract pg glucose-equivalent nutrient pool with
**three discrete Own-model substrate classes**, each nutrient patch (`FoodNode`) carrying a
`substrate_id` $s \in \lbrace 0, 1, 2 \rbrace$. This is the production/legacy categorical food-drop
path used by `Simulation`; it is not the empirically sourced recipe translator or reduced molar API
described in Section 10.4.<sup>[48](#ref-48)</sup> Published work supports the qualitative premise that
microbiota-accessible carbohydrate composition can alter substrate utilization and community
ecology, but it does not supply these three categories or their numerical calibrations.<sup>[10](#ref-10), [11](#ref-11), [27](#ref-27)</sup>
Raman et al. separately supports the qualitative context of impaired microbiota development in
undernourished children, not a numerical diet calibration.<sup>[24](#ref-24)</sup>

- $s = 0$ — **generic calorie**: a simple, universally fermentable sugar (e.g. rice/starch hydrolysate).
- $s = 1$ — **complex glycan**: dietary fibre only partially fermentable by the gut community.
- $s = 2$ — **MDCF prebiotic**: a Microbiota-Directed Complementary Food substrate that selectively rewards competent taxa.

Necromass, decomposition, and all legacy food default to $s = 0$, so simulations with no diet
configured reproduce the original dynamics exactly.<sup>[48](#ref-48)</sup>

**Dimensionless throughput/capacity multipliers.** Each species declares a multiplier vector
$\mathbf{a} = (a_0, a_1, a_2)$ with $a_s \ge 0$ (the columns `aff_calorie`, `aff_glycan`,
`aff_mdcf`). Despite the historical `aff_*` names, these values are not measured molecular binding
affinities. In the production pg path, each value gates chemotactic targeting and feeding
throughput for its substrate class; a value of $0$ makes that class inaccessible. In the opt-in
reduced molar API, the same dimensionless value instead scales the Own-model base $V^{\max}$ for
that class, as documented in Section 10.4.<sup>[47](#ref-47), [48](#ref-48)</sup> Chemotaxis weights
the score with $a_s$:

$$
\mathrm{score}_s = g_{\text{eff}} \thinspace a_s \thinspace \ln\negthinspace\left(1 + \frac{C}{C_{\text{ref}}}\right)
$$

where $C$ is the local patch concentration and $C_{\text{ref}}$ the receptor reference
concentration. Substrates with $a_s = 0$ are skipped entirely in both targeting and feeding.

**Multiplier-gated uptake with energy yield.** Within one time step a cell has a shared transporter
budget $U_{\max} = v_c \thinspace \Delta t \thinspace \phi_T$ (mass, pg), where $v_c$ is `consumption_speed` and
$\phi_T$ the temperature factor. For each reachable patch of class $s$, the ingested **mass** is

$$
\Delta m = \min\negthinspace\big(R \thinspace a_s, R, m_{\text{node}}, h_s\big) \qquad R = U_{\max} - \textstyle\sum \Delta m_{\text{prior}}
$$

with $R$ the remaining budget and $h_s$ the storage headroom expressed back in substrate mass.
Mass is conserved on the patch ($m_{\text{node}} \leftarrow m_{\text{node}} - \Delta m$), but the
**energy** added to the cell's internal store is scaled by a substrate-specific yield $Y_s$:

$$
\Delta E = \Delta m \thinspace Y_s
$$

The fractions $Y_0 = 1.0$, $Y_1 = 0.7$, and $Y_2 = 0.85$ are **Own (BacSandbox
model)** energy-yield calibrations, not measured substrate-specific extraction efficiencies; all
satisfy $0 \lt Y_s \le 1$. They operationalize qualitative differences in carbohydrate
accessibility while the multiplier vector determines which taxa can benefit from each class.<sup>[10](#ref-10), [27](#ref-27), [48](#ref-48)</sup>

**Own-model categorical nutrient-drop scenarios.** The environment's `diet_profile` selects a
categorical probability distribution over substrate classes for each fresh nutrient drop. The
implementation labels `malnutrition`, `standard`, and `mdcf` are scenario names, not empirical
trial formulations. Their per-profile fractions sum to $1$:<sup>[48](#ref-48)</sup>

| Profile | calorie $p_0$ | glycan $p_1$ | MDCF $p_2$ |
|---|---|---|---|
| `malnutrition` | $0.85$ | $0.10$ | $0.05$ |
| `standard` | $0.60$ | $0.30$ | $0.10$ |
| `mdcf` | $0.40$ | $0.25$ | $0.35$ |

These ratio vectors and the three energy yields above are **Own (BacSandbox model)** calibrations;
they are not ingredient percentages from the Chen et al. trial.<sup>[42](#ref-42), [48](#ref-48)</sup>
They are also distinct from the versioned `malnutrition`, `rusf`, and `mdcf_2` dry-recipe profiles
in `bacteria/models/diet_profiles.json` and Section 10.4. Each production-path drop class is drawn
by inverse-CDF sampling of the selected ratio vector, making `malnutrition` calorie-dominated and
`mdcf` comparatively enriched in the modeled prebiotic class by construction. The active scenario
is selectable in the **Globals** tab: clicking the `diet_profile` cell (or pressing Enter on its row) cycles
`malnutrition` $\to$ `standard` $\to$ `mdcf` $\to$ `malnutrition`, applied on the next reset/reload.


### 4.20 Spatial Heterogeneity: Lumen vs Mucus Layer

The circular petri dish is geometrically partitioned into two concentric zones that abstract the
radial organisation of the intestinal cross-section<sup>[6](#ref-6)</sup>: a central, well-mixed
**lumen** and a thin outer **mucus** annulus. The concentric geometry and all numerical zone
constants are project-specific mechanisms, not measurements of intestinal dimensions or flow.<sup>[48](#ref-48)</sup>
For a point at world position $(x, y)$, with petri centre $(c_x, c_y)$ and radius $R$, the radial fraction is

$$
\rho = \frac{\sqrt{(x - c_x)^2 + (y - c_y)^2}}{R}.
$$

A point lies in the mucus layer when $\rho \ge$ `MUCUS_LAYER_RADIUS_FRACTION` (default $0.75$), so
the inner $75\\%$ of the radius is lumen and the outer $25\\%$ is the mucus annulus. Four zone-specific
rules differentiate the modeled niches.<sup>[48](#ref-48)</sup>

**1. Mucus-associated motility reduction.** Published imaging supports spatial organization and
mucus association in the gut, but it does not determine this simulator's motility multiplier.<sup>[18](#ref-18)</sup>
Before forward translation, a cell whose current position satisfies
$\rho \ge$ `MUCUS_LAYER_RADIUS_FRACTION` has its per-step speed scaled by the Own-model mucus
motility factor $m =$ `MUCUS_MOTILITY_FACTOR` $= 0.3$, giving effective speed
$v_{\text{eff}} = v\thinspace m$. Slower modeled motion is intended to contribute to a candidate
low-clearance refuge; it is not evidence that persistence improves.<sup>[48](#ref-48)</sup>

**2. Nutrient washout (peristaltic-clearance abstraction).** Free food mass is reduced per real
second at a zone-specific Own-model rate. For a node with
mass $m_t$ in a zone of washout rate $k$,

$$
m_{t+\Delta t} = m_t\thinspace(1 - k\thinspace\Delta t)
$$

with $k =$ `LUMEN_FOOD_WASHOUT_RATE_PER_SEC` $(= 0.02)$ in the lumen and $k =$
`MUCUS_FOOD_WASHOUT_RATE_PER_SEC` $(= 0.002)$ in the mucus. The lumen is flushed an order of
magnitude faster in the implementation; the lower mucus rate is intended to retain modeled food
longer, not to reproduce a measured clearance rate.<sup>[48](#ref-48)</sup>

**3. Drop deposition.** Fresh dietary substrate is delivered chiefly to the well-mixed lumen by
luminal flow. A drop landing in the mucus annulus ($\rho \ge$ `MUCUS_LAYER_RADIUS_FRACTION`)
retains only a fraction `MUCUS_DROP_AMOUNT_FACTOR` $(= 0.5)$ of its sampled mass, while lumen
drops are deposited in full. This deposition rule and factor are Own-model choices.<sup>[48](#ref-48)</sup>

**4. Competition bias.** Under hard carrying-capacity pressure the well-mixed lumen imposes
an Own-model culling bias, so luminal cells are preferentially selected for removal. Culling uses
weighted sampling without replacement<sup>[7](#ref-7)</sup>: each
candidate cell draws a key $U^{1/w}$ with $U \sim \mathrm{Uniform}(0,1)$ and zone weight $w = 1$
(mucus) or $w =$ `LUMEN_COMPETITION_CULL_WEIGHT` $(= 3.0)$ (lumen); the cells with the largest keys
are removed. With equal weights this reduces exactly to uniform random culling, preserving legacy
behaviour, while the up-weighted lumen is thinned first. The zone weights are project-specific
mechanics rather than empirical competition coefficients.<sup>[48](#ref-48)</sup>

Together, slower mucus-associated motility, lower modeled washout, reduced mucus drop deposition,
and weighted carrying-capacity culling are implemented mechanisms intended to create a candidate
refuge.<sup>[48](#ref-48)</sup> Whether these mechanics reduce low-abundance extinction artifacts is
an explicit proposal hypothesis, not a demonstrated result. It must be tested against a homogeneous
control across fixed seed sets and predefined persistence and extinction metrics.

**Zone occupancy readout.** Every update the simulation classifies each live cell by its radial
fraction and records the lumen/mucus population split into the statistics arrays `lumen` and
`mucus` (which always sum to `total`). These are surfaced live in the **Stats** tab (Lumen Count /
Mucus Count) and exported as the `lumen` and `mucus` columns of the Global Stats CSV section, so
the compartmental dynamics can be analysed offline without changing any simulation behaviour.


### 4.21 Reproducible Baselines, Time-to-Extinction, and Survival Curves (Phase I)

To compare dietary and spatial conditions quantitatively, runs must be reproducible and extinction
events must be measurable.

**Seeded determinism.** The `Simulation` accepts an optional integer `seed`. When provided, the
RNG is seeded *before* any stochastic draw (including the random initial positions and sizes of
every cell), so two runs with the same seed, species, and parameters produce byte-identical
population trajectories. In the application the seed is read from `env_config.json` and applied in
`init_sim` prior to species construction; resetting the simulation re-applies the same seed.
Leaving the seed `null` reproduces the prior nondeterministic behaviour.

**Time-to-extinction.** Each species carries an `extinction_time`. Because a population cannot
regrow once it reaches zero in this model, the metric is the elapsed simulation time at the first
zero crossing. Let $N_s(t)$ be the count of species $s$ at elapsed time $t$; then

$$
T_s^{\text{ext}} = \min\lbrace\thinspace t : N_s(t) = 0 \thinspace\rbrace
$$

recorded once and never overwritten. Survivors retain $T_s^{\text{ext}} = \varnothing$ (blank in
the export). Elapsed time is the real-time accumulator $t_{\text{elapsed},n} = \sum_{i=0}^{n}\Delta t_i$,
so the metric is frame-rate independent.

**Export.** `simulation_stats.csv` gains an `--- Extinction Events ---` section with one row per
species: `Species`, `Time_To_Extinction_s`, `Extinct` (0/1), and `Final_Count`. The run's `seed`
is also written into the environment-configuration section so a baseline can be exactly reproduced
from its archived CSV.

**Survival curve.** From the per-species extinction times we build a discrete Kaplan-Meier-style
survival function<sup>[16](#ref-16)</sup> $S(t)$ — the fraction of the $M$ tracked species still extant at time $t$:

$$
S(t) = \frac{1}{M}\thinspace\Bigl\lvert \lbrace\thinspace s : T_s^{\text{ext}} \gt t \ \text{or}\ T_s^{\text{ext}} = \varnothing \thinspace\rbrace \Bigr\rvert
$$

a right-continuous step function that drops by $1/M$ at each extinction event.

**Plot deliverable.** The standalone script `plot_baselines.py` (repository root, depending only on
pandas and matplotlib) parses an exported CSV and writes up to four figures: `abundance_trajectories.png`
(each species' population versus elapsed time), `extinction_survival.png` (the survival curve
$S(t)$ alongside a per-species time-to-extinction bar chart), `zone_occupancy.png` (lumen vs mucus
population over time, see Section 4.20), and `substrate_consumption.png` (cumulative per-substrate
mass ingested, see Section 4.22). The last two are skipped automatically when the CSV predates those
columns. It imports none of the simulation code, so baselines can be regenerated from any archived
run via `python plot_baselines.py --input simulation_stats.csv --outdir .`.


### 4.22 Per-Substrate Consumption Analytics (Phase I)

To quantify *how a diet is actually metabolised* (not merely supplied), the simulation tracks the
mass of each substrate class ingested by the community over time. This closes the loop on the diet
model of Section 4.19: `diet_profile` controls the substrate **supply** mix, while these statistics
record the realised **uptake**, which depends on the resident taxa's affinity vectors.

**Per-cell attribution.** During the affinity-gated feeding step (Section 4.19), each ingestion
event of mass $\Delta m$ from a patch of class $s$ is attributed to that class in a per-cell ledger
that is reset at the start of every cell's update. Thus each cell reports only the mass it ingested
on the current tick, indexed by substrate id $s \in \lbrace 0, 1, 2 \rbrace$.

**Community accumulation.** Each update the simulation sums the per-cell ledgers into a per-step
total $\Delta M_s$ and carries forward the previous cumulative value, so the recorded series is the
total mass ever ingested per class:

$$
M_s(t_n) = M_s(t_{n-1}) + \Delta M_s(t_n) \qquad M_s(t_0) = \Delta M_s(t_0).
$$

Because ingested mass is never returned, each $M_s$ is monotonic non-decreasing. These cumulative
totals are surfaced live in the **Stats** tab (Calorie / Glycan / MDCF Consumed, in pg) and exported
as the `consumed_calorie`, `consumed_glycan`, and `consumed_mdcf` columns of the Global Stats CSV
section (Section 7). The tracking is observability-only and changes no simulation behaviour;
legacy CSVs that predate it import with blank consumption columns.

Comparing $M_s$ trajectories across diet arms reveals, for example, whether an MDCF-enriched diet is
being taken up (high $M_2$ slope) only when MDCF-competent taxa are present — the core selective
mechanism of the malnutrition model.


### 4.23 Initial Inoculum Spawning

Real microbial communities are not seeded as a uniform sheet of cells. An inoculum arrives as a
small number of founder cells that survive transit, attach, and establish localised
**microcolonies** (clonal patches) which then expand clonally by binary fission — and because
daughter cells are placed at the parent (Section 4.7), that clonal structure is preserved as the
population grows. The legacy initialiser instead scattered every cell uniformly in a fixed
$800 \times 600$ rectangle anchored at the top-left corner, which neither matched the circular
petri-dish geometry (Section 4.20) nor the lumen where fresh substrate is delivered. The spawning
model replaces that rectangle with a configurable, dish-aware, seed-reproducible founder placement<sup>[38](#ref-38)</sup>.

**Deferred, geometry-aware placement.** The petri dish centre $(c_x, c_y)$ and radius $R$ are only
known once the window size is available, so founder placement is deferred to the first simulation
update via `Simulation._seed_initial_population` (guarded by a one-shot flag). The dish centre is
anchored at the **world origin** $(c_x, c_y) = (0, 0)$ — the same pose the grid camera defaults to and
that the **Home** button restores (centre $(0,0)$ at $1\times$ zoom) — so the freshly seeded inoculum
is always framed without any camera offset. Each species stores its requested founder count
`init_count` and is placed once geometry exists; `reset` re-applies the same policy. Every founder is
constrained to lie within the dish, $\sqrt{(x-c_x)^2+(y-c_y)^2}\le R$.

**Area-uniform disc sampling.** Where a founder is drawn uniformly over an annular band
$[r_{\min}, r_{\max}]$, the radius uses

$$
r = \sqrt{U\thinspace(r_{\max}^2 - r_{\min}^2) + r_{\min}^2} \qquad U \sim \mathrm{Uniform}(0,1)
$$

so density is constant per unit **area** rather than per unit radius (a naive $r = U\thinspace R$ draw would
pile founders near the centre). The angle is uniform over the species' allowed sector.

**Lumen bias.** Fresh dietary substrate is delivered chiefly to the well-mixed lumen (Section 4.20),
so a gut-realistic inoculum is biased inward. Each founder is placed in the inner lumen band
($r \in [0, f R]$, with $f =$ `INITIAL_SPAWN_LUMEN_RADIUS_FRACTION`) with probability
`INITIAL_SPAWN_LUMEN_BIAS` (default $0.9$), otherwise in the outer mucus annulus ($r \in [fR, R]$).

**Spawn modes** (`INITIAL_SPAWN_MODE`):

- **`clustered`** (default): generate `INITIAL_SPAWN_CLUSTER_COUNT` founder nuclei (each lumen-biased),
  then scatter cells in Gaussian microcolonies about them. A cell offset from nucleus
  $(n_x, n_y)$ is $x = n_x + \mathcal{N}(0, \sigma^2)$, $y = n_y + \mathcal{N}(0, \sigma^2)$ with
  $\sigma =$ `INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION` $\times R$; offsets are resampled (up to a few
  attempts) so no cell escapes the dish. Founders are assigned round-robin across nuclei for an even,
  deterministic colony-size split. This seeds the spatial structure on which competition,
  cross-feeding, and HGT act.
- **`uniform_disc`**: area-uniform placement across the dish, still lumen-biased per cell.
- **`central`**: a single compact droplet within `INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION` $\times R$ of
  the centre (a point-inoculation assay).

**Inter-species arrangement** (`INITIAL_SPECIES_ARRANGEMENT`):

- **`independent`** (default): each species draws its own founder nuclei; colonies overlap only by
  chance — the generic well-mixed case.
- **`mixed`**: all species share one set of founder nuclei, co-locating colonies to maximise early
  contact-dependent interactions (T6SS, bacteriocins, HGT, cross-feeding).
- **`segregated`**: species $i$ of $n$ is confined to the angular sector
  $[\thinspace 2\pi i/n, 2\pi (i{+}1)/n\thinspace)$, separating founder origins to study invasion / range-expansion
  dynamics.

**Reproducibility.** All draws come from the shared RNG stream, which is seeded before placement
(Section 4.21), so a given seed regenerates a byte-identical starting layout. The relevant constants
are `INITIAL_SPAWN_MODE`, `INITIAL_SPAWN_LUMEN_BIAS`, `INITIAL_SPAWN_LUMEN_RADIUS_FRACTION`,
`INITIAL_SPAWN_CLUSTER_COUNT`, `INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION`,
`INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION`, and `INITIAL_SPECIES_ARRANGEMENT` (Section 4.13).

**Live configuration (Globals tab).** These spawn parameters are exposed as editable environment
keys in the Globals tab so an experiment can be re-framed without code changes. The categorical keys
`spawn_mode` (cycles `clustered` / `uniform_disc` / `central`) and `species_arrangement` (cycles
`independent` / `mixed` / `segregated`) advance through their option lists on click/Enter, while the
numeric keys `spawn_lumen_bias`, `spawn_cluster_count`, `spawn_cluster_spread`, and
`spawn_central_radius` are text-edited with validation (fractions clamped to $[0,1]$; the cluster
count is an integer $\ge 1$). Edits are persisted to `env_config.json` and applied by
`Simulation._apply_spawn_env` at the start of the next reset/re-seed (they do not relocate an
already-placed population mid-run). Blank or unrecognised entries fall back to the compiled-in
constant defaults, so the built-in behaviour is reproduced when nothing is overridden.


## 5. Species Parameters

Species are configured in `species_data.csv` and can also be edited in the Species tab.

* `name`: Species identifier displayed in the UI and statistics.
* `init`: Initial population count spawned at reset/start.
* `rep`: Base reproduction probability per second before carrying-capacity and environmental modifiers.
* `aging`: Base intrinsic aging multiplier used in the Gompertz death penalty.
* `so`: Oxygen sensitivity from -1 to 1. Negative values model anaerobes harmed by oxygen; positive values model aerobes needing oxygen.
* `sf`: Food sensitivity from 0 to 1. This is retained as a species-level nutrient responsiveness parameter for UI/data compatibility.
* `speed`: Base movement speed in micrometers per second before inverse-size drag scaling.
* `ph`: Preferred pH from 0 to 14. Survival fitness follows a Gaussian curve around this value.
* `food_capacity`: Base nutrient storage capacity in pg glucose-equivalent. Actual capacity scales by `base_size^3`.
* `metabolic_rate`: Basal nutrient use in pg/s. Actual metabolism scales by `base_size^2.25` and age.
* `consumption_speed`: Maximum nutrient uptake in pg/s. Actual uptake scales by `base_size^2`.
* `greediness`: Chemotaxis responsiveness from 0 to 1. Hunger raises the effective value toward 1.
* `t6ss_activity`: Contact-dependent Type VI secretion activity from 0 to 1. Damage scales with attacker surface area and elapsed time.
* `bacteriocin_production`: Diffusible antimicrobial production from 0 to 1. Production costs nutrient and emits bacteriocin into the ECM.
* `immunity_factor`: Resistance from 0 to 1 against bacteriocin exposure and T6SS contact damage.
* `rep_thresh`: Required fullness fraction from 0 to 1 before fission can occur.
* `share_thresh`: Fullness fraction from 0 to 1 above which same-species nutrient sharing may occur.
* `predation_thresh`: Fullness fraction from 0 to 1 below which larger cells may consume smaller competing cells.
* `base_size`: Structural cell size in micrometers. It controls capacity, metabolism, uptake, movement speed, turning, and necromass.
* `sensor_radius`: Food sensing radius in micrometers. Food outside this range does not contribute to chemotaxis or local feeding checks.
* `temp_min`: Ratkowsky minimum cardinal temperature $T_{\min}$ in degrees Celsius. No growth occurs at or below it (see Section 4.14)<sup>[23](#ref-23)</sup>.
* `temp_max`: Ratkowsky maximum cardinal temperature $T_{\max}$ in degrees Celsius. No growth occurs at or above it. It is clamped to be strictly above `temp_min`.
* `acetate_secretion`: Overflow acetate excreted per unit of fresh glucose uptake, from 0 to 1. Producers emit acetate into the ECM (see Section 4.15). Default 0 (no overflow)<sup>[30](#ref-30)</sup>.
* `acetate_uptake`: Acetate cross-feeding uptake activity from 0 to 1. Consumers oxidize local acetate as a secondary aerobic carbon source (see Section 4.15). Default 0 (no cross-feeding).
* `baranyi_q0`: Initial Baranyi physiological state $q_0$ of the fresh inoculum (dimensionless, $\gt 0$). A smaller value lengthens the lag phase via $h_0 = \ln(1 + 1/q_0)$ (see Section 4.16); it is clamped to $[10^{-4}, 10^3]$. Default 0.05 (moderate lag)<sup>[2](#ref-2)</sup>.
* `persister_switch`: Persister switching propensity from 0 to 1 (see Section 4.17). It scales the entry rate into the dormant persister state; 0 fully disables persistence (pre-M4 behaviour). Default 0<sup>[1](#ref-1), [20](#ref-20)</sup>.
* `persister_enter_thresh`: Starvation fullness threshold $S_1$ from 0 to 1. Below it, starvation-induced switching into persistence engages. Default 0.15.
* `persister_exit_thresh`: Recovery fullness threshold $S_2$ from 0 to 1. Above it, persisters resuscitate to the proliferative state. It is clamped to stay at least `PERSISTER_THRESHOLD_MIN_GAP` above $S_1$. Default 0.5.
* `phage_resistance`: Initial fraction of cells born resistant to bacteriophage adsorption, from 0 to 1 (see Section 4.18). Resistance is heritable; resistant cells survive an epidemic and drive the rebound. Default 0 (fully susceptible)<sup>[21](#ref-21)</sup>.
* `hgt_competence`: Horizontal gene transfer competence from 0 to 1 (see Section 4.18). Gates both donating and receiving phage resistance on cell contact; the recipient inherits a metabolic carriage cost. 0 disables HGT (pre-M5 behaviour). Default 0<sup>[22](#ref-22), [26](#ref-26)</sup>.
* `aff_calorie`: Nonnegative, dimensionless generic-calorie throughput/capacity multiplier (see Sections 4.19 and 10.4). Default 1.0 preserves full legacy pg uptake. It is not a measured molecular binding affinity.<sup>[47](#ref-47), [48](#ref-48)</sup>
* `aff_glycan`: Nonnegative, dimensionless complex-glycan throughput/capacity multiplier. A value of 0 disables uptake in both diet paths. Default 0; the production pg path's $0.7$ energy yield is a separate Own-model calibration.<sup>[47](#ref-47), [48](#ref-48)</sup>
* `aff_mdcf`: Nonnegative, dimensionless MDCF-prebiotic throughput/capacity multiplier. A value of 0 disables uptake in both diet paths. Default 0; the production pg path's $0.85$ energy yield is a separate Own-model calibration.<sup>[47](#ref-47), [48](#ref-48)</sup>

### 5.1 Parameter Inventory and Missing Values

`SpeciesDataManager.fieldnames` currently contains 33 fields: `name`, `init`, and 31
traits. These are 32 numeric fields, not 35 required numeric measurements. The loader
retains backward-compatible defaults for missing/empty cells and the legacy
`aggression`/`defense` aliases. The species constructor also normalizes selected
values. A successful load does not establish that a parameter was measured or is
biologically appropriate.<sup>[51](#ref-51)</sup>

Run a read-only inventory before interpreting a species configuration:

```bash
python audit_species_parameters.py --species-data species_data.csv --temperature 37 --output parameter_audit.json
python -m unittest bacteria.test_parameter_audit
```

The schema-1 JSON records the input SHA-256, field/profile/founder counts, raw and
loader-resolved values, explicit/default/alias origins, numeric/structural issues,
constructor transformations, and each profile's implemented temperature multiplier.
An explicit CSV value is not a measured value: per-value empirical provenance is
reported as unrecorded. `--require-explicit` makes compatibility-default use fail
the audit without changing the loader. Exit status is 0 for a completed audit
without detected issues, 1 for reported issues or the requested explicit-value
gate, and 2 for input/output errors. Existing output files are never overwritten.
The audit also checks the interactive app's integer-text parsing of `init`: values
such as `2.9`, `100.0`, or `1e2` can be converted by the sweep/constructor path but
are rejected by the app, which skips that profile. These inputs produce audit
issues and exit status 1. Constructor founder counts remain diagnostic values, not
executed populations; `interactive_app_init_count` records the app-compatible count
or null. `csv_row` numbers logical records including the header, while
`csv_line_end` locates the physical end of each record, including blank lines and
quoted multiline cells.<sup>[51](#ref-51)</sup>
This is a configuration audit, not an empirical imputation or calibration tool.
<sup>[51](#ref-51)</sup>

A constant environment temperature does not force the implemented species response
to one. Cardinal-temperature parameters still affect its value. Inspect that value
before proposing a sensitivity analysis or removing parameters; the audit does not
change the response model described in Section 4.14.<sup>[51](#ref-51)</sup>

## 6. Global Environment Parameters

Globals are stored in `env_config.json` and can be edited in the Globals tab.

* `o2`: Dissolved oxygen fraction from 0 to 1. Aerobes benefit from oxygen; anaerobes are harmed by it.
* `co2`: Dissolved carbon dioxide fraction from 0 to 1. Respiration increases it and buffering/pH logic responds to it.
* `ph`: Current environmental pH from 0 to 14.
* `temperature`: Environment temperature in degrees Celsius (slider range 0 to 60). It drives the Ratkowsky growth multiplier $f_{T}$ applied per species (see Section 4.14)<sup>[23](#ref-23)</sup>.
* `capacity`: Hard carrying capacity across all species. Excess population is culled if total cells exceed this value.
* `drop_frequency`: Expected food clusters per real second. The implementation uses an expected-rate spawning approximation $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$.
* `drop_range_min` and `drop_range_max`: Uniform random range for the amount in each new food cluster.
* `phage_inflow`: Expected free-phage seed packets introduced per real second (see Section 4.18). Each packet deposits `PHAGE_SEED_PACKET_VIRIONS` virions at a random petri location via the same Poisson rate $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$. Set above 0 to trigger an epidemic and back to 0 to halt seeding. Default 0 (no outbreak)<sup>[8](#ref-8), [19](#ref-19)</sup>.
* `diet_profile`: Active Own-model categorical nutrient-drop scenario (see Section 4.19), one of `malnutrition`, `standard`, or `mdcf`. It sets the production pg path's probability mix of generic calorie, complex glycan, and MDCF prebiotic classes; it does not select a Section 10.4 dry-recipe profile.<sup>[48](#ref-48)</sup> It is editable directly in the **Globals** tab — clicking its cell (or pressing Enter on the row) cycles to the next scenario. Like all Globals-tab values it defines starting conditions, so a change takes effect on the next simulation reset/reload. Default `standard`.
* `seed`: Optional integer RNG seed for reproducible extinction/abundance baselines (see Section 4.21)<sup>[16](#ref-16)</sup>. It is editable directly in the **Globals** tab — enter a non-negative integer to fix the run, or leave the field blank (shown as `(random)`) for a fresh nondeterministic run. When set, the random stream is seeded before species construction so the entire run — initial layout and every stochastic draw — is deterministic. Like all Globals-tab values it defines starting conditions, so a change takes effect on the next simulation reset/reload.
* `spawn_mode`: Initial inoculum geometry (see Section 4.23), one of `clustered`, `uniform_disc`, or `central`. Categorical — clicking its cell (or pressing Enter on the row) cycles to the next mode. Default `clustered`.
* `spawn_lumen_bias`: Probability (0 to 1) that a founder cell is placed in the inner lumen band rather than the outer mucus annulus. Default 0.9.
* `spawn_cluster_count`: Number of founder microcolonies per species in `clustered` mode (integer $\ge 1$). Default 3.
* `spawn_cluster_spread`: Microcolony Gaussian spread as a fraction (0 to 1) of the dish radius. Default 0.06.
* `spawn_central_radius`: Droplet radius as a fraction (0 to 1) of the dish radius in `central` mode. Default 0.15.
* `species_arrangement`: Inter-species layout (see Section 4.23), one of `independent`, `mixed`, or `segregated`. Categorical — clicking its cell (or pressing Enter on the row) cycles to the next option. Default `independent`. All six spawn keys define starting conditions, so changes take effect on the next simulation reset/reload.

## 7. Data Export And Import

`simulation_stats.csv` stores these logical sections when exported from the app:

* Environment configuration (now including the reproducibility `seed`, see Section 4.21).
* Species configuration.
* Global and per-species time-series statistics.
* Extinction events (Phase I): one row per species with header `Species,Time_To_Extinction_s,Extinct,Final_Count`, where `Time_To_Extinction_s` is the elapsed sim-time at the first zero-population crossing (blank for survivors) and `Extinct` is 0/1. This section feeds the survival-curve plot in `plot_baselines.py`.

The Global Stats header is:

```csv
time,total,births,deaths,elapsed_seconds,dt_seconds,tickrate_hz,overall_births,overall_deaths,t6ss_deaths,hgt_transfers,overall_t6ss_deaths,overall_hgt_transfers,lumen,mucus,consumed_calorie,consumed_glycan,consumed_mdcf,fba_biomass_gdw,fba_acetate_mol,fba_consumed_calorie_mol,fba_consumed_glycan_mol,fba_consumed_mdcf_mol
```

Column meanings:

* `time`: integer simulation step index retained for backward compatibility.
* `total`: total live bacteria after that update.
* `births`: births during that update only.
* `deaths`: deaths during that update only, including capacity culls.
* `elapsed_seconds`: cumulative real simulated seconds.
* `dt_seconds`: elapsed seconds represented by this update.
* `tickrate_hz`: instantaneous update frequency for this row.
* `overall_births`: cumulative global births through this row.
* `overall_deaths`: cumulative global deaths through this row.
* `t6ss_deaths`: T6SS positive-to-nonpositive food transitions during this update (see Section 4.12).
* `hgt_transfers`: false-to-true phage-resistance acquisitions through HGT during this update (see Section 4.18).
* `overall_t6ss_deaths`: cumulative `t6ss_deaths` through this row.
* `overall_hgt_transfers`: cumulative `hgt_transfers` through this row.
* `lumen`: live bacteria in the central lumen (radial fraction below the mucus threshold) after that update (see Section 4.20).
* `mucus`: live bacteria in the outer mucus annulus (radial fraction at or above the mucus threshold) after that update. `lumen` and `mucus` always sum to `total`. These columns are blank for stats produced by legacy imports that predate zone tracking.
* `consumed_calorie`: cumulative generic-calorie substrate mass (pg) ingested by the whole community through this row (see Section 4.22).
* `consumed_glycan`: cumulative complex-glycan substrate mass (pg) ingested through this row.
* `consumed_mdcf`: cumulative MDCF-prebiotic substrate mass (pg) ingested through this row. Each `consumed_*` series is monotonic non-decreasing; the three columns are blank for stats produced by legacy imports that predate diet analytics.
* `fba_biomass_gdw`: total committed sidecar biomass of the live cells after that live-FBA update, in gDW.
* `fba_acetate_mol`: acetate bookkeeping output reported for that live-FBA update, in moles.
* `fba_consumed_calorie_mol`, `fba_consumed_glycan_mol`, and `fba_consumed_mdcf_mol`: represented substrate consumed during that live-FBA update, in moles; these are per-tick values, not cumulative pg-path totals.<sup>[41](#ref-41)</sup>

Legacy-path exports leave the five FBA-only columns blank. Import accepts files that omit them and
initializes each absent FBA series to zero.<sup>[41](#ref-41)</sup>

The time fields are computed as:

$$
t_{\text{elapsed},n}=\sum_{i=0}^{n}\Delta t_i
$$

$$
f_{\text{tick},n}=\frac{1}{\Delta t_n}
$$

The cumulative global events are running sums:

$$
B_{\text{overall},n}=B_{\text{overall},n-1}+B_n
$$

$$
D_{\text{overall},n}=D_{\text{overall},n-1}+D_n
$$

Import remains backward-compatible with older four-column global stats files. If a legacy row only has `time,total,births,deaths`, the importer treats `time` as a completed update index and infers:

$$
\Delta t = \frac{1}{\text{FPS}} \quad t_{\text{elapsed}} = \frac{\text{time}+1}{\text{FPS}} \quad f_{\text{tick}}=\text{FPS}
$$

and reconstructs cumulative global births and deaths by running sum. Legacy files that predate all
four T6SS/HGT event columns import zeros for all four series.

Per-species rows include species index, name, step, population, step births, step deaths, cumulative births, and cumulative deaths.

## 8. User Interface

The application is a real-time Pygame interface with five main tabs:

* Grid: Primary simulation view with panning, zooming, minimap navigation, and a scale indicator showing the current pixel-to-micrometer conversion.
* Graph: Population history charts with export/import controls for timeline stats.
* Stats: Live numerical population, step births/deaths, cumulative global births/deaths, elapsed seconds, tickrate, and per-species summaries.
* Species: Editable table for species parameters, including duplicate/delete and undo/redo support.
* Globals: Editable environment table for oxygen, carbon dioxide, pH, capacity, and food drop settings.

Environment sliders in the Grid, Graph, and Stats views allow live adjustments to O2, CO2, pH, food drop frequency, and food drop range.

## 9. Software Architecture

`Zach Summer 2026 Project Proposal (Grounded Version).md` is the scientific charter: it governs the
research questions, hypotheses, scope, timeline, and intended deliverables. The software inventory
below describes implemented components; it does not redefine those proposal phases.

* `bacteria/app.py`: Pygame application loop, rendering, tabs, UI state, and high-level simulation orchestration.
* `bacteria/core_sim.py`: Headless biological and physical simulation models: `FoodNode`, `Bacterium`, `BacteriaSpecies`, `FoodSpatialHash`, and `Simulation`, including the opt-in live-FBA transaction boundary.
* `bacteria/ecm.py`: Public sparse-grid ECM API and Python/native backend orchestration for the four signal layers, plus separate Python-transported continuous substrate maps for the opt-in live-FBA path.
* `bacteria/physics.py`: Reusable structure-of-arrays buffers, the reference contact implementation, and Python/native physics backend orchestration.
* `bacteria/metabolism.py`: Versioned reduced-model and prior loading, exact immutable NumPy model buffers, the Python reference solver/cache/adaptive batch implementation, and lazy `MetabolicOrchestrator` modes and backends.
* `bacteria/metabolic_integration.py`: Opt-in spatial substrate registration, diet-bolus deposition, deterministic live-agent batching, reduced-metabolism state ownership, lifecycle policy, represented-amount validation, and adapter-level transactions used within the complete `Simulation` transaction.
* `bacteria/diet_translator.py`: Strict loader and validator for the versioned diet-profile asset plus the pure dry-bolus-to-three-substrate molar translation API.
* `bacteria/models/diet_profiles.json`: Versioned dry-recipe, measured-component, and Own-model routing/calibration asset used only by `DietTranslationEngine` (Section 10.4).
* `bacteria/models/reduced_metabolism_v1.json` and `bacteria/models/metabolic_priors_v1.json`: Revision-1 reduced metabolic bookkeeping model and neutral, nonempirical prior rows.
* `bacteria/src/ecm_kernel.cpp`, `bacteria/src/ecm_kernel.h`, `bacteria/src/physics_kernel.cpp`, `bacteria/src/physics_kernel.h`, `bacteria/src/metabolic_kernel.cpp`, and `bacteria/src/metabolic_kernel.h`: Optional C++17 ECM, contact-physics, and reduced-metabolism implementations compiled into the unified top-level pybind11 module `bacteria_c_kernel`; the sole `PYBIND11_MODULE` definition remains in `ecm_kernel.cpp`.
* `bacteria/data.py`: CSV/JSON serialization for `species_data.csv` and `env_config.json`.
* `bacteria/ui_components.py`: Reusable buttons, sliders, and editable cells.
* `bacteria/constants.py`: UI constants plus documented biological and numerical model constants.
* `bacteria/__main__.py`: CLI entry point for `python -m bacteria`.
* `setup.py` and `pyproject.toml`: Unified ECM/physics/metabolism native-extension build configuration, all three C++ source files, model-asset packaging, and Python package metadata.
* `bacteria/test_kernel.py`, `bacteria/test_physics_kernel.py`, `bacteria/test_metabolic_architecture.py`, and `bacteria/test_metabolic_kernel.py`: ECM, contact, and reduced-metabolism architecture, parity, solver, backend, ABI, boundary, transaction, cache, scale, stress, and opt-in performance checks.
* `bacteria/test_diet_translator.py`: Focused diet-asset schema, translation formula, array contract, reduced-orchestrator injection, parity, rejection, and smoke coverage.
* `bacteria/test_metabolic_integration.py` and `bacteria/test_integration_stress.py`: Live-FBA coordinate, bolus, batch, lifecycle, conservation, full-transaction rollback, scale, ownership, and conditional Python/native parity coverage.
* `run_comparative_sweeps.py` and `bacteria/test_comparative_sweeps.py`: Paired-seed legacy-pg benchmark harness, additive run manifest, and spatial-control, censoring, configuration, schema, affinity-audit, RNG-isolation, and non-overwrite tests.
* `audit_species_parameters.py` and `bacteria/test_parameter_audit.py`: Read-only raw/default/alias parameter inventory and constructor/temperature diagnostics; no empirical parameter inference.
* `bacteria/test_bacteria.py`: Legacy production-simulation regression suite, including core biology, categorical pg diet drops, lumen/mucus mechanics, serialization, and contact behavior.
* `plot_baselines.py`: Standalone parser and plotting utility for exported `simulation_stats.csv`, producing abundance, extinction/survival, zone-occupancy, and substrate-consumption figures when the corresponding data are present.
* `.github/workflows/kernel.yml`: `actions/checkout@v7.0.1` and `actions/setup-python@v7.0.0` read the shared Python version from `.python-version` and define OpenMP Windows/Linux/macOS and serial Linux build/test jobs. A separate requirements-only Linux job verifies `pygame-ce`, absence of the native extension, and the full suite with native-only skips. Both build paths assert that the unified extension exposes `MetabolicKernel` and invoke the ECM, physics, reduced-metabolism, translator, live-integration, integration-stress, comparative-sweep, and legacy suites; the serial job also invokes the required metabolic and physics scale fixtures. The file separately defines a serial Linux Valgrind job. These are workflow contracts, not evidence that a hosted workflow or Valgrind run passed.

## 10. Hybrid C++ Computational Architecture, Empirical Parameterization, and Validation

This section is organized by implemented component, not by proposal phase. The C++ kernels are
optional accelerators and are not proposal Phase II; optimization remains secondary to biological
correctness, reproducibility, and interpretation. The production Pygame simulation retains the
categorical pg glucose-equivalent metabolism and nutrient-drop path described in Section 4.19. The
reduced molar translator and LP orchestrator remain opt-in: a caller may attach the live adapter to
`Simulation`, but the application has no user-facing FBA configuration and the default constructor
continues to select the legacy pg/`FoodNode` path.

Cross-cutting implementation notes:

* Food nodes and bacteria use small object pools to reduce allocation churn.
* Food lookup uses `FoodSpatialHash` so bacteria query nearby food instead of scanning the full node list.
* Bacterial physical interactions use a spatial hash over neighboring cells.
* ECM diffusion updates every active cell each tick, but only stores non-negligible concentrations.
* Viewport culling suppresses expensive local chemotaxis for off-viewport bacteria while retaining simple movement.
* The optional reduced metabolic API batches solver work, deduplicates identical misses, and uses a bounded per-kernel cache without changing the app's default metabolism.

### 10.1 Optional ECM Diffusion Kernel

The ECM can optionally run diffusion and decay through the top-level pybind11 extension `bacteria_c_kernel`. The unified extension is compiled from `bacteria/src/ecm_kernel.cpp`, `bacteria/src/physics_kernel.cpp`, and `bacteria/src/metabolic_kernel.cpp`; the existing `process_layers`, `process_interactions`, and `openmp_enabled` APIs remain available alongside the `MetabolicKernel` class described in Section 10.3. This is an acceleration backend, not a separate ECM model: `ECM(cell_size=..., backend='auto')` retains the original public layer dictionaries and the `add`, `sample`, and `tick` methods. After a native update, each public dictionary is cleared and updated in place so callers that retain a reference continue to observe the same object.

Backend selection has three modes:

* `auto` performs no native import during construction or empty ticks. On the first nonempty tick it lazily imports the extension, falling back to Python only when the module is absent or lacks the required `process_layers` symbol; native computation errors propagate instead of silently changing backends.
* `python` always uses the reference implementation and never imports the extension.
* `cpp` requires both the extension and `process_layers`, and reports either absence rather than falling back.

All sparse ECM fields use one world-to-cell convention. For finite world coordinates $x$ and $y$
in micrometers and ECM cell width $h$ in micrometers, Python applies mathematical floor on each
axis, including for negative coordinates:<sup>[39](#ref-39)</sup>

$$
(c_x,c_y)=\left(
\left\lfloor\frac{x}{h}\right\rfloor,
\left\lfloor\frac{y}{h}\right\rfloor
\right).
$$

Thus, with $h=10\thinspace\mu\mathrm{m}$, $x=-0.1\thinspace\mu\mathrm{m}$ maps to $c_x=-1$,
not zero. The four signal dictionaries and the three continuous FBA substrate maps share this
integer coordinate system.

The pybind11 binding exposes this fixed interface:

```text
process_layers(coordinates, values, offsets, diffusion, decay, dt,
         cell_size, alpha_max, min_concentration)
  -> (coordinates, values, offsets)
```

All array buffers use C-contiguous storage:

| Buffer | Dtype | Shape | Contract |
|---|---|---|---|
| `coordinates` | signed `int64` | `(N, 2)` | Open-grid integer coordinate pairs. |
| `values` | `float64` | `(N,)` | Concentrations corresponding one-to-one with the coordinate rows. |
| `offsets` | `int64` | `(5,)` | Monotonic boundaries from `0` through `N` for the four layers. |
| `diffusion` | `float64` | `(4,)` | One diffusion coefficient per layer. |
| `decay` | `float64` | `(4,)` | One decay coefficient per layer. |

The fixed layer order is `qs`, `bac`, `acetate`, then `phage`. The return tuple uses the same coordinate/value/offset layout, with `N` denoting the returned entry count. Python flattens only active dictionary entries into these buffers and reconstructs only retained entries, avoiding dense bounding boxes on the unbounded grid. The flat sparse-buffer ABI, lazy backend orchestration, and four-layer parallel dispatch are original BacSandbox interface designs.<sup>[39](#ref-39)</sup>

This native ABI remains signal-only: it always has the four layers above and an offsets array of
shape `(5,)`. The separately registered `generic_calorie`, `complex_glycan`, and
`mdcf_prebiotic` concentration maps are not inserted into `_layers`, are not flattened into
`process_layers`, and are not passed to the native ECM kernel. `ECM.tick_substrates` currently
advances each of those three maps with the Python Fickian reference stencil<sup>[9](#ref-9)</sup>
and the same stability-clamp principle<sup>[5](#ref-5)</sup>, regardless of the signal-layer
backend selection. There is therefore no implemented C++ substrate-diffusion path or associated
performance claim.<sup>[39](#ref-39)</sup>

While holding the GIL, the binding validates buffer layouts, shapes, offsets, and scalar parameters, then copies inputs into RAII-owned native vectors. It releases the GIL around native dispatch and reacquires it before constructing the NumPy outputs. The native and OpenMP regions make no Python/C API calls. RAII provides explicit ownership discipline; neither that style nor an unexecuted Valgrind workflow definition is proof of leak freedom.

The native path preserves the §4.9 Fickian stencil<sup>[9](#ref-9)</sup> and stability clamp<sup>[5](#ref-5)</sup> exactly. Python supplies every model constant, so the C++ source duplicates none. For each layer, native code applies the alpha clamp and linear decay clamp; visits center, left, right, down, and up contributions in that order; applies the per-contribution floor; aggregates equal coordinates; and performs final concentration pruning. It preserves open signed coordinates and deterministic first-insertion output order, and guards neighbor calculations against signed `int64` overflow.

When OpenMP is compiled in, a guarded `parallel for` with `schedule(static)` distributes the four independent layer updates. Each worker owns its vectors and hash table, captures failures in an `exception_ptr`, and rethrows after the parallel region. A serial build executes the identical dispatch path without the pragma.

The root `pyproject.toml` and `setup.py` support Python `>=3.14,<3.15`, NumPy 2.3.4+/pybind11 3.0.1+, C++17, and portable compiler optimization flags without CPU-specific options. `BACSANDBOX_OPENMP=auto` probes OpenMP and builds serially when it is unavailable; `on` requires OpenMP and fails clearly if probing fails; `off` deliberately builds the serial path. Windows uses MSVC, Linux supports GCC or Clang, and macOS uses Apple Clang with `libomp` probing.

Kernel QA in `bacteria/test_kernel.py` uses `unittest`. Its parity fixture runs 100 ticks and, after every tick, compares all four layers by exact key equality and an absolute value tolerance of `1e-5`. Focused tests cover backend behavior, ABI validation, clamping, pruning, coordinate overflow, and allocation behavior. The bounded-fixture performance benchmark reports repeated medians and is opt-in:

```bash
BACSANDBOX_RUN_BENCHMARK=1 python -m unittest bacteria.test_kernel
```

Performance measurements are environment- and workload-specific, not fixed guarantees. The kernel workflow is configured for OpenMP builds on Windows, Linux, and macOS, a Linux serial build, and Linux Valgrind stress checks for definite and indirect losses; no hosted result is asserted here.

### 10.2 Optional Contact-Physics Spatial Hash Kernel

The unified optional `bacteria_c_kernel` module also exposes the top-level `process_interactions`
binding for contact physics. The extension continues to expose the ECM `process_layers` and
`openmp_enabled` APIs and also owns the stateful reduced-model `MetabolicKernel` class. A current build
supports ECM work (Section 10.1), contact-pair processing, and the opt-in reduced metabolic direct
API (Section 10.3) without creating another native module.

**Simulation and backend contract.** `Simulation` accepts the backward-compatible
`physics_backend='auto'|'python'|'cpp'` setting. `auto` imports lazily and falls back only when the
module itself is absent or when an older ECM-only binary lacks `process_interactions`.
`python` never imports the native module, while `cpp` requires the `process_interactions` symbol. Other import
failures and all native runtime or validation errors propagate; they do not silently switch the
model to Python. The orchestrator reports `backend_used` after successful processing and exposes
`native_available` as `None` before a native load is attempted, then as the result of that attempt.

`_process_interactions` is a deterministic two-phase contract. First, Python performs acetate
secretion and uptake (Section 4.15) plus bacteriocin ECM emission and exposure in stable species and
member order. The selected contact backend then performs HGT (Section 4.18), same-species nutrient
sharing, T6SS, and phagocytosis (Section 4.12). This staged contract preserves Python/native parity
for the contact phase, but it intentionally makes no claim of full-tick trajectory identity with
the historical implementation that interleaved all of these operations.

**Flat buffers and ordering.** Flattening visits `Simulation.species` in list order and then each
`sp.members` list in order. `species_id` is the zero-based index in the species list. Reusable
`PhysicsBuffers` storage grows geometrically, and only the active `(N,)` and `(S,)` slices cross the
binding. Every array is one-dimensional and C-contiguous; no `Bacterium` pointer crosses the
Python/C++ boundary. The `event_flags` active slice is reset on every flatten so inactive capacity
tails cannot leak events into a later frame.

Agent buffers:

| Buffer | Dtype | Shape | Source and ownership |
|---|---|---|---|
| `x` | `float64` | `(N,)` | `b.x`; read-only during contact processing. |
| `y` | `float64` | `(N,)` | `b.y`; read-only during contact processing. |
| `interaction_radius` | `float64` | `(N,)` | `b.size`; read-only effective contact radius. |
| `base_radius` | `float64` | `(N,)` | `b.base_size`; read-only reference radius. |
| `food_level` | `float64` | `(N,)` | Mutable in place for sharing, T6SS, and predation. |
| `food_capacity` | `float64` | `(N,)` | Read-only storage capacity. |
| `metabolic_burden_factor` | `float64` | `(N,)` | Mutable in place after successful HGT. |
| `species_id` | `int32` | `(N,)` | Read-only zero-based species-list index. |
| `is_persister` | `uint8` | `(N,)` | Read-only binary status. |
| `phage_resistant` | `uint8` | `(N,)` | Mutable binary HGT donor/recipient status. |
| `event_flags` | `uint8` | `(N,)` | Mutable per-tick event bit mask, reset before dispatch. |

Species buffers:

| Buffer | Dtype | Shape | Ownership |
|---|---|---|---|
| `t6ss_activity` | `float64` | `(S,)` | Read-only. |
| `immunity_factor` | `float64` | `(S,)` | Read-only. |
| `hgt_competence` | `float64` | `(S,)` | Read-only. |
| `share_threshold` | `float64` | `(S,)` | Read-only. |
| `predation_threshold` | `float64` | `(S,)` | Read-only. |

**Binding, validation, and ownership.** The exact binding order is shown below. The arrays are
followed by Python-supplied elapsed time and every model scalar, and `random_draw` is last:

```text
process_interactions(
    x, y, interaction_radius, base_radius, food_level, food_capacity,
    metabolic_burden_factor, species_id, is_persister, phage_resistant,
    event_flags,
    t6ss_activity, immunity_factor, hgt_competence, share_threshold,
    predation_threshold,
    dt, hash_cell_size, axis_cutoff, contact_margin,
    hgt_transfer_rate, hgt_metabolic_burden,
    t6ss_damage_scaling, t6ss_max_immunity_defense,
    t6ss_immunity_size_epsilon, persister_resistance_factor,
    phagocytosis_size_ratio, cell_wet_mass, cell_dry_mass_fraction,
    labile_necromass_fraction, cytosol_recovery_fraction,
    consumed_prey_sentinel, random_draw)
  -> (contact_pairs, t6ss_hits, t6ss_deaths, hgt_transfers,
      share_events, predation_events)
```

`dt` is real elapsed time in seconds. All other scalar values originate in
`bacteria/constants.py` through Python's `ContactParameters`: hash-cell size, axis cutoff and contact
margin, HGT rate and burden, T6SS scaling/immunity/epsilon/persister factors, phagocytosis ratio,
wet/dry/labile/cytosol recovery factors, and the consumed-prey sentinel. C++ duplicates none of
these model constants.

The NumPy array arguments use `noconvert`. With the GIL held, the binding enforces exact dtypes,
one-dimensional shapes, C contiguity, matching agent count `N` and species count `S`, valid species
IDs, binary status bytes, defined event bits, finite values and allowed ranges, writable mutable
arrays, and positive hash-cell size and T6SS size epsilon. Processing then uses zero-copy raw views
over the NumPy storage. Before serial mutation, the binding snapshots all four mutable arrays:
`food_level`, `metabolic_burden_factor`, `phage_resistant`, and `event_flags`. If serial resolution
or the Python RNG callback fails, those snapshots are restored transactionally before the exception
propagates. Python writes object-backed fields back to the `Bacterium` instances only after native
processing and six-integer summary validation both succeed.

**Sparse grid and deterministic discovery.** The flat agent/species interface, CSR-like ordering,
staged OpenMP discovery, and ordered six-counter event contract are original BacSandbox
implementation designs. Coordinates are assigned to signed cells by mathematical floor, not
truncation.<sup>[40](#ref-40)</sup>

$$
c_x = \left\lfloor \frac{x}{h} \right\rfloor
\qquad
c_y = \left\lfloor \frac{y}{h} \right\rfloor
$$

Members are stably sorted by `(cell_x, cell_y, original_index)`. Unique sorted cell keys, offsets,
and member indices form a CSR-like sparse grid, so no dense world-sized allocation is required.
Candidate discovery preserves flattened source order, scans neighbor `nx` in the outer loop and
`ny` in the inner loop, and visits each bucket in original member order.

The existing Section 4.12 contact model first applies inclusive axis gates and then a strict disk
test; the native path deliberately preserves the Python criterion rather than replacing it with a
no-margin or inclusive-tangency test.<sup>[35](#ref-35)</sup>

$$
\lvert \Delta x \rvert \le 15\thinspace\mu\mathrm{m}
\qquad
\lvert \Delta y \rvert \le 15\thinspace\mu\mathrm{m}
$$

$$
(\Delta x)^2 + (\Delta y)^2
< \left(r_i + r_j + 2\thinspace\mu\mathrm{m}\right)^2
$$

**Parallel discovery and serial resolution.** Contact discovery uses two OpenMP
`schedule(static)` passes over sources. The first writes one count per source; a serial prefix sum
and pair allocation then assign disjoint output slices; the second fills each source's slice. No
Python API call or dynamic allocation occurs inside either OpenMP pass, and each worker writes only
its own count or preallocated slice. The GIL is released for this discovery work.

State modification is intentionally not parallel. With the GIL held, the serial resolver traverses
the discovered ordered pairs, rechecks mutable liveness, and requests a Python random draw only for
an HGT-eligible pair. Immediate ordered mutations preserve same-tick HGT donor cascades,
capacity-bounded sharing, T6SS lethal transitions, and predation consumption ordering.

**Events and serialized statistics.** The six returned non-negative integer counters have these
meanings:

| Counter | Meaning |
|---|---|
| `contact_pairs` | Ordered pairs still live when serial resolution reaches them. |
| `t6ss_hits` | Strictly positive T6SS doses applied to live targets. |
| `t6ss_deaths` | T6SS transitions from positive to nonpositive target food. |
| `hgt_transfers` | False-to-true phage-resistance acquisitions. |
| `share_events` | Successful capacity-bounded same-species transfers. |
| `predation_events` | Successful size- and hunger-gated prey consumptions. |

The Global Stats CSV section (Section 7) serializes per-update `t6ss_deaths` and `hgt_transfers`
plus cumulative `overall_t6ss_deaths` and `overall_hgt_transfers`. Legacy imports that lack these
fields initialize them to zero.

**QA, scale fixture, and benchmark.** `bacteria/test_physics_kernel.py` covers a medium 50-tick
per-tick Python/native parity sequence with grow/shrink/regrow transitions; direct malformed-ABI
calls; signed and overflow cell coordinates; transactional callback failure; buffer-tail and event
isolation; retained Python allocation; native stress; the required scale fixture; and an opt-in
benchmark. The kernel workflow is configured to run the physics suite with OpenMP on Windows,
Linux, and macOS and with a serial extension on Linux. The serial job also configures the
1,000-agent/500-tick scale test, and a separate workflow job defines serial Linux Valgrind checks
for definite or indirect losses. These statements describe configuration, not hosted results. The
required scale fixture can also be run directly:

```bash
BACSANDBOX_RUN_PHYSICS_SCALE=1 python -m unittest bacteria.test_physics_kernel.TestPhysicsScale.test_thousand_agents_five_hundred_ticks
```

The required scale test uses a fixed seed and exactly five labels: *Bacteroides thetaiotaomicron*
VPI-5482, *Escherichia coli* K-12, *Faecalibacterium prausnitzii* A2-165,
*Bifidobacterium longum* NCC 2705, and *Akkermansia muciniphila* MucT. Each label has 200 agents,
for 1,000 agents over 500 contact ticks. The traits are explicitly synthetic test inputs, not
empirical species-data rows. Movement, metabolism, and environmental stochasticity are disabled,
and deterministic replacement turnover isolates real contact dispatch.

The opt-in benchmark uses a contact-only 1,000-agent mixed dense/sparse fixture for eight ticks,
with warmups, parity checks before timing, and three measured repeats:

```bash
BACSANDBOX_RUN_PHYSICS_BENCHMARK=1 python -m unittest bacteria.test_physics_kernel
```

Its measurements cover contact dispatch only, excluding movement, metabolism, ECM, and
environment work. They are environment- and workload-specific, not performance guarantees, and
the synthetic fixture is not a biological prediction or a full-`Simulation` benchmark.

### 10.3 Reduced Constraint-Based Metabolic Kernel and Cache

**Scope and rollout boundary.** BacSandbox's production simulation continues to use the existing
phenomenological, pg-based metabolism path by default. This component provides a deliberately
reduced, three-substrate constraint-based bookkeeping model through the opt-in
`MetabolicOrchestrator` programmatic API. A caller may now attach an `FbaIntegrationAdapter` through
the keyword-only `Simulation(..., metabolic_adapter=adapter)` argument. On that explicit path,
`Simulation.update` dispatches a live molar transaction and bypasses legacy food-node generation,
legacy pg reserve changes, and `BacteriaSpecies.step`. With no adapter, behavior remains on the
legacy pg/`FoodNode` path. The application does not construct or configure an adapter through its
UI. This reduced integration is not a GEM, SBML reconstruction, genome-scale model, or full dFBA
implementation.<sup>[41](#ref-41)</sup>

The orchestrator constructor is:

```python
MetabolicOrchestrator(
  mode='legacy', backend='auto', model=None, priors=None, phenotype_config=None
)
```

Modes have the following contracts:

| Mode | Behavior |
|---|---|
| `legacy` | Default no-op. It does not import a solver backend and leaves the production pg-based path unchanged. |
| `shadow` | Runs the selected solver, restores the supplied concentrations and biomass after a successful batch, and retains the diagnostic growth, biomass-delta, substrate-consumption, and acetate-output arrays. |
| `fba` | Runs the selected solver and commits all supplied mutable arrays. The name denotes this reduced direct API, not genome-scale dFBA. |

Outside `legacy`, callers must explicitly supply the model, priors, and phenotype configuration.
Backend `python` uses the reference implementation and never imports the extension. Backend `cpp`
requires the native `MetabolicKernel`. Backend `auto` falls back to Python only when
`bacteria_c_kernel` itself is absent or when an installed unified module lacks `MetabolicKernel`;
broken dependencies, validation failures, solver failures, and other native runtime errors
propagate. These modes and backends are programmatic controls, not app settings.

**Live spatial configuration and state ownership.** `FbaSpatialConfig` makes the live boundary
explicit. `environment_volume_l` must be finite and strictly positive, and every occupied sparse
cell receives that declared volume in liters; the two-dimensional geometry does not infer a depth
or volume. The configuration also owns the three ordered `SubstrateLayerSpec` records,
low-substrate threshold, depletion safety, maximum substeps, cache switch, and concentration-space
conservation tolerance. The substrate names must exactly match `generic_calorie`,
`complex_glycan`, and `mdcf_prebiotic` in that order. `conservation_atol_m` defaults to
$10^{-15}\thinspace\mathrm{M}$. The requested ceiling was `1e-10 M`, so the implemented default is
stricter.<sup>[41](#ref-41)</sup>

Continuous-substrate snapshots include registration specifications, canonical layer order, and
all sparse values. Restoring a snapshot taken before first use therefore removes registration as
well as values if that first transaction fails. Single-cell deposits validate every resulting
concentration, and batch replacements validate every proposed value, before mutating any substrate
map. A positive nonzero value below its layer's retention floor raises an error for the whole
operation rather than being silently discarded; an exact zero remains the explicit
absence/removal value.
<sup>[41](#ref-41)</sup>

`FbaAgentStateStore` is a sidecar keyed by live `Bacterium` identity, so reduced state remains
separate from each cell's legacy pg reserve. It reconciles additions and removals against the
current member lists and commits immutable `FbaAgentState` records containing biomass, growth rate,
biomass delta, three reported substrate amounts, and acetate bookkeeping. It also snapshots and
restores those records without moving pooled or live agent objects.<sup>[41](#ref-41)</sup>

**Deterministic spatial batch.** Batch construction first flattens agents in
`Simulation.species` order and then in each species' member order. It maps each agent through the
Section 10.1 floor rule, sorts the unique coordinate pairs lexicographically, and assigns shared
environment IDs from that sorted coordinate list. Repeated construction from unchanged state
therefore preserves both agent order and environment numbering.<sup>[41](#ref-41)</sup> For $N$
agents and $E$ occupied coordinates, the builder allocates C-contiguous NumPy arrays with these
exact contracts:

| Buffer | NumPy dtype | Shape |
|---|---|---|
| `phenotype_id` | `int32` | `(N,)` |
| `environment_id` | `int32` | `(N,)` |
| `concentrations_m` | `float64` | `(E, 3)` |
| `environment_volume_l` | `float64` | `(E,)` |
| `biomass_gdw` | `float64` | `(N,)` |
| `growth_rate_s` | `float64` | `(N,)` |
| `biomass_delta_gdw` | `float64` | `(N,)` |
| `substrate_consumed_mol` | `float64` | `(N, 3)` |
| `acetate_secreted_mol` | `float64` | `(N,)` |

The sampled concentration matrix is fresh, owning, writable storage. The adapter passes these
arrays to `MetabolicOrchestrator.advance_batch` in the direct ABI order documented below and does
not infer or coerce an alternate unit system.

**Diet bolus deposition.** A `DietBolus` carries an event ID, one versioned translator diet name,
total dry feed mass in grams, and a nonempty footprint of distinct integer ECM coordinates with
finite positive weights. The weights must be normalized to one within an absolute tolerance of
$10^{-12}$. Let $F$ be the number of footprint cells, $V_{\mathrm{cell}}$ the configured volume
of each cell in liters, $C_k^{\mathrm{tr}}$ the translator's additive concentration for substrate
$k$, and $w_j$ the weight of footprint cell $j$. The implementation uses<sup>[41](#ref-41), [47](#ref-47)</sup>

$$
\sum_{j=1}^{F}w_j=1,
\qquad
V_{\mathrm{tr}}=F V_{\mathrm{cell}},
\qquad
n_k^{\mathrm{bolus}}=C_k^{\mathrm{tr}}V_{\mathrm{tr}}.
$$

It then deposits the concentration increment

$$
\Delta C_{j,k}
=\frac{n_k^{\mathrm{bolus}}w_j}{V_{\mathrm{cell}}},
\qquad
\sum_{j=1}^{F}V_{\mathrm{cell}}\Delta C_{j,k}
=n_k^{\mathrm{bolus}}.
$$

Here $C$ is in $\mathrm{mol}/\mathrm{L}$ and $n$ is in $\mathrm{mol}$ of represented
glucose-equivalent substrate. This distributes the translator-owned amount rather than applying
the full bolus independently to every cell. Before writing, the direct deposition function samples
the complete footprint, computes all weighted shares, and stages the complete updated
`(F, 3)` concentration matrix. One ECM batch replacement then validates every staged value,
including every layer's retention floor, before replacing any cell. A sub-floor share in a later
footprint row therefore rejects the direct call without retaining a write to any earlier row. The
bolus queue is consumed only after the enclosing live transaction commits.<sup>[41](#ref-41)</sup>

**Separate conservation checks.** After the solver returns and before sparse-map or sidecar commit,
the adapter requires finite, nonnegative amount arrays and nonincreasing substrate concentrations.
For every environment $e$ and substrate $k$, it independently compares represented depletion with
the agents' reported consumption:<sup>[41](#ref-41)</sup>

$$
n_{e,k}^{\mathrm{depleted}}
=\left(C_{e,k}^{\mathrm{before}}-C_{e,k}^{\mathrm{after}}\right)V_e,
\qquad
n_{e,k}^{\mathrm{reported}}
=\sum_{i\thinspace:\thinspace e_i=e}\Delta n_{i,k}.
$$

With configured concentration tolerance $\varepsilon_C$, acceptance requires

$$
\left\lvert
n_{e,k}^{\mathrm{depleted}}-n_{e,k}^{\mathrm{reported}}
\right\rvert
\le \varepsilon_C V_e.
$$

Thus the represented conservation allowance is an amount in moles obtained from the configured
concentration tolerance and each environment's explicit volume; it is not a fixed mole allowance
independent of volume.<sup>[41](#ref-41)</sup>

Biomass has its own per-agent identity and a fixed absolute tolerance of
$10^{-18}\thinspace\mathrm{gDW}$:

$$
B_i^{\mathrm{after}}-B_i^{\mathrm{before}}=\Delta B_i.
$$

Molar concentration, substrate moles, biomass in $\mathrm{gDW}$, and acetate moles are not added
into one cross-unit balance. `acetate_secreted_mol` is a separate reduced-network bookkeeping
output; the live path records it in FBA state/statistics and does not deposit it into the legacy pg
`ecm.acetate` signal map.<sup>[41](#ref-41)</sup>

**Growth-only lifecycle boundary.** `GrowthOnlyFbaLifecyclePolicy` requires the caller to supply a
finite positive `initial_biomass_gdw`; the implementation supplies no universal cell-biomass
calibration. Its division threshold is exactly twice that initial value:<sup>[41](#ref-41)</sup>

$$
B_{\mathrm{division}}=2B_{\mathrm{initial}}.
$$

At or above that threshold, and only while the carrying-capacity factor is positive, the policy
returns division probability one. Fission splits biomass and every cumulative reduced-state amount
equally between parent and child. Render, motility, and chemotactic scales remain one, and
`starved` is always false. Revision 1 has no maintenance, damage, starvation, or death reaction, so
this policy does not invent any of those processes. Hard carrying-capacity culling remains a
separate `Simulation` rule. Consequently, live-FBA extinction dynamics are not validated and must
not be interpreted using the legacy pg sweep's extinction results.<sup>[41](#ref-41)</sup>

**Complete live-tick transaction.** `FbaIntegrationAdapter.advance_tick` first snapshots substrate
registration and every continuous substrate map together with the sidecar state and an opaque
solver-cache snapshot; a failure in registration, reconciliation, translation, bolus deposition,
Python substrate diffusion, solve, validation, grid commit, or state commit restores them and
re-raises the original exception.
`Simulation._update_fba` wraps that narrower adapter transaction in a full pre-tick snapshot. It
restores all of the following on any exception or other `BaseException` from the live update:
<sup>[41](#ref-41)</sup>

* Python RNG state and the pending bolus queue.
* Registered substrate specifications, all three complete substrate maps, all four signal maps, and the last ECM backend marker.
* Sidecar states and the exact pre-tick solver-cache state.
* Species/member collections and attributes, and every pre-tick agent's attributes.
* Global statistics, per-species step statistics, environment values, `last_fba_result`, and `step_idx`.

For a `Simulation`-owned tick, adapter notifications are deferred. Exactly one success event is
emitted only after adapter commit, lifecycle changes, statistics, environment updates, and the step
index have all committed. Any outer failure restores the snapshot, emits exactly one rollback event
with phase `simulation`, and uses a bare re-raise to preserve the original exception object, type,
and message. Adapter rollback and the outer `Simulation` rollback each restore the exact cache state
captured at their transaction boundary, so a valid warmed cache is preserved rather than cleared.
If no kernel existed at the first snapshot, the opaque snapshot remains lazy as `None`; if that
failed transaction creates a kernel, restoring `None` returns the created cache and all counters to
their empty initial state without having initialized a backend merely to take the snapshot.
<sup>[41](#ref-41)</sup>

The Pygame loop is the presentation boundary. After rollback it catches ordinary exceptions only
when an adapter is attached, pauses the simulation, records a displayable error string, and logs
the traceback. That error remains visible while paused; only a later successful running update or
simulation reinitialization clears it. Legacy-mode exceptions still propagate from the application
loop.<sup>[41](#ref-41)</sup>

**Reset and serialization.** `Simulation.reset` clears pending boluses, `last_fba_result`, all
continuous-substrate values, all sidecar states, and any already-created solver cache deliberately
for the new run, without forcing cache initialization. It retains the registered substrate
definitions for the attached adapter. Reset also recreates empty statistics, including the five
FBA series documented in Section 7. CSV export/import includes `fba_biomass_gdw`,
`fba_acetate_mol`, and the three per-tick `fba_consumed_*_mol` fields.<sup>[41](#ref-41)</sup>

**Live-integration validation status.** `bacteria.test_metabolic_integration` covers floor mapping,
strict and atomic substrate mutation, registration-restoring first-use rollback, Python transport,
weighted deposition, a later-share retention-floor partial-write regression, deterministic
shared-cell batches, growth-only fission, commit, conservation rejection, adapter rollback,
reset/cache clearing, FBA CSV fields, and legacy-pg bypass. `bacteria.test_integration_stress`
injects failures from translation through late statistics commit, verifies deferred event
cardinality, exact first-use and warmed-cache rollback restoration, full owner-state restoration,
and original-error propagation, exercises exactly 5,001 agents, checks bounded Python retention,
and runs warmed-cache snapshot/restore round trips on Python and conditionally on the native kernel.
It also conditionally compares Python/native live outputs at `rtol=0`, `atol=1e-12` when the native
metabolic kernel is built. These are executable software contracts, not validation of biological
extinction behavior.

The OpenMP and serial jobs in `.github/workflows/kernel.yml` are configured to run these exact
commands in addition to the kernel suites:<sup>[41](#ref-41)</sup>

```bash
python -m unittest bacteria.test_diet_translator
python -m unittest bacteria.test_metabolic_integration
python -m unittest bacteria.test_integration_stress
python -m unittest bacteria.test_comparative_sweeps
```

Both build paths first assert that the unified extension exposes `MetabolicKernel`; the serial job
also configures the required metabolic and physics scale fixtures. A separate serial Valgrind job
is defined for focused native allocation stress. Workflow configuration does not establish that
hosted CI or Valgrind passed for a particular commit, and no such result is asserted here.

**Versioned assets and reduced topology.** The reduced-model interface and topology, flat/CSR ABI,
deterministic posterior-mean constraint adapter, tolerance-neighbor cache with bounded LRU, and
adaptive substep design are original, nonempirical BacSandbox designs.<sup>[41](#ref-41)</sup>

| Asset | Revision-1 contents | Interpretation |
|---|---|---|
| `bacteria/models/reduced_metabolism_v1.json` | One model, 5 metabolites, and 8 reactions with finite lower and upper bounds. | Three exchanges feed three corresponding routes; each route balances one exchange unit to biomass precursor plus `0.25` acetate bookkeeping. Biomass is the sole objective. |
| `bacteria/models/metabolic_priors_v1.json` | One neutral row for each of the three substrates; every `alpha` and `beta` is positive and every `max_inhibition` is zero. | The rows make no empirical interaction claim and produce no inhibition in the shipped configuration. |

The reduced model declares schema version 1 and model revision 1.
The fixed substrate order is `generic_calorie`, `complex_glycan`, then `mdcf_prebiotic`.
Concentrations are molar glucose-equivalent monomer concentrations. In particular, glycan and MDCF
values are never literal glycan/polymer molarity. Every batch requires an explicit environment
volume in $\mathrm{L}$; there is no implicit voxel volume.

Exchange, route, and acetate fluxes use specific units of
$\mathrm{mol}/\mathrm{gDW}/\mathrm{s}$ of glucose-equivalent monomer. The biomass pseudo-flux uses
$\mathrm{s}^{-1}$. At full finite uptake in the shipped model, each exchange and route flux is
`1`, biomass flux is `3`, and acetate flux is `0.75`. These are analytical properties of the
nonempirical bookkeeping network, not biological predictions.

`bacteria/metabolism.py` validates both versioned assets and converts sparse stoichiometry into
exact-dtype, C-contiguous, read-only NumPy arrays. Multiple reduced models may be concatenated in
the same storage. Let `M` be model count, `R` total reaction count, `K` total metabolite count, and
`Z` stored stoichiometric nonzero count. The native constructor is:

```text
MetabolicKernel(
  model_reaction_offsets, model_metabolite_offsets, stoich_row_offsets,
  stoich_column_indices, stoich_values, lower_bounds, upper_bounds,
  objective_coefficients, exchange_reaction_indices,
  biomass_reaction_indices, acetate_reaction_indices,
  cache_tolerance_molar, cache_capacity)
```

Constructor buffers are copied into RAII-owned native storage:

| Buffer | Exact dtype | Exact shape | Meaning |
|---|---|---|---|
| `model_reaction_offsets` | `int64` | `(M + 1,)` | Concatenated reaction boundaries by model. |
| `model_metabolite_offsets` | `int64` | `(M + 1,)` | Concatenated metabolite boundaries by model. |
| `stoich_row_offsets` | `int64` | `(K + 1,)` | CSR row boundaries. |
| `stoich_column_indices` | `int32` | `(Z,)` | Global reaction columns for CSR entries. |
| `stoich_values` | `float64` | `(Z,)` | CSR coefficients. |
| `lower_bounds` | `float64` | `(R,)` | Finite reaction lower bounds. |
| `upper_bounds` | `float64` | `(R,)` | Finite reaction upper bounds. |
| `objective_coefficients` | `float64` | `(R,)` | Concatenated objective vector. |
| `exchange_reaction_indices` | `int32` | `(M, 3)` | Three global exchange indices per model in fixed substrate order. |
| `biomass_reaction_indices` | `int32` | `(M,)` | Global biomass pseudo-reaction index per model. |
| `acetate_reaction_indices` | `int32` | `(M,)` | Global acetate bookkeeping index per model. |

`cache_tolerance_molar` is finite and positive, and `cache_capacity` is positive and bounded by the
Python asset/configuration validation. Native construction precomputes one dense LP matrix from
each model's CSR rows. Phenotypes are then configured exactly as follows:

```text
configure_phenotypes(
  model_indices, uptake_vmax, uptake_km, inhibition_alpha,
  inhibition_beta, max_inhibition, model_revision, prior_revision)
```

| Argument | Exact dtype | Exact shape or scalar contract |
|---|---|---|
| `model_indices` | `int32` | `(P,)`; one valid model index per phenotype. |
| `uptake_vmax` | `float64` | `(P, 3)`; finite and nonnegative. |
| `uptake_km` | `float64` | `(P, 3)`; finite and positive. |
| `inhibition_alpha` | `float64` | `(P, 3)`; finite and positive. |
| `inhibition_beta` | `float64` | `(P, 3)`; finite and positive. |
| `max_inhibition` | `float64` | `(P, 3)`; finite values from zero through one. |
| `model_revision` | `uint64` | Nonnegative revision scalar. |
| `prior_revision` | `uint64` | Nonnegative revision scalar. |

Reconfiguration copies the phenotype data and resets every cache entry, counter, and LRU sequence.
`clear_cache()` performs the same cache reset. `cache_stats()` returns
`(entries, hits, misses, evictions, lp_solves)`.

Both the Python and C++ kernels expose `snapshot_cache()` and `restore_cache(snapshot)`. Their
backend-specific snapshot objects are opaque to callers and preserve all entries and full flux
payloads, hit/miss/eviction/LP-solve counters, each entry's insertion and last-used sequences, and
the global LRU sequence. Restoration therefore recovers eviction order and future hit behavior as
well as the visible statistics. `MetabolicOrchestrator.snapshot_cache()` does not force lazy kernel
construction and returns `None` before first use; `restore_cache(None)` clears a kernel created
during the failed transaction, while a non-`None` snapshot restores the warmed backend exactly.
<sup>[41](#ref-41)</sup>

The direct solve binding is:

```text
solve_flux(phenotype_index, concentrations_m, cache_enabled=True)
  -> (flux_values, cache_hit)
```

`concentrations_m` must be exact `float64`, C-contiguous, shape `(3,)`, finite, and nonnegative.
The returned `flux_values` is an independent `float64` vector of shape `(R,)` across all
concatenated models; reactions outside the phenotype's selected model are zero. `cache_hit` is a
Boolean. The orchestrator wrapper defaults `phenotype_index` to zero and requires
`concentrations_m` outside `legacy` mode.

The exact batch binding order is:

```text
advance_batch(
    phenotype_id, environment_id, concentrations_m, environment_volume_l,
    biomass_gdw, growth_rate_s, biomass_delta_gdw,
    substrate_consumed_mol, acetate_secreted_mol, dt_seconds, *,
    low_substrate_threshold_m=1e-4, depletion_safety,
    max_substeps, cache_enabled)
  -> (substeps, lp_solves, cache_hits, cache_misses, cache_evictions,
      minimum_substep_seconds)
```

The `MetabolicOrchestrator.advance_batch` wrapper uses the same buffer order and defaults
`low_substrate_threshold_m=1e-4`, `depletion_safety=0.9`, `max_substeps=256`, and
`cache_enabled=True`.

| Buffer | Exact dtype | Exact shape | Ownership during a successful call |
|---|---|---|---|
| `phenotype_id` | `int32` | `(N,)` | Read-only configured phenotype IDs. |
| `environment_id` | `int32` | `(N,)` | Read-only environment IDs; multiple agents may share an environment. |
| `concentrations_m` | `float64` | `(E, 3)` | Writable; committed concentrations by environment and fixed substrate order. |
| `environment_volume_l` | `float64` | `(E,)` | Read-only, finite, positive explicit volumes in $\mathrm{L}$. |
| `biomass_gdw` | `float64` | `(N,)` | Writable, nonnegative biomass in $\mathrm{gDW}$. |
| `growth_rate_s` | `float64` | `(N,)` | Writable output reset to the time-weighted mean growth rate in $\mathrm{s}^{-1}$. |
| `biomass_delta_gdw` | `float64` | `(N,)` | Writable output reset and accumulated in $\mathrm{gDW}$. |
| `substrate_consumed_mol` | `float64` | `(N, 3)` | Writable output reset and accumulated in moles of glucose-equivalent monomer. |
| `acetate_secreted_mol` | `float64` | `(N,)` | Writable output reset and accumulated in moles of acetate bookkeeping. |

`dt_seconds` is real elapsed time in seconds and must be finite and nonnegative;
`low_substrate_threshold_m` is finite and nonnegative; `depletion_safety` is finite in `(0, 1]`;
`max_substeps` is a positive integer; and `cache_enabled` is Boolean. Across native construction,
configuration, solve, and batch calls, arrays require exact dtype, dimension, shape, and
C-contiguity under pybind11 `noconvert`; mutable arrays must be writable. Validation rejects
nonfinite or out-of-range values, invalid IDs, invalid scalars, and any unsafe memory overlap
involving a mutable array.

**Bounded LP solver.** For the selected reduced model, the deterministic embedded two-phase dense
simplex maximizes the supplied objective under steady-state balance and finite reaction bounds:

$$
\max_{v}\quad c^{\mathsf{T}}v
\qquad \text{subject to}\qquad
Sv=0,\qquad \ell\le v\le u
$$

The implementation shifts lower bounds with

$$
x=v-\ell\ge 0,\qquad
Sx=-S\ell,\qquad
x\le u-\ell
$$

and represents each equality as paired inequalities plus one finite upper constraint per shifted
reaction. The solver is dependency-free C++17 and uses the lowest stable index for deterministic
ties. Its tolerance is $10^{-10}$ and is strictly a numerical simplex tolerance, not a biological
parameter. It reports explicit `optimal`, `infeasible`, `unbounded`, and `numerical_error` states;
the production asset validator rejects nonfinite bounds, so shipped reduced models are bounded.
This small dense solver is tested for the reduced assets and does not claim genome-scale
robustness.

**Deterministic posterior-mean exchange bounds.** There is no random sampling. For phenotype
$p$ and substrate $k$, the inhibition-probability mean and retained uptake multiplier are

$$
\bar{p}_{pk}=\frac{\alpha_{pk}}{\alpha_{pk}+\beta_{pk}},\qquad
m_{pk}=1-\bar{p}_{pk}e_{pk}
$$

where $e_{pk}$ is `max_inhibition`. At concentration $C_k$, the effective exchange upper bound is

$$
u_{pk}^{\mathrm{eff}}(C_k)=\min\left(
u_{k}^{\mathrm{model}},
V_{pk}^{\max}\frac{C_k}{K_{pk}+C_k}m_{pk}
\right)
$$

Zero concentration therefore gives a zero exchange bound. In the shipped neutral priors,
$e_{pk}=0$ and hence $m_{pk}=1$ for every phenotype/substrate row.

**Tolerant cache and batch ownership.** Each `MetabolicKernel` owns a private
`std::unordered_map`; cache data never crosses kernel instances or runs. With
$\delta_C=10^{-5}\thinspace\mathrm{M}$, each concentration component is quantized as

$$
q_k=\left\lfloor\frac{C_k}{10^{-5}\thinspace\mathrm{M}}\right\rfloor
$$

The key is `(model_revision, prior_revision, phenotype_index, q_0, q_1, q_2)`. Lookup probes all
27 adjacent three-dimensional bin combinations and accepts a stored concentration only when

$$
\max_{k\in\lbrace 0,1,2\rbrace}
\lvert C'_k-C_k\rvert\le 10^{-5}\thinspace\mathrm{M}
$$

The stored full flux vector is intentionally reused even though continuously derived effective
bounds can differ within that tolerance. If several entries match, lookup chooses the oldest
insertion sequence and then updates its last-used sequence. Capacity is bounded; eviction chooses
the least recently used entry, breaking a last-use tie by oldest insertion sequence. Each entry
stores exact concentrations, effective exchange bounds, the full concatenated flux vector, and
both sequences. Model/prior revisions are part of the key, while reconfiguration and
`clear_cache()` reset entries and counters.

Within a batch substep, misses with exactly identical phenotype, concentrations, and effective
bounds are deduplicated. `cache_misses` still counts missed agents, so `lp_solves` can be smaller.
Cache lookup and insertion remain serial. Unique misses are solved with OpenMP
`schedule(static)` using preallocated per-request workspaces; no Python API call, shared-cache
mutation, or dynamic allocation occurs inside that parallel region. The GIL is released around
pure native solve and batch work.

**Adaptive metabolic integration.** Adaptive substeps apply only inside the requested metabolic
batch interval. They do not change the global movement, contact, or ECM elapsed-time step. For
environment $e$ and substrate $k$, fixed-order reduction computes the depletion rate

$$
D_{ek}=\frac{\sum_{i\mid e_i=e}u_{ik}B_i}{V_e}
$$

where $u_{ik}$ is specific uptake in
$\mathrm{mol}/\mathrm{gDW}/\mathrm{s}$, $B_i$ is biomass in $\mathrm{gDW}$, $V_e$ is the
explicit environment volume in $\mathrm{L}$, and $D_{ek}$ is in $\mathrm{M}/\mathrm{s}$.

Let $\sigma$ be `depletion_safety`. The per-environment safety multiplier is

$$
\lambda_{ek}=\sigma\quad\mathrm{when}\quad
k\in\lbrace 1,2\rbrace\text{ and }
C_{ek}\lt 10^{-4}\thinspace\mathrm{M}=0.1\thinspace\mathrm{mM},
\qquad \lambda_{ek}=1\quad\text{otherwise}
$$

Thus glycan and MDCF, indices 1 and 2, receive the additional factor only when strictly below the
threshold. For every positive depletion rate, the nonnegative concentration limit always applies:

$$
\Delta t_{\mathrm{sub}}=\min\left(
\Delta t_{\mathrm{remaining}},
\min_{e,k\thinspace:\thinspace D_{ek}\gt 0}
\lambda_{ek}\frac{C_{ek}}{D_{ek}}
\right)
$$

An empty set of positive rates contributes no tighter limit. Zero concentration produces zero
uptake. There is no lower time-step clamp that could override depletion safety; failure to make
numerical progress or exhaustion of `max_substeps` raises an error.

For substep $r$, explicit Euler updates per-agent substrate use, biomass, and acetate as

$$
\Delta n_{ik}^{(r)}=u_{ik}^{(r)}B_i^{(r)}\Delta t_r,\qquad
\Delta B_i^{(r)}=g_i^{(r)}B_i^{(r)}\Delta t_r,\qquad
\Delta a_i^{(r)}=a_i^{\mathrm{flux},(r)}B_i^{(r)}\Delta t_r
$$

and shared environments and biomass advance by

$$
C_{ek}^{(r+1)}=\max\left(
0,
C_{ek}^{(r)}-\frac{\sum_{i\mid e_i=e}\Delta n_{ik}^{(r)}}{V_e}
\right),\qquad
B_i^{(r+1)}=B_i^{(r)}+\Delta B_i^{(r)}
$$

If floating-point summation requests slightly more substrate than an environment contains, all
requests for that environment/substrate are scaled proportionally to availability as a
deterministic roundoff backstop. `substrate_consumed_mol`, `biomass_delta_gdw`, and
`acetate_secreted_mol` accumulate across substeps. The final growth output is the time-weighted
mean

$$
\bar{g}_i=\frac{1}{\Delta t}\sum_r g_i^{(r)}\Delta t_r
$$

Native calls snapshot all six mutable arrays together with the complete cache, counters, and LRU
sequence. Any validation, solver, substep, or allocation failure restores that state
transactionally before propagating. Successful `shadow` orchestration performs its separate
documented restoration of concentrations and biomass while preserving diagnostics.

**QA, scale fixture, and reduced-model cache benchmark.**
`bacteria/test_metabolic_architecture.py` and `bacteria/test_metabolic_kernel.py` cover versioned
asset validation, seeded Python/native references, $Sv=0$ residuals, unsaturated cache-tolerance
boundaries, deterministic LRU/revision/deduplication behavior, adaptive shared-environment
depletion, exact ABI and aliases, transactional arrays plus cache state, ownership, native stress,
scale parity, and opt-in performance measurement. Batch parity compares mutable arrays with the
subsystem's declared tolerances and checks represented-substrate depletion against explicit
environment volumes; these are model-level bookkeeping contracts, not whole-simulation mass
conservation.

The direct-solve parity fixture combines 20 fixed cache/adaptive-boundary states with 100
seed-3103 log-sweep states. Python and native fluxes are checked against the reduced analytical
reference with `rtol=0`, `atol=1e-9`, and against each other with `rtol=0`, `atol=1e-10`.

The required synthetic scale fixture runs 1,000 agents for 12 repeated adaptive ticks and compares
the Python and native summaries and mutable arrays after every tick. It is not a full-`Simulation`
benchmark, genome-scale workload, or biological prediction. Run it with:

```bash
BACSANDBOX_RUN_METABOLIC_SCALE=1 python -m unittest bacteria.test_metabolic_kernel.TestMetabolicScale.test_thousand_agents_repeated_adaptive_ticks
```

The opt-in reduced-model cache benchmark uses 64 tolerance-separated concentration states,
repeats each state 32 times, warms both uncached and cached paths, verifies output parity, and
compares medians across five measured samples. The environment-variable name retains the existing
reduced-API `FBA` label; it does not denote genome-scale dFBA:

```bash
BACSANDBOX_RUN_FBA_BENCHMARK=1 python -m unittest bacteria.test_metabolic_kernel.TestMetabolicCacheBenchmark
```

Performance measurements are environment- and workload-specific, not guarantees. The kernel
workflow is configured to build required OpenMP variants on Windows, Linux, and macOS plus a
forced-serial Linux variant, run the metabolic suite and required scale fixture, and define
2,000-iteration metabolic, physics, and ECM stress workloads under serial Linux Valgrind checks for
definite and indirect losses. No hosted or Valgrind execution result is asserted here.

### 10.4 Empirically Sourced Dietary Parameterization and Own-Model Molar Translation

#### 10.4.1 Versioned asset, engine, and injection contract

Revision 1 of the dietary translation component ships the versioned asset
`bacteria/models/diet_profiles.json`, with schema version 1, profile ID `diet_profiles_v1`, profile
revision 1, and provenance `Own (BacSandbox model)`.
`bacteria.diet_translator.DietTranslationEngine` eagerly loads and validates this asset. These
schema and API choices are original BacSandbox design.<sup>[47](#ref-47)</sup>

The public bolus API is

```text
get_molar_substrate_bounds(diet_name, feed_mass_g, volume_L, dt_s)
  -> float64[3]
```

Every successful call returns a fresh, writable, owning, C-contiguous NumPy `float64` array with
shape `(3,)`. Its exact ABI order is `generic_calorie`, `complex_glycan`, then
`mdcf_prebiotic`.<sup>[47](#ref-47)</sup> The vector is one environment's additive concentration
increment in $\mathrm{M}$; it is not an exchange flux and is not itself a valid `concentrations_m`
argument to `advance_batch`. A caller allocates writable, C-contiguous `(E, 3)` concentration
storage and adds the returned vector into the selected environment row before advancing the batch:

```python
increment_m = engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, dt_s)
concentrations_m = np.zeros((environment_count, 3), dtype=np.float64, order='C')
concentrations_m[environment_index] += increment_m
```

`volume_L` is explicit because the two-dimensional simulation has no physical depth-to-liter
model. `feed_mass_g` is the total dry bolus delivered by that call, not a rate. `dt_s` must be a
finite positive real number, but it is validation-only and is absent from the bolus equation. No
continuous-feed API currently exists; such an API would need an explicitly named dry-mass rate and
would multiply that rate by elapsed seconds before translation.<sup>[47](#ref-47)</sup>

The schema is intentionally strict. It rejects duplicate JSON keys at every depth; missing or
unknown fields; wrong schema, profile, revision, provenance, diet order, component order,
ingredient order, or substrate order; unresolved source hooks and cross-references; and nonfinite
or out-of-range numbers. Every component and ingredient slot must be present, including explicit
zeros. There are no profile aliases, clipping rules, fallback profiles, inferred components, or
automatic normalization. The sole normalization divides source-reported recipe parts by that
recipe's explicit reported total, absorbing the small rounding offsets in the published tables
while preserving their stated parts.<sup>[42](#ref-42), [47](#ref-47)</sup>

#### 10.4.2 Ingredient-to-substrate equations and units

Let $r_{d,i}$ be the source-reported dry recipe part for diet $d$ and ingredient $i$, $R_d$ the
corresponding reported total, $f_{i,j}$ the mass fraction of measured or explicitly modeled
component $j$ in ingredient $i$, and $A_{j,k}$ the one-hot indicator routing component $j$ to ABI
substrate $k$. The normalized dry ingredient fraction is

$$
q_{d,i} = \frac{r_{d,i}}{R_d}
$$

<sup>[42](#ref-42), [47](#ref-47)</sup>

and the represented mass fraction routed to substrate $k$ is

$$
\omega_{d,k} = \sum_i q_{d,i}\sum_j f_{i,j}A_{j,k}.
$$

<sup>[43](#ref-43), [47](#ref-47)</sup>

For dry feed mass $m_{\mathrm{feed}}$ and explicit volume $V$, the additive concentration is

$$
C_{d,k} = \frac{m_{\mathrm{feed}}\thinspace\omega_{d,k}}
{M_{\mathrm{GE}}\thinspace V}.
$$

<sup>[46](#ref-46), [47](#ref-47)</sup>

Here $m_{\mathrm{feed}}$ is in $\mathrm{g}$, $V$ is in $\mathrm{L}$, and the asset fixes the
D-glucose-equivalent molar mass at
$M_{\mathrm{GE}} = 180.156\thinspace\mathrm{g}/\mathrm{mol}$ by rounding NIST's
$180.1559\thinspace\mathrm{g}/\mathrm{mol}$ molecular weight to three decimal places.<sup>[46](#ref-46)</sup> The
output is therefore $\mathrm{mol}/\mathrm{L} = \mathrm{M}$ glucose-equivalent monomer. Glycan
values are monomer-equivalent bookkeeping concentrations, not literal polymer molarities. The
complete recipe-component-routing translator, including use of one molar mass for all represented
classes, is **Own (BacSandbox model)**.<sup>[47](#ref-47)</sup>

The source totals $R_d = 99.94\thinspace\mathrm{g}/100\thinspace\mathrm{g}$ for MDCF-2 and
$R_d = 100.04\thinspace\mathrm{g}/100\thinspace\mathrm{g}$ for RUSF are retained rather than
silently replaced by 100.<sup>[42](#ref-42)</sup> Dividing every reported part by its own total
gives

$$
\sum_i q_{d,i} = 1,
\qquad
\sum_k A_{j,k} = 1.
$$

<sup>[42](#ref-42), [47](#ref-47)</sup>

At the model-definition level, the one-hot map routes every represented component gram once, so
$\sum_k\omega_{d,k} = \sum_i q_{d,i}\sum_j f_{i,j}$. In floating-point implementation, frozen
component/substrate fixtures and linear mass/volume scaling are verified with `rtol=0` and
`atol=1e-15` $\mathrm{M}$. This is model-level, tolerance-aware represented-component
conservation.<sup>[47](#ref-47)</sup> Oil, micronutrients, and unrepresented ingredient mass remain
outside the three carbohydrate pools, so the represented substrate fractions need not sum to one.
The contract is neither whole-bolus nor whole-simulation mass conservation.

The revision-1 component routing is fixed as follows:

| Component $j$ | ABI substrate $k$ | Provenance |
|---|---|---|
| `simple_hexose` | `generic_calorie` | Glucose-equivalent bookkeeping, Own.<sup>[47](#ref-47)</sup> |
| `digestible_starch` | `generic_calorie` | FITDOG total-starch allocation, Own.<sup>[47](#ref-47)</sup> |
| `resistant_starch` | `complex_glycan` | Explicit reserved class; revision-1 values are zero, Own.<sup>[47](#ref-47)</sup> |
| `cellulose` | `complex_glycan` | Measured class from Supplementary Table 9d; routing Own.<sup>[43](#ref-43), [47](#ref-47)</sup> |
| `arabinan` | `complex_glycan` | Measured class from Supplementary Table 9d; routing Own.<sup>[43](#ref-43), [47](#ref-47)</sup> |
| `xylan` | `complex_glycan` | Measured class from Supplementary Table 9d; routing Own.<sup>[43](#ref-43), [47](#ref-47)</sup> |
| `beta_1_4_galactan` | `mdcf_prebiotic` | Measured class from Supplementary Table 9d; routing Own.<sup>[43](#ref-43), [47](#ref-47)</sup> |
| `beta_1_4_mannan` | `mdcf_prebiotic` | Measured class from Supplementary Table 9d; routing Own.<sup>[43](#ref-43), [47](#ref-47)</sup> |

#### 10.4.3 Empirical recipe and ingredient profiles

The MDCF-2 and RUSF parts come from Chen et al. Supplementary Table S1.<sup>[42](#ref-42)</sup>
Zeros below are explicit asset values, not omitted data:

| Dry ingredient ($\mathrm{g}/100\thinspace\mathrm{g}$) | MDCF-2 | RUSF | Source |
|---|---:|---:|---|
| Micronutrient mix | $3.14$ | $3.14$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Chickpea | $10.0$ | $0.0$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Peanut | $10.0$ | $0.0$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Soybean | $8.0$ | $0.0$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Green banana | $19.0$ | $0.0$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Rice | $0.0$ | $18.9$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Lentil | $0.0$ | $21.5$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Powdered skimmed milk | $0.0$ | $10.5$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Added sugar | $29.8$ | $17.0$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| Soybean oil | $20.0$ | $29.0$ | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |
| **Reported total** | **$99.94$** | **$100.04$** | Chen et al. Table S1.<sup>[42](#ref-42)</sup> |

The `malnutrition` profile is not a clinical formulation. It is an **Own (BacSandbox model)**
low-diversity baseline containing
$3.14\thinspace\mathrm{g}$ micronutrient mix,
$18.9\thinspace\mathrm{g}$ rice,
$49.0\thinspace\mathrm{g}$ added sugar, and
$29.0\thinspace\mathrm{g}$ soybean oil per source-style recipe, with reported total
$100.04\thinspace\mathrm{g}$; every other ingredient is explicitly zero.<sup>[47](#ref-47)</sup>

For chickpea, peanut, soybean, green banana, rice, lentil, and milk powder, the asset's ingredient
component fractions are replicate means derived from Hibberd et al. Supplementary Table 9d.
That table reports micrograms of polysaccharide class per milligram of dry ingredient; the numeric
conversion used by the asset is

$$
f_{i,j}\quad[\mathrm{g}/\mathrm{g}]
= 10^{-3}\thinspace\bar{x}_{i,j}\quad[\mu\mathrm{g}/\mathrm{mg}].
$$

<sup>[43](#ref-43)</sup>

The measured classes represented in revision 1 are FITDOG total starch, cellulose, arabinan,
xylan/hemicellulose backbone, beta-1,4-linked galactan, and beta-1,4-linked mannan.<sup>[43](#ref-43)</sup>
The JSON asset preserves the full machine-readable replicate-mean coefficients; prose rounds no
coefficient into a second source of truth.

The source quantified total starch but did not establish a validated digestible-versus-resistant
partition. Revision 1 therefore assigns the measured FITDOG total-starch value to
`digestible_starch` and sets `resistant_starch` to exactly zero. This is an Own allocation and is
not evidence that resistant starch is absent from an ingredient.<sup>[43](#ref-43), [47](#ref-47)</sup>
Added sugar contributes one gram of simple-hexose glucose equivalent per gram of sugar, whereas
soybean oil and micronutrient mix contribute zero to all three carbohydrate pools. Those are also
Own bookkeeping choices, not measurements reported by the recipe source.<sup>[47](#ref-47)</sup>

#### 10.4.4 Frozen MDCF-2 translation fixture

For a $50\thinspace\mathrm{g}$ dry MDCF-2 bolus in $1\thinspace\mathrm{L}$, the verified
revision-1 fixture is:<sup>[47](#ref-47)</sup>

| Output | Concentration ($\mathrm{M}$) | Concentration ($\mathrm{mM}$) |
|---|---:|---:|
| Beta-1,4 galactan component | $0.0002780913974742124$ | $0.2780913974742124$ |
| Beta-1,4 mannan component | $0.000020245970462459088$ | $0.020245970462459088$ |
| Combined `mdcf_prebiotic` | $0.00029833736793667145$ | $0.29833736793667145$ |
| `generic_calorie` | $0.11754291369757508$ | $117.54291369757508$ |
| `complex_glycan` | $0.0012267241245826435$ | $1.2267241245826435$ |

The acceptance test explicitly checks combined `mdcf_prebiotic` error
$\le 10^{-6}\thinspace\mathrm{M}$ and also freezes the component and substrate values with
$10^{-15}\thinspace\mathrm{M}$ absolute checks. These are regression tolerances, not estimates of
experimental uncertainty.<sup>[47](#ref-47)</sup>

#### 10.4.5 Species capacities and evidence boundary

`build_species_phenotype_config` reads the nonnegative, dimensionless CSV values `aff_calorie`,
`aff_glycan`, and `aff_mdcf` in ABI order. For phenotype $p$ and substrate $k$, it constructs

$$
V_{p,k}^{\max} = a_{p,k}V_k^{\mathrm{base}}.
$$

<sup>[47](#ref-47)</sup>

The affinity $a_{p,k}$ scales uptake capacity only. The half-saturation concentration
$K_{p,k}$ is independent of that multiplier.<sup>[47](#ref-47)</sup> The shipped calibration is:

| ABI substrate $k$ | $V_k^{\mathrm{base}}\thinspace(\mathrm{mol}\thinspace\mathrm{GE}/(\mathrm{gDW}\thinspace\mathrm{s}))$ | $K_{p,k}\thinspace(\mathrm{M})$ | Provenance |
|---|---:|---:|---|
| `generic_calorie` | $2.8 \times 10^{-6}$ | $1.0 \times 10^{-4}$ | Own calibration.<sup>[47](#ref-47)</sup> |
| `complex_glycan` | $1.4 \times 10^{-6}$ | $2.0 \times 10^{-4}$ | Own calibration.<sup>[47](#ref-47)</sup> |
| `mdcf_prebiotic` | $1.8 \times 10^{-6}$ | $1.5 \times 10^{-4}$ | Own calibration.<sup>[47](#ref-47)</sup> |

These values and every affinity tier below are **Own (BacSandbox model)** calibrations, not
measured whole-cell kinetic constants.<sup>[47](#ref-47)</sup> With the notation established in
Section 10.3, the kernel applies the actual upper-bound semantics

$$
u_{p,k}^{\mathrm{eff}}(C_k) = \min\left(
u_k^{\mathrm{model}},
V_{p,k}^{\max}\frac{C_k}{K_{p,k}+C_k}m_{p,k}
\right).
$$

<sup>[41](#ref-41), [47](#ref-47)</sup>

Both $u_k^{\mathrm{model}}$ and $V_{p,k}^{\max}$ have units
$\mathrm{mol}\thinspace\mathrm{GE}/(\mathrm{gDW}\thinspace\mathrm{s})$; $C_k$ and $K_{p,k}$
are in $\mathrm{M}$. The final factor $m_{p,k}$ is the dimensionless retained uptake fraction
after deterministic posterior-mean inhibition,

$$
m_{p,k} = 1 - \frac{\alpha_{p,k}}{\alpha_{p,k}+\beta_{p,k}}e_{p,k},
$$

<sup>[41](#ref-41)</sup>

where $e_{p,k}$ is `max_inhibition`. It is not biomass, elapsed time, diet affinity, or an
additional species-capacity tier. The shipped neutral priors have $e_{p,k}=0$ and therefore
$m_{p,k}=1$.<sup>[41](#ref-41), [47](#ref-47)</sup> An affinity of zero instead makes
$V_{p,k}^{\max}$ and the effective uptake bound exactly zero.<sup>[47](#ref-47)</sup>

The four strain-evidence profiles in `species_data.csv` are opt-in with `init=0`:

| CSV profile | Profile type | $(a_{p,0},a_{p,1},a_{p,2})$ | Evidence scope and calibration |
|---|---|---:|---|
| Bg0018-like | Growth-associated MAG, not a cultured BgF5_2 isolate | $(1,2,3)$ | The MDCF tier 3 is an Own qualitative encoding tied narrowly to Bg0018 strain-specific MAG association and conserved-PUL evidence; it is not a kinetic measurement.<sup>[43](#ref-43), [47](#ref-47)</sup> |
| Bg0019-like | Growth-associated MAG, not a cultured BgF5_2 isolate | $(1,2,3)$ | The MDCF tier 3 is an Own qualitative encoding tied narrowly to Bg0019 strain-specific MAG association and conserved-PUL evidence; it is not a kinetic measurement.<sup>[43](#ref-43), [47](#ref-47)</sup> |
| *Segatella copri* BgF5_2 | Cultured representative | $(1,2,4)$ | Direct BgF5_2 glucomannan machinery evidence from Zhou et al.; tier values Own.<sup>[45](#ref-45), [47](#ref-47)</sup> |
| *Segatella copri* DSM 18205 | Cultured MDCF non-utilizer | $(1,1,0)$ | Direct inability to grow on glucomannan from Zhou et al.; tier values Own.<sup>[45](#ref-45), [47](#ref-47)</sup> |

Source-era papers and the MAG labels may use *Prevotella copri*; this document uses the current
name *Segatella copri* for the cultured strains while retaining source identifiers.<sup>[44](#ref-44), [45](#ref-45)</sup>
Utilization is strain-specific. Hibberd et al. supports the growth-associated Bg0018/Bg0019 MAGs
and conservation of relevant polysaccharide-utilization loci.<sup>[43](#ref-43)</sup> Chang et al.
supports diet-dependent MDCF-2 degradation and weight effects mediated by *P. copri* and other
microbiota members; it does not establish universal utilization across the species.<sup>[44](#ref-44)</sup>
Zhou et al. directly supports BgF5_2 PUL7/PUL8 beta-mannan and glucomannan processing machinery and
DSM 18205's inability to grow on glucomannan.<sup>[45](#ref-45)</sup> PUL7/PUL8 are not described
here as generic dual galactan/mannan loci: other conserved PULs provide the evidence boundary for
pectic galactan and arabinogalactan processing.<sup>[43](#ref-43), [45](#ref-45)</sup>

The separate *Escherichia coli* K-12 row is also configured with `init=0` and affinity vector
$(1,0,0)$. Its `aff_mdcf=0` makes MDCF-prebiotic uptake exactly zero in both implemented diet
paths for this configuration. This is an **Own (BacSandbox model)** zero-founder/configuration
choice, not evidence of a measured K-12 kinetic constant or a species-wide utilization
claim.<sup>[47](#ref-47)</sup>

#### 10.4.6 QA, memory layout, and native floating-point contract

`bacteria/test_diet_translator.py` covers malformed assets, exact translation formulas, array ABI
and ownership, package-resource inclusion, Python/native row injection, broadcast-view rejection,
exact zero utilization, and bounded repeated-call and performance smoke behavior. When the native
metabolic extension is available, translated-batch parity uses `rtol=1e-12` and `atol=1e-15`.
Focused translator validation and repository-wide discovery are executable with:

```bash
python -m unittest bacteria.test_diet_translator
python -m unittest discover
```

`.github/workflows/kernel.yml` invokes `bacteria.test_diet_translator`,
`bacteria.test_metabolic_integration`, `bacteria.test_integration_stress`, and
`bacteria.test_comparative_sweeps` in both its OpenMP matrix and serial Linux job. Both jobs assert
that the built extension exposes `MetabolicKernel`. A configured command is not a claim that a
particular hosted workflow run passed, and no observed test count or timing is a versioned
acceptance contract.

`np.broadcast_to(increment, (E, 3))` is rejected because the result is a read-only, stride-zero
view rather than writable owning batch storage. An allocated writable, C-contiguous `(E, 3)` array
accepts row-wise additive injection and satisfies both Python and native batch contracts.

Native nonfinite-value rejection depends on IEEE-respecting floating-point behavior. The extension
build retains `-O3` on non-MSVC compilers and `/O2` on MSVC, removes `-ffast-math`, and uses
`/fp:precise` rather than `/fp:fast`. This preserves the finite checks exercised by the serial and
native parity tests.

### 10.5 Comparative Spatial-Control Sweep Harness

**Purpose and model boundary.** `run_comparative_sweeps.py` is a headless paired-seed harness for
the implemented legacy pg production model. It compares the categorical production diet profiles
`malnutrition`, `standard`, and `mdcf` under an explicit homogeneous control and the spatial-niche
control. It does not activate the reduced model, translator boluses, or
`GrowthOnlyFbaLifecyclePolicy`; its extinction records therefore do not validate live-FBA
extinction dynamics.<sup>[50](#ref-50)</sup>

The default configuration runs every diet for 1,000 updates per condition, with
$\Delta t=0.1\thinspace\mathrm{s}$, seed 0, Python backends, and the checked-in species data. A
typical multi-seed invocation that retains those other defaults is

```bash
python run_comparative_sweeps.py --seeds 0 1 2
```

The default output base is `/workspace/scratch`. By default, each invocation derives a fresh run
directory from its timestamped sweep ID; an existing name is never replaced, and explicit name
collisions receive a numeric suffix such as `-001`. Each CSV is opened in exclusive-create mode
inside that new directory.<sup>[50](#ref-50)</sup>

**Paired seed isolation.** For every diet/seed pair, the harness runs the homogeneous and spatial
conditions sequentially from independently reconstructed state. `isolated_rng` saves the caller's
Python and NumPy global RNG states, seeds both streams for one condition, and restores both states
on exit. Before accepting a pair, the harness requires equal hashes for the complete founder state
and for the post-founder RNG state. Its temporary spatial-global overrides are also restored after
each condition. This is paired stochastic control, not a claim that one seed demonstrates
robustness.<sup>[50](#ref-50)</sup>

The homogeneous control disables the partition, sets the mucus radius fraction and motility factor
to one, uses the lumen food-washout rate in both zones, deposits the full food-node amount, and uses
uniform carrying-capacity culling. The spatial-niche control enables the implemented
Section 4.20 settings: radius fraction $0.75$, mucus motility factor $0.3$, lumen and mucus
food-node washout rates $0.02\thinspace\mathrm{s}^{-1}$ and
$0.002\thinspace\mathrm{s}^{-1}$, mucus drop factor $0.5$, and lumen cull weight $3.0$.
These are **Own (BacSandbox model)** controls for testing the proposal's refuge hypothesis, not
measured gut-flow or persistence parameters.<sup>[48](#ref-48), [50](#ref-50)</sup>

**Presence/absence distance.** At each paired time point, let $A_t$ and $B_t$ be the sets of
species with positive counts in the homogeneous and spatial runs. The emitted Jaccard community
distance is<sup>[49](#ref-49)</sup>

$$
d_J(A_t,B_t)
=1-\frac{\lvert A_t\cap B_t\rvert}{\lvert A_t\cup B_t\rvert}.
$$

The implementation defines $d_J(\varnothing,\varnothing)=0$. Counts affect set membership only;
the abundance magnitudes do not enter this metric.<sup>[50](#ref-50)</sup>

**Extinction and right-censoring.** Extinction rows include only configured low-abundance founders,
defined by $0\lt N_s(0)\le L$, where the default threshold is $L=10$ cells. Let
$T_s^{\mathrm{ext}}$ be the first recorded zero crossing and $T_{\mathrm{end}}$ the final observed
time. The exported event indicator and analysis time are<sup>[16](#ref-16), [50](#ref-50)</sup>

$$
\delta_s=1,
\qquad
\widetilde{T}_s=T_s^{\mathrm{ext}}
\quad\text{when }T_s^{\mathrm{ext}}\text{ is observed by }T_{\mathrm{end}}.
$$

$$
\delta_s=0,
\qquad
\widetilde{T}_s=T_{\mathrm{end}}
\quad\text{otherwise}.
$$

Rows with $\delta_s=0$ are labeled `right_censored`; absence of an observed event is not recorded
as infinite survival or proof of persistence.

**Exact CSV schemas.** The four files have stable, ordered headers; tests reject drift and verify
that a second write cannot overwrite a prior run.<sup>[50](#ref-50)</sup>

`relative_abundance.csv` stores every species count and within-community proportion at every
timepoint in both conditions:

```text
sweep_id,pair_id,seed,diet,condition,step,elapsed_seconds,species,abundance_count,total_count,relative_abundance
```

`jaccard_distance.csv` stores one paired presence/absence comparison per timepoint:

```text
sweep_id,pair_id,seed,diet,step,elapsed_seconds,left_condition,right_condition,left_present_species_count,right_present_species_count,intersection_size,union_size,jaccard_community_distance
```

`extinction_censoring.csv` stores event/censoring records only for the predefined low-abundance
founders:

```text
sweep_id,pair_id,seed,diet,condition,species,initial_count,final_count,low_abundance_max_initial_count,low_abundance_definition,low_abundance_evidence_class,event_observed,status,time_to_event_seconds,censor_time_seconds
```

`run_metadata.csv` stores one reproducibility/provenance row per condition:

```text
sweep_id,pair_id,seed,diet,condition,dt_seconds,step_count,timepoint_count,simulation_engine,backend,physics_backend_used,ecm_backend_used,metabolic_model_id,metabolic_model_revision,metabolic_prior_id,metabolic_prior_revision,reduced_model_active,available_fba_lifecycle_policy,fba_extinction_dynamics_validated,fba_lifecycle_limitation,species_data_path,species_data_sha256,founder_state_sha256,post_founder_rng_state_sha256,world_width_um,world_height_um,carrying_capacity_cells,drop_frequency_s_inv,drop_range_min_pg,drop_range_max_pg,low_abundance_max_initial_count,low_abundance_definition,low_abundance_evidence_class,spatial_parameter_evidence_class,spatial_partition_enabled,mucus_layer_radius_fraction,mucus_motility_anchoring_factor,lumen_food_washout_rate_s_inv,mucus_food_washout_rate_s_inv,mucus_drop_amount_factor,lumen_competition_cull_weight,initial_spawn_mode,initial_species_arrangement,initial_spawn_lumen_bias,initial_spawn_lumen_radius_fraction,initial_spawn_cluster_count,initial_spawn_cluster_spread_fraction,initial_spawn_central_radius_fraction,petri_radius_padding_um,ecm_cell_size_um
```

The `backend` column records the requested selection. `physics_backend_used` and
`ecm_backend_used` record only a backend actually selected by successful dispatch; either field is
`not_exercised` when that subsystem never selected a backend, rather than repeating the requested
value.<sup>[50](#ref-50)</sup>

**Additive run manifest.** Every new sweep also writes `run_manifest.json` (schema 1)
beside the four unchanged CSVs. It records the species-file SHA-256, realized initial
counts, configured versus inoculated profile counts, initial harness environment,
actual backend selections, nominal and observed modeled-second horizons, and CSV
row counts. The harness constructs its own environment and does not load
`env_config.json`; pH, oxygen and carbon dioxide subsequently evolve under the
legacy rules. Initial values must not be described as clamped throughout the run.
Seeds and diets must be unique; integer-typed configuration values must be valid
integers and the modeled horizon must be finite.<sup>[50](#ref-50)</sup>

Supported NumPy real-valued configuration inputs are converted to built-in numbers
only in the manifest; simulation inputs retain their original numerical types.
The manifest is serialized before any output directory or CSV is created, so a
JSON serialization failure cannot leave an apparently complete CSV-only run.
Filesystem write failures can still leave incomplete outputs.<sup>[50](#ref-50)</sup>

The manifest explicitly records the harness's interpretive boundaries: no selected
anatomical region, effective physical depth/volume, empirical comparison,
experimental detection threshold, Kaplan–Meier curve calculation, or compartment
count export. Its event/censoring rows are inputs to possible later analysis.
Simulation seed variation is distinct from biological or measurement variation.
These fields describe new harness outputs; they are not retroactively present in
historical run packages.<sup>[50](#ref-50)</sup>

`bacteria.test_comparative_sweeps` verifies the 1,000-step and scratch-directory defaults, paired
RNG restoration and founder equivalence, distinct controls, Jaccard bounds, observed and censored
events, exact schemas, unexercised-backend metadata, non-overwriting writes, and the checked-in
affinity audit. It is configured in both the OpenMP and serial workflow jobs and can be run
directly with:

```bash
python -m unittest bacteria.test_comparative_sweeps
```

**Interpretive limitation.** In these production-model sweeps, lumen washout decreases and removes
free pg `FoodNode` mass. It does not advect, flush, or otherwise transport bacterial cells, and it
does not operate on the three live-FBA substrate maps. The benchmark produces analysis-ready
paired outputs, not a refuge result. A mucus-associated refuge remains a proposal hypothesis until
predefined multi-seed analyses demonstrate a robust difference against the homogeneous control.
<sup>[48](#ref-48), [50](#ref-50)</sup>

## Appendix A: Scientific References

Every biological or mathematical model, formula, named equation, scaling law, or empirically-derived
constant in this document carries an inline superscript marker (rendered as <sup>[N](#ref-N)</sup>)
that resolves to the correspondingly numbered entry below. Entries labeled **Own (BacSandbox model)**
are original formulations or calibrations with no external source. Every listed DOI was verified to
resolve; works without a registered DOI are marked *(no DOI)*.

<a id="ref-1"></a>**[1]** Balaban, Merrin, Chait, Kowalik & Leibler (2004). *Bacterial Persistence as a Phenotypic Switch*. Science, 305(5690):1622–1625. [https://doi.org/10.1126/science.1099390](https://doi.org/10.1126/science.1099390)

<a id="ref-2"></a>**[2]** Baranyi & Roberts (1994). *A dynamic approach to predicting bacterial growth in food*. International Journal of Food Microbiology, 23(3–4):277–294. [https://doi.org/10.1016/0168-1605(94)90157-0](https://doi.org/10.1016/0168-1605%2894%2990157-0)

<a id="ref-3"></a>**[3]** Berg & Brown (1972). *Chemotaxis in Escherichia coli analysed by three-dimensional tracking*. Nature, 239(5374):500–504. [https://doi.org/10.1038/239500a0](https://doi.org/10.1038/239500a0)

<a id="ref-4"></a>**[4]** Berg (1993). *Random Walks in Biology* (Expanded ed.). Princeton University Press. ISBN 978-0691000640 (no DOI)

<a id="ref-5"></a>**[5]** Courant, Friedrichs & Lewy (1928). *Über die partiellen Differenzengleichungen der mathematischen Physik*. Mathematische Annalen, 100(1):32–74. [https://doi.org/10.1007/BF01448839](https://doi.org/10.1007/BF01448839)

<a id="ref-6"></a>**[6]** Donaldson, Lee & Mazmanian (2016). *Gut biogeography of the bacterial microbiota*. Nature Reviews Microbiology, 14(1):20–32. [https://doi.org/10.1038/nrmicro3552](https://doi.org/10.1038/nrmicro3552)

<a id="ref-7"></a>**[7]** Efraimidis & Spirakis (2006). *Weighted random sampling with a reservoir*. Information Processing Letters, 97(5):181–185. [https://doi.org/10.1016/j.ipl.2005.11.003](https://doi.org/10.1016/j.ipl.2005.11.003)

<a id="ref-8"></a>**[8]** Ellis & Delbrück (1939). *The Growth of Bacteriophage*. Journal of General Physiology, 22(3):365–384. [https://doi.org/10.1085/jgp.22.3.365](https://doi.org/10.1085/jgp.22.3.365)

<a id="ref-9"></a>**[9]** Fick (1855). *Ueber Diffusion*. Annalen der Physik, 170(1):59–86. [https://doi.org/10.1002/andp.18551700105](https://doi.org/10.1002/andp.18551700105)

<a id="ref-10"></a>**[10]** Flint, Scott, Duncan, Louis & Forano (2012). *Microbial degradation of complex carbohydrates in the gut*. Gut Microbes, 3(4):289–306. [https://doi.org/10.4161/gmic.19897](https://doi.org/10.4161/gmic.19897)

<a id="ref-11"></a>**[11]** Gehrig, Venkatesh, Chang, Hibberd, Kung, Cheng, Chen, Subramanian, Cowardin, Meier et al. (2019). *Effects of microbiota-directed foods in gnotobiotic animals and undernourished children*. Science, 365(6449):eaau4732. [https://doi.org/10.1126/science.aau4732](https://doi.org/10.1126/science.aau4732)

<a id="ref-12"></a>**[12]** Gompertz (1825). *XXIV. On the nature of the function expressive of the law of human mortality, and on a new mode of determining the value of life contingencies. In a letter to Francis Baily, Esq. F. R. S. &c.* Philosophical Transactions of the Royal Society of London, 115:513–583. [https://doi.org/10.1098/rstl.1825.0026](https://doi.org/10.1098/rstl.1825.0026)

<a id="ref-13"></a>**[13]** Harman (1956). *Aging: a theory based on free radical and radiation chemistry*. Journal of Gerontology, 11(3):298–300. [https://doi.org/10.1093/geronj/11.3.298](https://doi.org/10.1093/geronj/11.3.298)

<a id="ref-14"></a>**[14]** Henderson (1908). *Concerning the relationship between the strength of acids and their capacity to preserve neutrality*. American Journal of Physiology, 21(2):173–179. [https://doi.org/10.1152/ajplegacy.1908.21.2.173](https://doi.org/10.1152/ajplegacy.1908.21.2.173)

<a id="ref-15"></a>**[15]** Kalinin, Jiang, Tu & Wu (2009). *Logarithmic sensing in Escherichia coli bacterial chemotaxis*. Biophysical Journal, 96(6):2439–2448. [https://doi.org/10.1016/j.bpj.2008.10.027](https://doi.org/10.1016/j.bpj.2008.10.027)

<a id="ref-16"></a>**[16]** Kaplan & Meier (1958). *Nonparametric estimation from incomplete observations*. Journal of the American Statistical Association, 53(282):457–481. [https://doi.org/10.1080/01621459.1958.10501452](https://doi.org/10.1080/01621459.1958.10501452)

<a id="ref-17"></a>**[17]** Kleiber (1932). *Body size and metabolism*. Hilgardia, 6(11):315–353. [https://doi.org/10.3733/hilg.v06n11p315](https://doi.org/10.3733/hilg.v06n11p315)

<a id="ref-18"></a>**[18]** Earle, Billings, Sigal, Lichtman, Hansson, Elias, Amieva, Huang & Sonnenburg (2015). *Quantitative Imaging of Gut Microbiota Spatial Organization*. Cell Host & Microbe, 18(4):478–488. [https://doi.org/10.1016/j.chom.2015.09.002](https://doi.org/10.1016/j.chom.2015.09.002)

<a id="ref-19"></a>**[19]** Levin, Stewart & Chao (1977). *Resource-Limited Growth, Competition, and Predation: A Model and Experimental Studies with Bacteria and Bacteriophage*. The American Naturalist, 111(977):3–24. [https://doi.org/10.1086/283134](https://doi.org/10.1086/283134)

<a id="ref-20"></a>**[20]** Lewis (2010). *Persister Cells*. Annual Review of Microbiology, 64:357–372. [https://doi.org/10.1146/annurev.micro.112408.134306](https://doi.org/10.1146/annurev.micro.112408.134306)

<a id="ref-21"></a>**[21]** Luria & Delbrück (1943). *Mutations of Bacteria from Virus Sensitivity to Virus Resistance*. Genetics, 28(6):491–511. [https://doi.org/10.1093/genetics/28.6.491](https://doi.org/10.1093/genetics/28.6.491)

<a id="ref-22"></a>**[22]** Ochman, Lawrence & Groisman (2000). *Lateral gene transfer and the nature of bacterial innovation*. Nature, 405(6784):299–304. [https://doi.org/10.1038/35012500](https://doi.org/10.1038/35012500)

<a id="ref-23"></a>**[23]** Ratkowsky, Lowry, McMeekin, Stokes & Chandler (1983). *Model for bacterial culture growth rate throughout the entire biokinetic temperature range*. Journal of Bacteriology, 154(3):1222–1226. [https://doi.org/10.1128/jb.154.3.1222-1226.1983](https://doi.org/10.1128/jb.154.3.1222-1226.1983)

<a id="ref-24"></a>**[24]** Raman, Gehrig, Venkatesh, Chang, Hibberd, Subramanian, Kang, Bessong, Lima, Kosek et al. (2019). *A sparse covarying unit that describes healthy and impaired human gut microbiota development*. Science, 365(6449):eaau4735. [https://doi.org/10.1126/science.aau4735](https://doi.org/10.1126/science.aau4735)

<a id="ref-25"></a>**[25]** Rich (2003). *The molecular machinery of Keilin's respiratory chain*. Biochemical Society Transactions, 31(6):1095–1105. [https://doi.org/10.1042/bst0311095](https://doi.org/10.1042/bst0311095)

<a id="ref-26"></a>**[26]** San Millán & MacLean (2017). *Fitness Costs of Plasmids: a Limit to Plasmid Transmission*. Microbiology Spectrum, 5(5):MTBP-0016-2017. [https://doi.org/10.1128/microbiolspec.MTBP-0016-2017](https://doi.org/10.1128/microbiolspec.MTBP-0016-2017)

<a id="ref-27"></a>**[27]** Sonnenburg & Sonnenburg (2014). *Starving our microbial self: the deleterious consequences of a diet deficient in microbiota-accessible carbohydrates*. Cell Metabolism, 20(5):779–786. [https://doi.org/10.1016/j.cmet.2014.07.003](https://doi.org/10.1016/j.cmet.2014.07.003)

<a id="ref-28"></a>**[28]** Stouthamer (1973). *A theoretical study on the amount of ATP required for synthesis of microbial cell material*. Antonie van Leeuwenhoek, 39(1):545–565. [https://doi.org/10.1007/BF02578899](https://doi.org/10.1007/BF02578899)

<a id="ref-29"></a>**[29]** Verhulst (1838). *Notice sur la loi que la population suit dans son accroissement*. Correspondance Mathématique et Physique, 10:113–121. (no DOI)

<a id="ref-30"></a>**[30]** Wolfe (2005). *The Acetate Switch*. Microbiology and Molecular Biology Reviews, 69(1):12–50. [https://doi.org/10.1128/MMBR.69.1.12-50.2005](https://doi.org/10.1128/MMBR.69.1.12-50.2005)

<a id="ref-31"></a>**[31]** Nelson & Cox (2017). *Lehninger Principles of Biochemistry* (7th ed.). W.H. Freeman. ISBN 978-1319117689 (no DOI)

<a id="ref-32"></a>**[32]** **Own (BacSandbox model).** Softened Euclidean nutrient-concentration proxy and gradient weighting (chemotaxis field regularization), §4.3.

<a id="ref-33"></a>**[33]** **Own (BacSandbox model).** Rendered-size starvation cue and volumetric nutrient-capacity scaling (geometric), §4.1.

<a id="ref-34"></a>**[34]** **Own (BacSandbox model).** Necromass recovery fractions (cytosol/dry/labile) and predation digestibility, §4.8.

<a id="ref-35"></a>**[35]** **Own (BacSandbox model).** Bacteriocin, contact-geometry, nutrient-sharing, T6SS damage/immunity, and phagocytosis formulations and constants, §4.12.

<a id="ref-36"></a>**[36]** **Own (BacSandbox model).** Age-dependent maintenance factor $A(t) = 1 + t/(180 \cdot L)$, §4.11.

<a id="ref-37"></a>**[37]** **Own (BacSandbox model).** Modeled dissolved-gas saturation capacities and exponential gas/pH relaxation rates, §4.10.

<a id="ref-38"></a>**[38]** **Own (BacSandbox model).** Initial-inoculum spawning policy (lumen-bias, clustered microcolonies) and elastic boundary reflection, §4.6/§4.23.

<a id="ref-39"></a>**[39]** **Own (BacSandbox model).** Floor-based world-to-sparse-cell mapping; fixed four-signal flat-buffer ABI; separate Python-transported three-substrate maps; lazy backend orchestration; and signal-layer parallel dispatch, §10.1.

<a id="ref-40"></a>**[40]** **Own (BacSandbox model).** Flat agent/species array interface, CSR ordering, staged OpenMP design, and six-counter event contract, §10.2.

<a id="ref-41"></a>**[41]** **Own (BacSandbox model).** Reduced-model interface and topology; flat/CSR ABI; deterministic posterior-mean constraint adapter; tolerance-neighbor cache/LRU policy and exact cache snapshots; adaptive substeps; live spatial configuration and registration-aware substrate snapshots; atomic retention-floor enforcement and complete weighted-bolus batch staging; state store; deterministic coordinate batching; normalized bolus deposition; represented-substrate and biomass checks; growth-only lifecycle; deferred transaction events; reset and serialization contracts; and full live-tick rollback, §10.3.

<a id="ref-42"></a>**[42]** Chen et al. (2021). *A Microbiota-Directed Food Intervention for Undernourished Children*. New England Journal of Medicine, 384:1517-1528. [https://doi.org/10.1056/NEJMoa2023294](https://doi.org/10.1056/NEJMoa2023294)

<a id="ref-43"></a>**[43]** Hibberd et al. (2024). *Bioactive glycans in a microbiome-directed food for children with malnutrition*. Nature, 625:157-165. [https://doi.org/10.1038/s41586-023-06838-3](https://doi.org/10.1038/s41586-023-06838-3)

<a id="ref-44"></a>**[44]** Chang et al. (2024). *Prevotella copri and microbiota members mediate the beneficial effects of a therapeutic food for malnutrition*. Nature Microbiology, 9:922-937. [https://doi.org/10.1038/s41564-024-01628-7](https://doi.org/10.1038/s41564-024-01628-7)

<a id="ref-45"></a>**[45]** Zhou et al. (2025). *Glycoside hydrolase-mediated glucomannan catabolism in Segatella copri, a target of microbiota-directed foods for malnourished children*. Proceedings of the National Academy of Sciences, 122(49):e2521522122. [https://doi.org/10.1073/pnas.2521522122](https://doi.org/10.1073/pnas.2521522122)

<a id="ref-46"></a>**[46]** National Institute of Standards and Technology (2026). *Glucose*. NIST Chemistry WebBook, SRD 69. Molecular weight $180.1559\thinspace\mathrm{g}/\mathrm{mol}$, rounded to $180.156\thinspace\mathrm{g}/\mathrm{mol}$ in the versioned asset. [https://webbook.nist.gov/cgi/cbook.cgi?ID=C50997&Units=SI](https://webbook.nist.gov/cgi/cbook.cgi?ID=C50997&Units=SI) *(no DOI)*

<a id="ref-47"></a>**[47]** **Own (BacSandbox model).** Versioned dietary-profile schema and bolus API; multi-substrate recipe, component, and routing translator; source-total normalization and represented-component mass conservation; total-starch allocation; sugar, oil, and micronutrient bookkeeping; nonclinical malnutrition baseline; molar-mass precision; species $V^{\max}$, $K$, and affinity-tier calibration; and translation QA fixtures, §10.4.

<a id="ref-48"></a>**[48]** **Own (BacSandbox model).** Production pg glucose-equivalent substrate classes, categorical nutrient-drop ratios, energy-yield calibrations, and dimensionless throughput/capacity rules, §4.19; lumen/mucus geometry, motility multiplier, food-washout rates, drop-retention factor, and carrying-capacity cull weighting used to construct a candidate-refuge hypothesis, §4.20.

<a id="ref-49"></a>**[49]** Real & Vargas (1996). *The Probabilistic Basis of Jaccard's Index of Similarity*. Systematic Biology, 45(3):380–385. [https://doi.org/10.1093/sysbio/45.3.380](https://doi.org/10.1093/sysbio/45.3.380)

<a id="ref-50"></a>**[50]** **Own (BacSandbox model).** Paired-seed legacy-pg comparative harness; RNG isolation and founder-state checks; homogeneous/spatial control construction; low-abundance definition; empty-community distance convention; CSV schemas and additive run manifest; actual-versus-unexercised backend metadata; configuration validation; metadata provenance; and non-overwriting output policy, §10.5.

<a id="ref-51"></a>**[51]** **Own (BacSandbox model).** Read-only species-parameter audit, §5.1: loader schema and compatibility-default/alias inventory, constructor normalization diagnostics, explicit-value gate, fixed-temperature evaluation, and non-overwriting JSON output. Configuration origin does not establish empirical provenance or biological validity.
