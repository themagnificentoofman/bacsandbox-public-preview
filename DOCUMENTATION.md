# BacSandbox: Formal Documentation

## 1. Introduction

BacSandbox is a micro-scale bacterial environment simulator built with Python and Pygame. It models independent bacterial species in a 2D petri-dish-like environment where cells move, sense organic nutrient gradients, consume resources, compete, reproduce, age, and die under configurable environmental pressures.

The simulation is real-time: biological rates are expressed per second and multiplied by elapsed `dt`. The `FPS` constant controls visual smoothness, not the amount of biology performed per update.

Food nodes represent an abstract organic nutrient pool measured in pg glucose-equivalent. This lets respiration, fermentation, uptake, and biomass synthesis share one energy/carbon accounting unit without implying that all environmental matter is literal glucose.

## 2. Setup And Usage

1. Install Python 3.
2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

3. Run the application:

   ```bash
   python -m bacteria
   ```

The app records population and event history in `simulation_stats.csv` for spreadsheet or external analysis workflows.

## 3. Core Units

* Distance: 1 world unit equals approximately 1 micrometer.
* Time: all simulation rates use elapsed seconds.
* Nutrient amount: 1 food unit equals 1 pg glucose-equivalent organic nutrient.
* Cell size: `base_size` is a structural cell size in micrometers. The rendered `size` only varies from 90% to 100% with fullness, so stored food does not inflate storage capacity.
* ECM concentration: quorum-sensing molecules and bacteriocins are stored on a coarse spatial grid and updated every simulation tick.
* Screen scale: `zoom` is pixels per micrometer in the Grid tab.

The Grid scale indicator uses the same transform as the renderer:

$$
x_{\text{screen}} = x_{\text{world}} z
$$

Therefore a bar of length $p$ pixels represents:

$$
d_{\mu\text{m}} = \frac{p}{z}
$$

The current UI displays this as `100 px = ...` next to the minimap.

## 4. Mathematical And Biological Models

### 4.1 Nutrient Storage And Cell Size

Nutrient capacity scales with structural volume, not temporary fullness<sup>[33](#ref-33)</sup>:
Nutrient capacity scales with structural volume, not temporary fullness<sup>[33](#ref-33)</sup>:

$$
C_{\text{food}} = k_{\text{cap}} \cdot s_{\text{base}}^3
$$

The rendered cell radius uses fullness only as a small starvation cue<sup>[33](#ref-33)</sup>:
The rendered cell radius uses fullness only as a small starvation cue<sup>[33](#ref-33)</sup>:

$$
s_{\text{rendered}} = s_{\text{base}} \left(0.9 + 0.1 \cdot \frac{F}{C_{\text{food}}}\right)
$$

This avoids the unrealistic behavior where eating more food directly increases the maximum storage cap.

### 4.2 Metabolism And Uptake

Basal metabolism follows an allometric Kleiber-style relationship<sup>[17](#ref-17)</sup>. Because cell mass scales with volume, a mass exponent of 3/4 becomes a size exponent of 2.25:
Basal metabolism follows an allometric Kleiber-style relationship<sup>[17](#ref-17)</sup>. Because cell mass scales with volume, a mass exponent of 3/4 becomes a size exponent of 2.25:

$$
B = \beta \cdot s_{\text{base}}^{2.25}
$$

Nutrient uptake is surface-area limited<sup>[33](#ref-33)</sup>:
Nutrient uptake is surface-area limited<sup>[33](#ref-33)</sup>:

$$
U_{\max} = u \cdot s_{\text{base}}^2
$$

During each elapsed-time update:

$$
F_{t+\Delta t} = F_t - B \cdot \Delta t \cdot A(t)
$$

where `A(t)` is the senescence maintenance factor.

### 4.3 Food Sensing And Chemotaxis

Bacteria do not teleport or directly seek a target. They sense local nutrient fields and bias their swimming direction and run duration. Each food node contributes a softened Euclidean concentration proxy<sup>[32](#ref-32)</sup>:
Bacteria do not teleport or directly seek a target. They sense local nutrient fields and bias their swimming direction and run duration. Each food node contributes a softened Euclidean concentration proxy<sup>[32](#ref-32)</sup>:

$$
C_i = \frac{A_i}{\pi(d_i^2 + a^2)}
$$

where $A_i$ is node amount, $d_i$ is Euclidean distance, and $a = 5.0\,\mu\text{m}$ prevents a zero-distance singularity.

Chemoreceptor response is log-compressed to approximate saturation<sup>[15](#ref-15)</sup>:
Chemoreceptor response is log-compressed to approximate saturation<sup>[15](#ref-15)</sup>:

$$
S_i = g_{\text{eff}} \log\left(1 + \frac{C_i}{C_{\text{ref}}}\right)
$$

where $C_{\text{ref}} = 10^{-3}\,\text{pg}/\mu\text{m}^2$ and $g_{\text{eff}}$ increases as the bacterium becomes hungry.

The local gradient direction is a weighted sum of node direction vectors<sup>[32](#ref-32)</sup>:
The local gradient direction is a weighted sum of node direction vectors<sup>[32](#ref-32)</sup>:

$$
\vec{G} = \sum_i S_i \frac{\vec{x}_i - \vec{x}_{\text{cell}}}{\lVert \vec{x}_i - \vec{x}_{\text{cell}} \rVert}
$$

### 4.4 Run-And-Tumble Motion

BacSandbox uses an E. coli-inspired run-and-tumble model<sup>[3](#ref-3)</sup>. Runs are sampled around 1 second and clamped to biological bounds:
BacSandbox uses an E. coli-inspired run-and-tumble model<sup>[3](#ref-3)</sup>. Runs are sampled around 1 second and clamped to biological bounds:

$$
t_{\text{run}} \sim \mathrm{clamp}(\mathcal{N}(1.0, 0.2^2), 0.25, 3.0)
$$

Tumbles last about 0.1 seconds:

$$
t_{\text{tumble}} = 0.1\,\text{s}
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
During a tumble, the cell performs one event-scale reorientation rather than a tiny frame-scaled wobble. The unbiased magnitude is centered near the classic observed E. coli reorientation angle of 68 degrees<sup>[3](#ref-3)</sup>:

$$
|\Delta \theta_{\text{tumble}}| \sim \mathcal{N}(1.19, 0.45^2)\,\text{rad}
$$

Translation during tumbles is multiplied by 0.05, representing motor activity being spent on reorientation rather than propulsion. Speed is never negative; apparent reversal is a new heading after a tumble, not backward swimming.

### 4.5 Rotational Diffusion

Brownian rotational noise follows the low-Reynolds scaling for small particles in fluid<sup>[4](#ref-4)</sup>:
Brownian rotational noise follows the low-Reynolds scaling for small particles in fluid<sup>[4](#ref-4)</sup>:

$$
D_{\text{rot}} \propto \frac{1}{s_{\text{base}}^3}
$$

The code uses a reference value of about `0.16 rad^2/s` for a 1 micrometer sphere in water.

Because $D_{\text{rot}}$ is an angular variance rate, the sampled angular standard deviation is:

$$
\sigma_{\theta} = \sqrt{2D_{\text{rot}}\Delta t}
$$

This keeps Brownian heading noise tied to elapsed seconds instead of tick count.

### 4.6 Boundary Reflection

The petri dish is circular. If a bacterium reaches the dish boundary, it is moved slightly inside the rim and its velocity heading reflects off the outward normal only when it was moving outward<sup>[38](#ref-38)</sup>:
The petri dish is circular. If a bacterium reaches the dish boundary, it is moved slightly inside the rim and its velocity heading reflects off the outward normal only when it was moving outward<sup>[38](#ref-38)</sup>:

$$
\vec{v}' = \vec{v} - 2(\vec{v} \cdot \hat{n})\hat{n}
$$

This is applied only when:

$$
\vec{v} \cdot \hat{n} > 0
$$

### 4.7 Reproduction

Population-level reproduction probability follows Verhulst logistic pressure through carrying capacity<sup>[29](#ref-29)</sup>:
Population-level reproduction probability follows Verhulst logistic pressure through carrying capacity<sup>[29](#ref-29)</sup>:

$$
\frac{dN}{dt} = rN\left(1 - \frac{N}{K}\right)
$$

At the individual level, binary fission requires stored nutrient above the reproduction threshold and a synthesis cost inspired by Stouthamer-style ATP/biomass costs<sup>[28](#ref-28)</sup>:
At the individual level, binary fission requires stored nutrient above the reproduction threshold and a synthesis cost inspired by Stouthamer-style ATP/biomass costs<sup>[28](#ref-28)</sup>:

$$
F_{\text{cost}} = 40.0 \cdot s_{\text{base}}^3
$$

The parent divides only if both parent and child retain enough nutrient to survive the next update.

The implementation exposes this conversion as `FISSION_COST_PG_PER_UM3 = 40.0` so tests and model code share the same biological assumption.

### 4.8 Death, Predation, And Necromass Recycling

Dead bacteria do not become pure glucose. Death releases a lossy glucose-equivalent nutrient value from cytosolic leakage plus labile dry necromass<sup>[34](#ref-34)</sup>:
Dead bacteria do not become pure glucose. Death releases a lossy glucose-equivalent nutrient value from cytosolic leakage plus labile dry necromass<sup>[34](#ref-34)</sup>:

$$
F_{\text{recovered}} = \max(0, F_{\text{stored}}) f_{\text{cytosol}} + V_{\text{cell}} m_{\text{wet}} f_{\text{dry}} f_{\text{labile}}
$$

The implemented constants are:

$$
f_{\text{cytosol}} = 0.75, \quad m_{\text{wet}} = 1.0\,\text{pg}/\mu\text{m}^3, \quad f_{\text{dry}} = 0.30, \quad f_{\text{labile}} = 0.35
$$

Most structural biomass is therefore unavailable, refractory, respired, or lost to the environment. Predation uses the same digestible nutrient logic instead of consuming a fixed fraction of the prey's storage capacity.

### 4.9 ECM Diffusion Grid

The extracellular matrix tracks quorum-sensing molecules, bacteriocins, and cross-feeding acetate (Section 4.15) on a spatial hash grid. ECM diffusion now processes every active grid cell each simulation tick. It uses an explicit Fickian 2D stencil with a stability clamp:
The extracellular matrix tracks quorum-sensing molecules, bacteriocins, and cross-feeding acetate (Section 4.15) on a spatial hash grid. ECM diffusion now processes every active grid cell each simulation tick. It uses an explicit Fickian 2D stencil<sup>[9](#ref-9)</sup> with a stability clamp<sup>[5](#ref-5)</sup>:
The extracellular matrix tracks quorum-sensing molecules, bacteriocins, and cross-feeding acetate (Section 4.15) on a spatial hash grid. ECM diffusion now processes every active grid cell each simulation tick. It uses an explicit Fickian 2D stencil<sup>[9](#ref-9)</sup> with a stability clamp<sup>[5](#ref-5)</sup>:

$$
\alpha = \min\left(\alpha_{\max}, \frac{D\Delta t}{\Delta x^2}\right), \quad \alpha_{\max}=0.24
$$

Decay is applied as:

$$
m = \mathrm{clamp}(1 - \lambda\Delta t, 0, 1)
$$

The finite-difference update is:

$$
C'_{i,j}=m\left((1-4\alpha)C_{i,j}+\alpha(C_{i-1,j}+C_{i+1,j}+C_{i,j-1}+C_{i,j+1})\right)
$$

The simulation uses a 20 micrometer ECM cell size. Quorum molecules use a higher diffusion coefficient than bacteriocins because the latter are larger antimicrobial peptides/proteins.

Concentrations below $10^{-3}$ are pruned from the sparse grid. This threshold is numerical housekeeping: those values are below the simulator's biological resolution and would otherwise accumulate as dead sparse-grid entries.

### 4.10 Gas Exchange, Fermentation, And pH

Aerobic and anaerobic environmental effects are based on each cell's actual size- and age-scaled nutrient consumption during the current update. Each bacterium stores `last_metabolic_consumption_pg` immediately after its metabolism step:

$$
G_{\Delta t}=\sum_i \min(F_i, B_i\Delta t A_i)
$$

where $G_{\Delta t}$ is consumed glucose-equivalent nutrient, $F_i$ is the cell's available stored nutrient before metabolism, $B_i$ is its size-scaled basal metabolism, and $A_i$ is the age maintenance factor.

For aerobic species with positive oxygen sensitivity $s_{O_2}$, respiration follows<sup>[31](#ref-31)</sup>:
For aerobic species with positive oxygen sensitivity $s_{O_2}$, respiration follows<sup>[31](#ref-31)</sup>:

$$
\text{C}_6\text{H}_{12}\text{O}_6 + 6\text{O}_2 \rightarrow 6\text{CO}_2 + 6\text{H}_2\text{O}
$$

The implementation converts glucose-equivalent pg into dissolved gas fractions with molar mass ratios and model capacities<sup>[37](#ref-37)</sup>:
The implementation converts glucose-equivalent pg into dissolved gas fractions with molar mass ratios and model capacities<sup>[37](#ref-37)</sup>:

$$
\Delta O_2 = \frac{G_{\Delta t}\left(\frac{192.0}{180.16}\right)s_{O_2}}{4.0\times10^7}
$$

$$
\Delta CO_2 = \frac{G_{\Delta t}\left(\frac{264.0}{180.16}\right)s_{O_2}}{1.5\times10^9}
$$

Carbonic acid contributes a small pH decrement:

$$
\Delta \text{acid}_{CO_2}=0.02\Delta CO_2
$$

For anaerobic species with negative oxygen sensitivity, the simulator uses a simplified lactic-acid pathway plus trace heterolactic carbon dioxide<sup>[31](#ref-31)</sup>:
For anaerobic species with negative oxygen sensitivity, the simulator uses a simplified lactic-acid pathway plus trace heterolactic carbon dioxide<sup>[31](#ref-31)</sup>:

$$
\Delta CO_{2,\text{anaerobic}} = \frac{0.1G_{\Delta t}|s_{O_2}|}{1.5\times10^9}
$$

$$
\Delta \text{acid}_{\text{lactic}} = \frac{2.5G_{\Delta t}|s_{O_2}|}{4.0\times10^7}
$$

The update sequence is: subtract oxygen use, add carbon dioxide production, subtract acid from pH, equilibrate gases toward an open dish<sup>[9](#ref-9)</sup>, then buffer pH toward neutral<sup>[14](#ref-14)</sup>. Gas and buffer restoration are<sup>[37](#ref-37)</sup>:
The update sequence is: subtract oxygen use, add carbon dioxide production, subtract acid from pH, equilibrate gases toward an open dish<sup>[9](#ref-9)</sup>, then buffer pH toward neutral<sup>[14](#ref-14)</sup>. Gas and buffer restoration are<sup>[37](#ref-37)</sup>:

$$
O_{2,t+\Delta t}=O_{2,t}+(0.21-O_{2,t})k_{\text{gas}}\Delta t
$$

$$
CO_{2,t+\Delta t}=CO_{2,t}+(0-CO_{2,t})k_{\text{gas}}\Delta t
$$

$$
pH_{t+\Delta t}=pH_t+(7.0-pH_t)k_{\mathrm{buffer}}\Delta t
$$

### 4.11 Senescence

Maintenance cost rises with age<sup>[36](#ref-36)</sup>:
Maintenance cost rises with age<sup>[36](#ref-36)</sup>:

$$
A(t) = 1 + \frac{t}{180 \cdot L}
$$

Intrinsic death risk follows a Gompertz-style exponential penalty<sup>[12](#ref-12), [13](#ref-13)</sup>:
Intrinsic death risk follows a Gompertz-style exponential penalty<sup>[12](#ref-12), [13](#ref-13)</sup>:

$$
\mu(t) = r_{\text{aging}} \Delta t \cdot e^{\alpha t}
$$

### 4.12 Competition, Bacteriocins, And T6SS

Diffusible bacteriocin production consumes stored nutrient and emits antimicrobial material into the ECM<sup>[35](#ref-35)</sup>:
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
Lag phase (Baranyi & Roberts dynamic state, see Section 4.16)<sup>[2](#ref-2)</sup>:

* `DEFAULT_BARANYI_Q0 = 0.05` (fresh-inoculum physiological state; $h_0 = \ln(1 + 1/q_0) \approx 3.04$)
* `BARANYI_Q_STATE_MAX` $= 1.0 \times 10^{6}$ (numerical cap where $\alpha > 0.999999$)
* `BARANYI_Q0_MIN` $= 1.0 \times 10^{-4}$ (floor; guarantees a finite lag)
* `BARANYI_Q0_MAX` $= 1.0 \times 10^{3}$ (ceiling; effectively disables the lag)

Phenotypic switching (persister bet-hedging, see Section 4.17)<sup>[1](#ref-1), [20](#ref-20)</sup>:
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
derived from the four-parameter Ratkowsky square-root model<sup>[23](#ref-23)</sup>, which is biologically validated across the entire biokinetic temperature range:

$$
\sqrt{\mu(T)} = b\,(T - T_{min})\,\bigl(1 - e^{\,c\,(T - T_{max})}\bigr)
$$

where $T$ is the environment temperature in degrees Celsius, $T_{min}$ and $T_{max}$ are the
species' theoretical minimum and maximum growth temperatures (its *cardinal* temperatures),
$b$ is the low-temperature square-root slope, and $c$ governs the steepness of high-temperature
thermal inactivation. Outside the open interval $(T_{min}, T_{max})$ growth is zero.

The raw rate $\mu(T) = \bigl(\sqrt{\mu(T)}\bigr)^2$ is converted into a dimensionless multiplier in
$[0, 1]$ by normalizing against each species' own optimum $\sqrt{\mu_{\text{peak}}}$, located once at
construction by uniform sampling across $[T_{min}, T_{max}]$:

$$
f_{T}(T) = \left( \frac{\sqrt{\mu(T)}}{\sqrt{\mu_{\text{peak}}}} \right)^{2}, \qquad f_{T} \in [0, 1]
$$

Because the multiplier is normalized to the optimum, the slope parameter $b$ cancels and only the
curve *shape* (set by $T_{min}$, $T_{max}$, and $c$) matters. The factor reaches $1.0$ at the
species optimum and falls to $0.0$ at its cardinal limits. It is applied as follows:

* **Reproduction:** $f_{T}$ multiplies the intrinsic Verhulst growth rate $r$ alongside the
  carrying-capacity term and environmental fitness (see Section 4.7).
* **Metabolism:** basal nutrient demand is scaled by $f_{T}$, so chilled cultures deplete their
  nutrient reserves far more slowly, reproducing the long shelf life of refrigerated cells.
* **Nutrient uptake:** membrane transporter throughput (`consumption_speed`) is scaled by $f_{T}$.
* **Motility:** flagellar swimming speed is scaled by $f_{T}$, since motor output is enzyme-driven.

Near and beyond the cardinal limits $f_{T} \to 0$, driving cells toward thermal dormancy rather than
immediate death. Different cardinal temperatures yield distinct thermal niches: a psychrophile
(low $T_{min}$, $T_{max}$) peaks at a colder optimum than a mesophile, enabling temperature-based
niche separation between species.


### 4.15 Metabolic Cross-Feeding (Acetate Overflow)

Real microbial communities are knit together by *cross-feeding*: the metabolic waste of one
organism is the food of another. The most studied example in *Escherichia coli* is the **acetate
switch**<sup>[30](#ref-30)</sup>: during fast aerobic growth on
switch**<sup>[30](#ref-30)</sup>: during fast aerobic growth on
glucose the cell's central metabolism outpaces its respiratory capacity and excretes *overflow*
acetate; when glucose runs low (or in neighbouring taxa) that acetate is re-imported and oxidized
as a secondary carbon source. BacSandbox models acetate as a third diffusible ECM layer (alongside
quorum-sensing inducers and bacteriocins) coupling acetate **producers** to acetate **consumers**.

**Secretion.** A producer with secretion trait $s \in [0, 1]$ converts a fraction of the glucose it
took up on the previous tick, $U_{\text{glc}}$ (pg), into excreted acetate, bounded by the
stoichiometric ceiling of one glucose splitting into two acetate molecules by mass:

$$
A_{\text{sec}} = \min\!\Bigl(\, s \, U_{\text{glc}} \, Y_{a/g},\; F \,\Bigr), \qquad
Y_{a/g} = \frac{2\,M_{\text{ace}}}{M_{\text{glc}}} = \frac{2 \times 60.05}{180.16} \approx 0.667
$$

where $F$ is the cell's current internal nutrient store (acetate cannot be excreted from carbon the
cell does not hold) and $M_{\text{ace}}, M_{\text{glc}}$ are the molar masses of acetate and glucose.
The secreted mass is removed from the cell's store and deposited into the local acetate ECM cell.

**Uptake.** A consumer with uptake trait $u \in [0, 1]$ oxidizes locally available acetate at a rate
limited by its transporter budget, by local availability $A_{\text{avail}}$, and by remaining storage
room. Because acetate oxidation is aerobic and enzyme-driven, uptake is gated by both the species'
Ratkowsky temperature factor $f_T$ (Section 4.14) and a dissolved-oxygen gate:

$$
g_{O_2} = \min\!\left(1,\; \frac{[\mathrm{O_2}]}{[\mathrm{O_2}]_{\text{sat}}}\right), \qquad
B = c_{\text{eat}} \, u \, \Delta t \, f_T \, g_{O_2}
$$

$$
A_{\text{con}} = \max\!\Bigl(0,\; \min\bigl(A_{\text{avail}},\; B,\; R \bigr)\Bigr), \qquad
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

so $\Delta F = A_{\text{con}} \, Y_E$. The ATP yields ($\sim 32$ ATP/glucose<sup>[25](#ref-25)</sup>; $\sim 10$
so $\Delta F = A_{\text{con}} \, Y_E$. The ATP yields ($\sim 32$ ATP/glucose<sup>[25](#ref-25)</sup>; $\sim 10$
ATP/acetate via acetyl-CoA → TCA cycle + electron transport chain) make acetate a slightly less
energy-dense carbon source per unit mass, as observed experimentally.

**Gas coupling.** Acetate oxidation follows $\mathrm{C_2H_4O_2 + 2\,O_2 \rightarrow 2\,CO_2 + 2\,H_2O}$<sup>[31](#ref-31)</sup>.
**Gas coupling.** Acetate oxidation follows $\mathrm{C_2H_4O_2 + 2\,O_2 \rightarrow 2\,CO_2 + 2\,H_2O}$<sup>[31](#ref-31)</sup>.
Each picogram of acetate consumed therefore draws down dissolved $\mathrm{O_2}$ and releases
$\mathrm{CO_2}$ (with a corresponding carbonic-acid pH effect) independently of the glucose-driven
respiration in Section 4.10:

$$
\Delta[\mathrm{O_2}] = -\,A_{\text{con}}\,\frac{64.0}{M_{\text{ace}}}\Big/ C_{O_2}, \qquad
\Delta[\mathrm{CO_2}] = +\,A_{\text{con}}\,\frac{88.02}{M_{\text{ace}}}\Big/ C_{CO_2}
$$

where $64.0$ g and $88.02$ g are the masses of 2 mol $\mathrm{O_2}$ and 2 mol $\mathrm{CO_2}$ per mol
acetate, and $C_{O_2}, C_{CO_2}$ are the modeled dissolved-gas saturation capacities.

**Numerical note.** Picogram acetate masses are scaled into ECM grid units by a constant factor
(`ACETATE_ECM_UNITS_PER_PG`) so traces survive the sparse-grid pruning floor; the identical factor is
inverted on uptake so secreted and consumed mass are conserved. The acetate layer diffuses on the
same explicit Fickian stencil as the other ECM layers (Section 4.9). With the default traits
($s = u = 0$) cross-feeding is inert and pre-M2 behaviour is reproduced exactly.


### 4.16 Lag Phase (Baranyi & Roberts Dynamic Physiological State)

Newly inoculated cells do not divide immediately: they must first synthesize the ribosomes, enzymes,
and membrane components needed for growth. This delay is the **lag phase**. BacSandbox models it with
the dynamic approach of Baranyi & Roberts<sup>[2](#ref-2)</sup>, which augments each cell with a dimensionless
the dynamic approach of Baranyi & Roberts<sup>[2](#ref-2)</sup>, which augments each cell with a dimensionless
*physiological state* variable $q$ that captures the "work still to be done" before division.

The realized growth rate is the maximum (temperature-scaled) rate multiplied by an **adjustment
function** that depends only on $q$:

$$\alpha(q) = \frac{q}{1 + q}$$

When $q$ is small (fresh inoculum) $\alpha \approx 0$ and reproduction is suppressed; as the cell
adapts $q \to \infty$ and $\alpha \to 1$, restoring the full rate. The physiological state itself
grows in proportion to the specific growth rate,

$$\frac{dq}{dt} = \mu_{max}\, \phi(T)\, q,$$

where $\mu_{max}$ is the species' base reproduction rate (`rep`) and $\phi(T)$ is the Ratkowsky
temperature factor (Section 4.14). Because the right-hand side is linear in $q$, the equation is
integrated **analytically** over each tick rather than with forward Euler, so the result is exact and
frame-rate independent:

$$q_{t+\Delta t} = q_t \, e^{\mu_{max}\,\phi(T)\,\Delta t}.$$

This temperature coupling reproduces the observed fact that cold cultures lag far longer than
cultures at their thermal optimum, since $\phi(T) \to 0$ slows adaptation.

**Lag duration.** Integrating the adjustment function gives the classic Baranyi lag time

$$\lambda = \frac{h_0}{\mu_{max}}, \qquad h_0 = \ln\!\left(1 + \frac{1}{q_0}\right),$$

where $q_0$ is the initial physiological state (`baranyi_q0`) and $h_0$ is the dimensionless "work to
be done". A smaller $q_0$ means a longer lag. The reproduction probability per tick (Section 4.14)
therefore becomes

$$P_{rep} = (\text{rep} \cdot \Delta t)\, (1 - N/K)\, f_{env}\, \phi(T)\, \alpha(q).$$

**Inheritance and reset.** Only the **initial inoculum** starts at $q_0$; daughter cells produced by
fission inherit the parent's already-elevated $q$, so an actively growing population does not re-enter
lag on every division. A fresh inoculum (or a culture seeded into a new patch via reset) therefore
shows the characteristic flat-then-exponential growth curve.

**Numerical bounds.** $q_0$ is clamped to $[\,10^{-4}, 10^{3}\,]$ so the lag is always finite and
non-infinite, and the running state $q$ is capped at `BARANYI_Q_STATE_MAX` $= 10^6$ (where
$\alpha > 0.999999$) to prevent floating-point overflow on long runs. Setting `baranyi_q0` to its
ceiling makes $\alpha \approx 1$ from the outset, reproducing the pre-M3 (no-lag) growth behaviour.


### 4.17 Phenotypic Switching (Persister Bet-Hedging)

Even in clonal, genetically identical populations a small fraction of cells stochastically enter a
**persister** state: a dormant, slow-growing, stress-tolerant phenotype that survives antimicrobial
exposure lethal to actively dividing cells. When conditions improve the persisters resuscitate back
to the proliferative state. This reversible bistable switch is the canonical explanation for
antibiotic *tolerance* (as distinct from inherited resistance)<sup>[1](#ref-1), [20](#ref-20)</sup>. BacSandbox models it as a two-state
antibiotic *tolerance* (as distinct from inherited resistance)<sup>[1](#ref-1), [20](#ref-20)</sup>. BacSandbox models it as a two-state
Markov process whose transition rates depend on the cell's internal **fullness**
$f = \min(1, \max(0, \text{food}/\text{capacity}))$, used as a self-consistent proxy for sustained
local nutrient availability.

**Entry (proliferative $\rightarrow$ persister).** Below the starvation threshold $S_1$ the entry
rate combines a constant bet-hedging term with a stress term that grows as reserves empty:

$$
k^{+}(f) = k_0^{+} + \mathbb{1}[f < S_1]\, k_s^{+}\left(1 - \frac{f}{\max(S_1, \varepsilon)}\right)
$$

**Exit (persister $\rightarrow$ proliferative).** Above the recovery threshold $S_2$ the resuscitation
rate combines a constant reversion term with a recovery term that grows as reserves refill:

$$
k^{-}(f) = k_0^{-} + \mathbb{1}[f > S_2]\, k_r^{-}\,\frac{f - S_2}{\max(1 - S_2, \varepsilon)}
$$

**Per-tick probability.** Each transition is realized with the exact exponential survival
complement over the real elapsed time $\Delta t$, so the dynamics are frame-rate independent:

$$
P_{\text{switch}} = 1 - e^{-\,k\,\Delta t}.
$$

The entry rate is additionally scaled by the per-species propensity $s \in [0,1]$
(`persister_switch`), i.e. $k_{\text{eff}}^{+} = s\,k^{+}(f)$, so $s = 0$ disables persistence
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
and keeps the model lightweight while reproducing the classic predator–prey boom–bust oscillation<sup>[19](#ref-19)</sup>.

**Adsorption (susceptible $\rightarrow$ infected).** A susceptible (non-resistant, uninfected) cell
adsorbs phage by mass action against the local free-virion concentration $[P]$ sampled at its
position. Over the real elapsed time $\Delta t$ the infection probability is the exponential survival
complement, so the dynamics are frame-rate independent:

$$
P_{\text{adsorb}} = 1 - e^{-\,k\,[P]\,\Delta t},
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
structure (a fixed latent period followed by a discrete burst) follows Ellis & Delbrück<sup>[8](#ref-8)</sup>.

**Free-virion field.** Released virions diffuse (`ECM_PHAGE_DIFFUSION_UM2_S`) and lose infectivity
(`ECM_PHAGE_DECAY_S`) via the same Fickian stencil used for the other ECM layers (Section 4.9),
spreading the epidemic spatially and decaying it where hosts are absent.

**Resistance, mutation, and rebound.** Each cell is born resistant with probability equal to its
species' `phage_resistance` trait (a heritable, chromosomal fraction). Resistant cells never adsorb
phage, so they survive the crash. On every division a susceptible lineage can spontaneously mutate to
resistant with probability `PHAGE_RESISTANCE_MUTATION_RATE` (a Luria–Delbrück-style fluctuation)<sup>[21](#ref-21)</sup>,
resistant with probability `PHAGE_RESISTANCE_MUTATION_RATE` (a Luria–Delbrück-style fluctuation)<sup>[21](#ref-21)</sup>,
seeding the post-crash recovery. Daughters inherit the parent's resistance state and carriage cost.
The net result is the textbook sequence: phage outbreak $\rightarrow$ susceptible population crash
$\rightarrow$ resistant-only survivors $\rightarrow$ rebound.

**Horizontal gene transfer (HGT).** Resistance also spreads laterally<sup>[22](#ref-22)</sup>. When a resistant donor and a
**Horizontal gene transfer (HGT).** Resistance also spreads laterally<sup>[22](#ref-22)</sup>. When a resistant donor and a
susceptible recipient are in physical contact (the same contact test used for sharing/T6SS), and both
species are HGT-competent (`hgt_competence` $> 0$), the trait transfers by a conjugation-like process
with probability

$$
P_{\text{HGT}} = 1 - e^{-\,k_H\, c_{\text{donor}}\, c_{\text{recip}}\,\Delta t},
$$

where $k_H$ is `HGT_TRANSFER_RATE_PER_SEC` and $c$ are the two cells' competences. The recipient
becomes resistant and inherits a **plasmid carriage cost**: its basal maintenance metabolism is
multiplied by $1 +$ `HGT_METABOLIC_BURDEN`, modelling the well-documented fitness penalty of
maintaining a resistance plasmid<sup>[26](#ref-26)</sup>. This cost is inherited by daughters, so
maintaining a resistance plasmid<sup>[26](#ref-26)</sup>. This cost is inherited by daughters, so
in the absence of phage pressure the burdened lineage is slightly outcompeted — capturing the
cost–benefit tension that governs the spread of mobile resistance elements.

**Outbreak control.** The global `phage_inflow` environment parameter seeds free virions via a
Poisson process ($\lambda\,\Delta t$ packets of `PHAGE_SEED_PACKET_VIRIONS` at random petri
locations), so a user can trigger or halt an epidemic live from the Globals tab. A value of $0$
(the default) means no outbreak and exactly reproduces the pre-M5 behaviour.


### 4.19 Diet-Driven Substrate Metabolism (Phase I Malnutrition Model)

To ground BacSandbox in malnutrition research<sup>[27](#ref-27)</sup>, the single abstract glucose pool is replaced by
To ground BacSandbox in malnutrition research<sup>[27](#ref-27)</sup>, the single abstract glucose pool is replaced by
**three discrete dietary substrate classes**, each nutrient patch (`FoodNode`) carrying a
`substrate_id` $s \in \{0, 1, 2\}$:

- $s = 0$ — **generic calorie**: a simple, universally fermentable sugar (e.g. rice/starch hydrolysate).
- $s = 1$ — **complex glycan**: dietary fibre only partially fermentable by the gut community.
- $s = 2$ — **MDCF prebiotic**: a Microbiota-Directed Complementary Food substrate that selectively rewards competent taxa.

Necromass, decomposition, and all legacy food default to $s = 0$, so simulations with no diet
configured reproduce the original dynamics exactly.

**Per-substrate affinity.** Each species declares an affinity vector
$\mathbf{a} = (a_0, a_1, a_2)$ with $a_s \ge 0$ (the columns `aff_calorie`, `aff_glycan`,
`aff_mdcf`). Affinity multiplies the cell's transporter throughput for that substrate class; an
affinity of $0$ means the substrate is metabolically inaccessible. Chemotaxis is likewise
selective — cells are not attracted to substrates they cannot consume — by weighting the
chemotactic score with $a_s$:

$$\text{score}_s = g_{\text{eff}} \, a_s \, \ln\!\left(1 + \frac{C}{C_{\text{ref}}}\right)$$

where $C$ is the local patch concentration and $C_{\text{ref}}$ the receptor reference
concentration. Substrates with $a_s = 0$ are skipped entirely in both targeting and feeding.

**Affinity-gated uptake with energy yield.** Within one time step a cell has a shared transporter
budget $U_{\max} = v_c \, \Delta t \, \phi_T$ (mass, pg), where $v_c$ is `consumption_speed` and
$\phi_T$ the temperature factor. For each reachable patch of class $s$, the ingested **mass** is

$$\Delta m = \min\!\big(R \, a_s,\; R,\; m_{\text{node}},\; h_s\big), \qquad R = U_{\max} - \textstyle\sum \Delta m_{\text{prior}}$$

with $R$ the remaining budget and $h_s$ the storage headroom expressed back in substrate mass.
Mass is conserved on the patch ($m_{\text{node}} \leftarrow m_{\text{node}} - \Delta m$), but the
**energy** added to the cell's internal store is scaled by a substrate-specific yield $Y_s$:

$$\Delta E = \Delta m \, Y_s$$

The yields are fractions of pure glucose energy: $Y_0 = 1.0$ (simple sugar, fully oxidised),
$Y_1 = 0.7$ (complex glycan, energy lost to incomplete fermentation and indigestible residue),
$Y_2 = 0.85$ (MDCF prebiotic, a targeted, efficiently fermentable substrate). All yields satisfy
$0 < Y_s \le 1$. This makes fibre and MDCF intrinsically less calorically dense than refined sugar
while still selectively favouring the taxa equipped to use them<sup>[10](#ref-10)</sup>.
while still selectively favouring the taxa equipped to use them<sup>[10](#ref-10)</sup>.

**Diet configuration profiles.** The environment's `diet_profile` selects a categorical
probability distribution over substrate classes for each fresh nutrient drop. Three feeding arms
mirror the MDCF trial design<sup>[11](#ref-11), [24](#ref-24)</sup>; the per-profile fractions sum to $1$:
mirror the MDCF trial design<sup>[11](#ref-11), [24](#ref-24)</sup>; the per-profile fractions sum to $1$:

| Profile | calorie $p_0$ | glycan $p_1$ | MDCF $p_2$ |
|---|---|---|---|
| `malnutrition` | $0.85$ | $0.10$ | $0.05$ |
| `standard` | $0.60$ | $0.30$ | $0.10$ |
| `mdcf` | $0.40$ | $0.25$ | $0.35$ |

Each drop's class is drawn by inverse-CDF sampling of these ratios. A malnutrition-like diet is
low-diversity and calorie-dominated; the MDCF arm is enriched for the targeted prebiotic, so taxa
with high `aff_mdcf` gain a selective advantage only under that diet. These ratios are calibration
estimates flagged for refinement against empirical intake data. The active arm is selectable in the
**Globals** tab: clicking the `diet_profile` cell (or pressing Enter on its row) cycles
`malnutrition` $\to$ `standard` $\to$ `mdcf` $\to$ `malnutrition`, applied on the next reset/reload.


### 4.20 Spatial Heterogeneity: Lumen vs Mucus Layer (Phase I)

The circular petri dish is geometrically partitioned into two concentric zones that abstract the
radial organisation of the intestinal cross-section<sup>[6](#ref-6), [18](#ref-18)</sup>: a central, well-mixed
**lumen** and a thin outer **mucus** annulus. For a point at world position $(x, y)$, with petri
centre $(c_x, c_y)$ and radius $R$, the radial fraction is

$$\rho = \frac{\sqrt{(x - c_x)^2 + (y - c_y)^2}}{R}.$$

A point lies in the mucus layer when $\rho \ge$ `MUCUS_LAYER_RADIUS_FRACTION` (default $0.75$), so
the inner $75\%$ of the radius is lumen and the outer $25\%$ is the mucus annulus. Four zone-specific
rules differentiate the niches.

**1. Motility (anchoring refuge).** Mucus-associated cells adhere to the glycoprotein matrix and
move much more slowly than free luminal swimmers<sup>[18](#ref-18)</sup>. Before the forward translation,
move much more slowly than free luminal swimmers<sup>[18](#ref-18)</sup>. Before the forward translation,
a cell whose current position satisfies $\rho \ge$ `MUCUS_LAYER_RADIUS_FRACTION` has its per-step
speed scaled by the mucus motility factor $m$ (constant `MUCUS_MOTILITY_FACTOR` $= 0.3$), giving
effective speed $v_{\text{eff}} = v\,m$. This makes the mucus a low-clearance refuge where cells
persist rather than being swept downstream.

**2. Nutrient washout (peristaltic clearance).** Free food mass decays exponentially per real
second at a zone-specific rate, modelling luminal-flow / peristaltic clearance. For a node with
mass $m_t$ in a zone of washout rate $k$,

$$m_{t+\Delta t} = m_t\,(1 - k\,\Delta t),$$

with $k =$ `LUMEN_FOOD_WASHOUT_RATE_PER_SEC` $(= 0.02)$ in the lumen and $k =$
`MUCUS_FOOD_WASHOUT_RATE_PER_SEC` $(= 0.002)$ in the mucus. The lumen is flushed an order of
magnitude faster, so the mucus stabilises a nutrient reservoir.

**3. Drop deposition.** Fresh dietary substrate is delivered chiefly to the well-mixed lumen by
luminal flow. A drop landing in the mucus annulus ($\rho \ge$ `MUCUS_LAYER_RADIUS_FRACTION`)
retains only a fraction `MUCUS_DROP_AMOUNT_FACTOR` $(= 0.5)$ of its sampled mass, while lumen
drops are deposited in full.

**4. Competition bias.** Under hard carrying-capacity pressure the well-mixed lumen imposes
stronger resource competition than the stabilised mucus refuge, so luminal cells are preferentially
lost. Culling uses weighted sampling without replacement<sup>[7](#ref-7)</sup>: each
lost. Culling uses weighted sampling without replacement<sup>[7](#ref-7)</sup>: each
candidate cell draws a key $U^{1/w}$ with $U \sim \mathrm{Uniform}(0,1)$ and zone weight $w = 1$
(mucus) or $w =$ `LUMEN_COMPETITION_CULL_WEIGHT` $(= 3.0)$ (lumen); the cells with the largest keys
are removed. With equal weights this reduces exactly to uniform random culling, preserving legacy
behaviour, while the up-weighted lumen is thinned first.

Together these rules make the lumen a high-throughput, high-turnover compartment and the mucus a
slow, protected refuge — the spatial substrate for distinct extinction/abundance dynamics studied
in Phase I.

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

$$T_s^{\text{ext}} = \min\{\, t : N_s(t) = 0 \,\},$$

recorded once and never overwritten. Survivors retain $T_s^{\text{ext}} = \varnothing$ (blank in
the export). Elapsed time is the real-time accumulator $t_{\text{elapsed},n} = \sum_{i=0}^{n}\Delta t_i$,
so the metric is frame-rate independent.

**Export.** `simulation_stats.csv` gains an `--- Extinction Events ---` section with one row per
species: `Species`, `Time_To_Extinction_s`, `Extinct` (0/1), and `Final_Count`. The run's `seed`
is also written into the environment-configuration section so a baseline can be exactly reproduced
from its archived CSV.

**Survival curve.** From the per-species extinction times we build a discrete Kaplan-Meier-style
survival function<sup>[16](#ref-16)</sup> $S(t)$ — the fraction of the $M$ tracked species still extant at time $t$:
survival function<sup>[16](#ref-16)</sup> $S(t)$ — the fraction of the $M$ tracked species still extant at time $t$:

$$S(t) = \frac{1}{M}\,\bigl|\{\, s : T_s^{\text{ext}} > t \ \text{or}\ T_s^{\text{ext}} = \varnothing \,\}\bigr|,$$

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
on the current tick, indexed by substrate id $s \in \{0, 1, 2\}$.

**Community accumulation.** Each update the simulation sums the per-cell ledgers into a per-step
total $\Delta M_s$ and carries forward the previous cumulative value, so the recorded series is the
total mass ever ingested per class:

$$M_s(t_n) = M_s(t_{n-1}) + \Delta M_s(t_n), \qquad M_s(t_0) = \Delta M_s(t_0).$$

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

$$r = \sqrt{U\,(r_{\max}^2 - r_{\min}^2) + r_{\min}^2}, \qquad U \sim \mathrm{Uniform}(0,1),$$

so density is constant per unit **area** rather than per unit radius (a naive $r = U\,R$ draw would
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
  $[\,2\pi i/n,\; 2\pi (i{+}1)/n\,)$, separating founder origins to study invasion / range-expansion
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
* `temp_min`: Ratkowsky minimum cardinal temperature $T_{min}$ in degrees Celsius. No growth occurs at or below it (see Section 4.14)<sup>[23](#ref-23)</sup>.
* `temp_min`: Ratkowsky minimum cardinal temperature $T_{min}$ in degrees Celsius. No growth occurs at or below it (see Section 4.14)<sup>[23](#ref-23)</sup>.
* `temp_max`: Ratkowsky maximum cardinal temperature $T_{max}$ in degrees Celsius. No growth occurs at or above it. It is clamped to be strictly above `temp_min`.
* `acetate_secretion`: Overflow acetate excreted per unit of fresh glucose uptake, from 0 to 1. Producers emit acetate into the ECM (see Section 4.15). Default 0 (no overflow)<sup>[30](#ref-30)</sup>.
* `acetate_secretion`: Overflow acetate excreted per unit of fresh glucose uptake, from 0 to 1. Producers emit acetate into the ECM (see Section 4.15). Default 0 (no overflow)<sup>[30](#ref-30)</sup>.
* `acetate_uptake`: Acetate cross-feeding uptake activity from 0 to 1. Consumers oxidize local acetate as a secondary aerobic carbon source (see Section 4.15). Default 0 (no cross-feeding).
* `baranyi_q0`: Initial Baranyi physiological state $q_0$ of the fresh inoculum (dimensionless, $> 0$). A smaller value lengthens the lag phase via $h_0 = \ln(1 + 1/q_0)$ (see Section 4.16); it is clamped to $[10^{-4}, 10^3]$. Default 0.05 (moderate lag)<sup>[2](#ref-2)</sup>.
* `persister_switch`: Persister switching propensity from 0 to 1 (see Section 4.17). It scales the entry rate into the dormant persister state; 0 fully disables persistence (pre-M4 behaviour). Default 0<sup>[1](#ref-1), [20](#ref-20)</sup>.
* `baranyi_q0`: Initial Baranyi physiological state $q_0$ of the fresh inoculum (dimensionless, $> 0$). A smaller value lengthens the lag phase via $h_0 = \ln(1 + 1/q_0)$ (see Section 4.16); it is clamped to $[10^{-4}, 10^3]$. Default 0.05 (moderate lag)<sup>[2](#ref-2)</sup>.
* `persister_switch`: Persister switching propensity from 0 to 1 (see Section 4.17). It scales the entry rate into the dormant persister state; 0 fully disables persistence (pre-M4 behaviour). Default 0<sup>[1](#ref-1), [20](#ref-20)</sup>.
* `persister_enter_thresh`: Starvation fullness threshold $S_1$ from 0 to 1. Below it, starvation-induced switching into persistence engages. Default 0.15.
* `persister_exit_thresh`: Recovery fullness threshold $S_2$ from 0 to 1. Above it, persisters resuscitate to the proliferative state. It is clamped to stay at least `PERSISTER_THRESHOLD_MIN_GAP` above $S_1$. Default 0.5.
* `phage_resistance`: Initial fraction of cells born resistant to bacteriophage adsorption, from 0 to 1 (see Section 4.18). Resistance is heritable; resistant cells survive an epidemic and drive the rebound. Default 0 (fully susceptible)<sup>[21](#ref-21)</sup>.
* `hgt_competence`: Horizontal gene transfer competence from 0 to 1 (see Section 4.18). Gates both donating and receiving phage resistance on cell contact; the recipient inherits a metabolic carriage cost. 0 disables HGT (pre-M5 behaviour). Default 0<sup>[22](#ref-22), [26](#ref-26)</sup>.
* `phage_resistance`: Initial fraction of cells born resistant to bacteriophage adsorption, from 0 to 1 (see Section 4.18). Resistance is heritable; resistant cells survive an epidemic and drive the rebound. Default 0 (fully susceptible)<sup>[21](#ref-21)</sup>.
* `hgt_competence`: Horizontal gene transfer competence from 0 to 1 (see Section 4.18). Gates both donating and receiving phage resistance on cell contact; the recipient inherits a metabolic carriage cost. 0 disables HGT (pre-M5 behaviour). Default 0<sup>[22](#ref-22), [26](#ref-26)</sup>.
* `temp_min`: Ratkowsky minimum cardinal temperature $T_{min}$ in degrees Celsius. No growth occurs at or below it (see Section 4.14).
* `temp_max`: Ratkowsky maximum cardinal temperature $T_{max}$ in degrees Celsius. No growth occurs at or above it. It is clamped to be strictly above `temp_min`.
* `acetate_secretion`: Overflow acetate excreted per unit of fresh glucose uptake, from 0 to 1. Producers emit acetate into the ECM (see Section 4.15). Default 0 (no overflow).
* `acetate_uptake`: Acetate cross-feeding uptake activity from 0 to 1. Consumers oxidize local acetate as a secondary aerobic carbon source (see Section 4.15). Default 0 (no cross-feeding).
* `baranyi_q0`: Initial Baranyi physiological state $q_0$ of the fresh inoculum (dimensionless, $> 0$). A smaller value lengthens the lag phase via $h_0 = \ln(1 + 1/q_0)$ (see Section 4.16); it is clamped to $[10^{-4}, 10^3]$. Default 0.05 (moderate lag).
* `persister_switch`: Persister switching propensity from 0 to 1 (see Section 4.17). It scales the entry rate into the dormant persister state; 0 fully disables persistence (pre-M4 behaviour). Default 0.
* `persister_enter_thresh`: Starvation fullness threshold $S_1$ from 0 to 1. Below it, starvation-induced switching into persistence engages. Default 0.15.
* `persister_exit_thresh`: Recovery fullness threshold $S_2$ from 0 to 1. Above it, persisters resuscitate to the proliferative state. It is clamped to stay at least `PERSISTER_THRESHOLD_MIN_GAP` above $S_1$. Default 0.5.
* `phage_resistance`: Initial fraction of cells born resistant to bacteriophage adsorption, from 0 to 1 (see Section 4.18). Resistance is heritable; resistant cells survive an epidemic and drive the rebound. Default 0 (fully susceptible).
* `hgt_competence`: Horizontal gene transfer competence from 0 to 1 (see Section 4.18). Gates both donating and receiving phage resistance on cell contact; the recipient inherits a metabolic carriage cost. 0 disables HGT (pre-M5 behaviour). Default 0.
* `aff_calorie`: Generic-calorie substrate uptake affinity ($\ge 0$, see Section 4.19). Multiplies transporter throughput for simple sugars; energy yield is $1.0$. Default 1.0 (full legacy uptake)<sup>[10](#ref-10), [11](#ref-11), [24](#ref-24)</sup>.
* `aff_calorie`: Generic-calorie substrate uptake affinity ($\ge 0$, see Section 4.19). Multiplies transporter throughput for simple sugars; energy yield is $1.0$. Default 1.0 (full legacy uptake)<sup>[10](#ref-10), [11](#ref-11), [24](#ref-24)</sup>.
* `aff_glycan`: Complex-glycan (dietary fibre) uptake affinity ($\ge 0$, see Section 4.19). 0 means the species cannot metabolise fibre; ingested mass yields $0.7$ of glucose energy. Default 0.
* `aff_mdcf`: MDCF prebiotic uptake affinity ($\ge 0$, see Section 4.19). 0 means the species cannot use the targeted prebiotic; ingested mass yields $0.85$ of glucose energy. High values are selectively rewarded under the `mdcf` diet. Default 0.

## 6. Global Environment Parameters

Globals are stored in `env_config.json` and can be edited in the Globals tab.

* `o2`: Dissolved oxygen fraction from 0 to 1. Aerobes benefit from oxygen; anaerobes are harmed by it.
* `co2`: Dissolved carbon dioxide fraction from 0 to 1. Respiration increases it and buffering/pH logic responds to it.
* `ph`: Current environmental pH from 0 to 14.
* `temperature`: Environment temperature in degrees Celsius (slider range 0 to 60). It drives the Ratkowsky growth multiplier $f_{T}$ applied per species (see Section 4.14)<sup>[23](#ref-23)</sup>.
* `temperature`: Environment temperature in degrees Celsius (slider range 0 to 60). It drives the Ratkowsky growth multiplier $f_{T}$ applied per species (see Section 4.14)<sup>[23](#ref-23)</sup>.
* `temperature`: Environment temperature in degrees Celsius (slider range 0 to 60). It drives the Ratkowsky growth multiplier $f_{T}$ applied per species (see Section 4.14).
* `capacity`: Hard carrying capacity across all species. Excess population is culled if total cells exceed this value.
* `drop_frequency`: Expected food clusters per real second. The implementation uses an expected-rate spawning approximation $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$.
* `drop_range_min` and `drop_range_max`: Uniform random range for the amount in each new food cluster.
* `phage_inflow`: Expected free-phage seed packets introduced per real second (see Section 4.18). Each packet deposits `PHAGE_SEED_PACKET_VIRIONS` virions at a random petri location via the same Poisson rate $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$. Set above 0 to trigger an epidemic and back to 0 to halt seeding. Default 0 (no outbreak)<sup>[8](#ref-8), [19](#ref-19)</sup>.
* `phage_inflow`: Expected free-phage seed packets introduced per real second (see Section 4.18). Each packet deposits `PHAGE_SEED_PACKET_VIRIONS` virions at a random petri location via the same Poisson rate $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$. Set above 0 to trigger an epidemic and back to 0 to halt seeding. Default 0 (no outbreak)<sup>[8](#ref-8), [19](#ref-19)</sup>.
* `phage_inflow`: Expected free-phage seed packets introduced per real second (see Section 4.18). Each packet deposits `PHAGE_SEED_PACKET_VIRIONS` virions at a random petri location via the same Poisson rate $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$. Set above 0 to trigger an epidemic and back to 0 to halt seeding. Default 0 (no outbreak).
* `diet_profile`: Active dietary feeding arm (see Section 4.19), one of `malnutrition`, `standard`, or `mdcf`. It sets the categorical probability mix of substrate classes (generic calorie / complex glycan / MDCF prebiotic) assigned to each new food drop<sup>[10](#ref-10), [11](#ref-11), [24](#ref-24)</sup>. It is editable directly in the **Globals** tab — clicking its cell (or pressing Enter on the row) cycles to the next arm. Like all Globals-tab values it defines starting conditions, so a change takes effect on the next simulation reset/reload. Default `standard`.
* `seed`: Optional integer RNG seed for reproducible extinction/abundance baselines (see Section 4.21)<sup>[16](#ref-16)</sup>. It is editable directly in the **Globals** tab — enter a non-negative integer to fix the run, or leave the field blank (shown as `(random)`) for a fresh nondeterministic run. When set, the random stream is seeded before species construction so the entire run — initial layout and every stochastic draw — is deterministic. Like all Globals-tab values it defines starting conditions, so a change takes effect on the next simulation reset/reload.
* `diet_profile`: Active dietary feeding arm (see Section 4.19), one of `malnutrition`, `standard`, or `mdcf`. It sets the categorical probability mix of substrate classes (generic calorie / complex glycan / MDCF prebiotic) assigned to each new food drop<sup>[10](#ref-10), [11](#ref-11), [24](#ref-24)</sup>. It is editable directly in the **Globals** tab — clicking its cell (or pressing Enter on the row) cycles to the next arm. Like all Globals-tab values it defines starting conditions, so a change takes effect on the next simulation reset/reload. Default `standard`.
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
time,total,births,deaths,elapsed_seconds,dt_seconds,tickrate_hz,overall_births,overall_deaths,lumen,mucus,consumed_calorie,consumed_glycan,consumed_mdcf
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
* `lumen`: live bacteria in the central lumen (radial fraction below the mucus threshold) after that update (see Section 4.20).
* `mucus`: live bacteria in the outer mucus annulus (radial fraction at or above the mucus threshold) after that update. `lumen` and `mucus` always sum to `total`. These columns are blank for stats produced by legacy imports that predate zone tracking.
* `consumed_calorie`: cumulative generic-calorie substrate mass (pg) ingested by the whole community through this row (see Section 4.22).
* `consumed_glycan`: cumulative complex-glycan substrate mass (pg) ingested through this row.
* `consumed_mdcf`: cumulative MDCF-prebiotic substrate mass (pg) ingested through this row. Each `consumed_*` series is monotonic non-decreasing; the three columns are blank for stats produced by legacy imports that predate diet analytics.

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
\Delta t = \frac{1}{\text{FPS}}, \quad t_{\text{elapsed}} = \frac{\text{time}+1}{\text{FPS}}, \quad f_{\text{tick}}=\text{FPS}
$$

and reconstructs cumulative global births and deaths by running sum.

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

* `bacteria/app.py`: Pygame application loop, rendering, tabs, UI state, and high-level simulation orchestration.
* `bacteria/core_sim.py`: Headless biological and physical simulation models: `FoodNode`, `Bacterium`, `BacteriaSpecies`, `FoodSpatialHash`, and `Simulation`.
* `bacteria/ecm.py`: Extracellular matrix diffusion and decay for quorum-sensing molecules, bacteriocins, and cross-feeding acetate.
* `bacteria/ecm.py`: Extracellular matrix diffusion and decay for quorum-sensing molecules, bacteriocins, and cross-feeding acetate.
* `bacteria/data.py`: CSV/JSON serialization for `species_data.csv` and `env_config.json`.
* `bacteria/ui_components.py`: Reusable buttons, sliders, and editable cells.
* `bacteria/constants.py`: UI constants plus documented biological and numerical model constants.
* `bacteria/__main__.py`: CLI entry point for `python -m bacteria`.

## 10. Performance Notes

* Food nodes and bacteria use small object pools to reduce allocation churn.
* Food lookup uses `FoodSpatialHash` so bacteria query nearby food instead of scanning the full node list.
* Bacterial physical interactions use a spatial hash over neighboring cells.
* ECM diffusion updates every active cell each tick, but only stores non-negligible concentrations.
* Viewport culling suppresses expensive local chemotaxis for off-viewport bacteria while retaining simple movement.

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

<a id="ref-12"></a>**[12]** Gompertz (1825). *On the nature of the function expressive of the law of human mortality*. Philosophical Transactions of the Royal Society of London, 115:513–583. [https://doi.org/10.1098/rstl.1825.0026](https://doi.org/10.1098/rstl.1825.0026)

<a id="ref-13"></a>**[13]** Harman (1956). *Aging: a theory based on free radical and radiation chemistry*. Journal of Gerontology, 11(3):298–300. [https://doi.org/10.1093/geronj/11.3.298](https://doi.org/10.1093/geronj/11.3.298)

<a id="ref-14"></a>**[14]** Henderson (1908). *Concerning the relationship between the strength of acids and their capacity to preserve neutrality*. American Journal of Physiology, 21(2):173–179. [https://doi.org/10.1152/ajplegacy.1908.21.2.173](https://doi.org/10.1152/ajplegacy.1908.21.2.173)

<a id="ref-15"></a>**[15]** Kalinin, Jiang, Tu & Wu (2009). *Logarithmic sensing in Escherichia coli bacterial chemotaxis*. Biophysical Journal, 96(6):2439–2448. [https://doi.org/10.1016/j.bpj.2008.10.027](https://doi.org/10.1016/j.bpj.2008.10.027)

<a id="ref-16"></a>**[16]** Kaplan & Meier (1958). *Nonparametric estimation from incomplete observations*. Journal of the American Statistical Association, 53(282):457–481. [https://doi.org/10.1080/01621459.1958.10501452](https://doi.org/10.1080/01621459.1958.10501452)

<a id="ref-17"></a>**[17]** Kleiber (1932). *Body size and metabolism*. Hilgardia, 6(11):315–353. [https://doi.org/10.3733/hilg.v06n11p315](https://doi.org/10.3733/hilg.v06n11p315)

<a id="ref-18"></a>**[18]** Earle, Billings, Sigal, Lichtman, Hansson, Elias, Amieva, Huang & Sonnenburg (2015). *Quantitative Imaging of Gut Microbiota Spatial Organization*. Cell Host & Microbe, 18(4):478–488. [https://doi.org/10.1016/j.chom.2015.09.002](https://doi.org/10.1016/j.chom.2015.09.002) — quantitatively maps the spatial distribution of the gut microbiota relative to the host mucus and epithelium, providing the empirical basis for the §4.20 dual-compartment radial (lumen vs mucus) model and its mucus-associated anchoring refuge.

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

<a id="ref-35"></a>**[35]** **Own (BacSandbox model).** Bacteriocin and T6SS damage/immunity formulations and constants, §4.12.

<a id="ref-36"></a>**[36]** **Own (BacSandbox model).** Age-dependent maintenance factor A(t)=1+t/(180·L), §4.11.

<a id="ref-37"></a>**[37]** **Own (BacSandbox model).** Modeled dissolved-gas saturation capacities and exponential gas/pH relaxation rates, §4.10.

<a id="ref-38"></a>**[38]** **Own (BacSandbox model).** Initial-inoculum spawning policy (lumen-bias, clustered microcolonies) and elastic boundary reflection, §4.6/§4.23.