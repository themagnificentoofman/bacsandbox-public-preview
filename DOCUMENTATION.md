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

Nutrient capacity scales with structural volume, not temporary fullness:

$$
C_{\text{food}} = k_{\text{cap}} \cdot s_{\text{base}}^3
$$

The rendered cell radius uses fullness only as a small starvation cue:

$$
s_{\text{rendered}} = s_{\text{base}} \left(0.9 + 0.1 \cdot \frac{F}{C_{\text{food}}}\right)
$$

This avoids the unrealistic behavior where eating more food directly increases the maximum storage cap.

### 4.2 Metabolism And Uptake

Basal metabolism follows an allometric Kleiber-style relationship. Because cell mass scales with volume, a mass exponent of 3/4 becomes a size exponent of 2.25:

$$
B = \beta \cdot s_{\text{base}}^{2.25}
$$

Nutrient uptake is surface-area limited:

$$
U_{\max} = u \cdot s_{\text{base}}^2
$$

During each elapsed-time update:

$$
F_{t+\Delta t} = F_t - B \cdot \Delta t \cdot A(t)
$$

where `A(t)` is the senescence maintenance factor.

### 4.3 Food Sensing And Chemotaxis

Bacteria do not teleport or directly seek a target. They sense local nutrient fields and bias their swimming direction and run duration. Each food node contributes a softened Euclidean concentration proxy:

$$
C_i = \frac{A_i}{\pi(d_i^2 + a^2)}
$$

where $A_i$ is node amount, $d_i$ is Euclidean distance, and $a = 5.0\,\mu\text{m}$ prevents a zero-distance singularity.

Chemoreceptor response is log-compressed to approximate saturation:

$$
S_i = g_{\text{eff}} \log\left(1 + \frac{C_i}{C_{\text{ref}}}\right)
$$

where $C_{\text{ref}} = 10^{-3}\,\text{pg}/\mu\text{m}^2$ and $g_{\text{eff}}$ increases as the bacterium becomes hungry.

The local gradient direction is a weighted sum of node direction vectors:

$$
\vec{G} = \sum_i S_i \frac{\vec{x}_i - \vec{x}_{\text{cell}}}{\lVert \vec{x}_i - \vec{x}_{\text{cell}} \rVert}
$$

### 4.4 Run-And-Tumble Motion

BacSandbox uses an E. coli-inspired run-and-tumble model. Runs are sampled around 1 second and clamped to biological bounds:

$$
t_{\text{run}} \sim \operatorname{clamp}(\mathcal{N}(1.0, 0.2^2), 0.25, 3.0)
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
\Delta \theta_{\text{run}} = \operatorname{clamp}(\theta_G - \theta, -\omega g_{\text{eff}} \Delta t, \omega g_{\text{eff}} \Delta t)
$$

During a tumble, the cell performs one event-scale reorientation rather than a tiny frame-scaled wobble. The unbiased magnitude is centered near the classic observed E. coli reorientation angle of 68 degrees:

$$
|\Delta \theta_{\text{tumble}}| \sim \mathcal{N}(1.19, 0.45^2)\,\text{rad}
$$

Translation during tumbles is multiplied by 0.05, representing motor activity being spent on reorientation rather than propulsion. Speed is never negative; apparent reversal is a new heading after a tumble, not backward swimming.

### 4.5 Rotational Diffusion

Brownian rotational noise follows the low-Reynolds scaling for small particles in fluid:

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

The petri dish is circular. If a bacterium reaches the dish boundary, it is moved slightly inside the rim and its velocity heading reflects off the outward normal only when it was moving outward:

$$
\vec{v}' = \vec{v} - 2(\vec{v} \cdot \hat{n})\hat{n}
$$

This is applied only when:

$$
\vec{v} \cdot \hat{n} > 0
$$

### 4.7 Reproduction

Population-level reproduction probability follows Verhulst logistic pressure through carrying capacity:

$$
\frac{dN}{dt} = rN\left(1 - \frac{N}{K}\right)
$$

At the individual level, binary fission requires stored nutrient above the reproduction threshold and a synthesis cost inspired by Stouthamer-style ATP/biomass costs:

$$
F_{\text{cost}} = 40.0 \cdot s_{\text{base}}^3
$$

The parent divides only if both parent and child retain enough nutrient to survive the next update.

The implementation exposes this conversion as `FISSION_COST_PG_PER_UM3 = 40.0` so tests and model code share the same biological assumption.

### 4.8 Death, Predation, And Necromass Recycling

Dead bacteria do not become pure glucose. Death releases a lossy glucose-equivalent nutrient value from cytosolic leakage plus labile dry necromass:

$$
F_{\text{recovered}} = \max(0, F_{\text{stored}}) f_{\text{cytosol}} + V_{\text{cell}} m_{\text{wet}} f_{\text{dry}} f_{\text{labile}}
$$

The implemented constants are:

$$
f_{\text{cytosol}} = 0.75, \quad m_{\text{wet}} = 1.0\,\text{pg}/\mu\text{m}^3, \quad f_{\text{dry}} = 0.30, \quad f_{\text{labile}} = 0.35
$$

Most structural biomass is therefore unavailable, refractory, respired, or lost to the environment. Predation uses the same digestible nutrient logic instead of consuming a fixed fraction of the prey's storage capacity.

### 4.9 ECM Diffusion Grid

The extracellular matrix tracks quorum-sensing molecules and bacteriocins on a spatial hash grid. ECM diffusion now processes every active grid cell each simulation tick. It uses an explicit Fickian 2D stencil with a stability clamp:

$$
\alpha = \min\left(\alpha_{\max}, \frac{D\Delta t}{\Delta x^2}\right), \quad \alpha_{\max}=0.24
$$

Decay is applied as:

$$
m = \operatorname{clamp}(1 - \lambda\Delta t, 0, 1)
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

For aerobic species with positive oxygen sensitivity $s_{O_2}$, respiration follows:

$$
\text{C}_6\text{H}_{12}\text{O}_6 + 6\text{O}_2 \rightarrow 6\text{CO}_2 + 6\text{H}_2\text{O}
$$

The implementation converts glucose-equivalent pg into dissolved gas fractions with molar mass ratios and model capacities:

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

For anaerobic species with negative oxygen sensitivity, the simulator uses a simplified lactic-acid pathway plus trace heterolactic carbon dioxide:

$$
\Delta CO_{2,\text{anaerobic}} = \frac{0.1G_{\Delta t}|s_{O_2}|}{1.5\times10^9}
$$

$$
\Delta \text{acid}_{\text{lactic}} = \frac{2.5G_{\Delta t}|s_{O_2}|}{4.0\times10^7}
$$

The update sequence is: subtract oxygen use, add carbon dioxide production, subtract acid from pH, equilibrate gases toward an open dish, then buffer pH toward neutral. Gas and buffer restoration are:

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

Maintenance cost rises with age:

$$
A(t) = 1 + \frac{t}{180 \cdot L}
$$

Intrinsic death risk follows a Gompertz-style exponential penalty:

$$
\mu(t) = r_{\text{aging}} \Delta t \cdot e^{\alpha t}
$$

### 4.12 Competition, Bacteriocins, And T6SS

Diffusible bacteriocin production consumes stored nutrient and emits antimicrobial material into the ECM:

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

Type VI secretion system contact damage scales with attacker surface area:

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

Movement, rendering, and senescence:

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

## 6. Global Environment Parameters

Globals are stored in `env_config.json` and can be edited in the Globals tab.

* `o2`: Dissolved oxygen fraction from 0 to 1. Aerobes benefit from oxygen; anaerobes are harmed by it.
* `co2`: Dissolved carbon dioxide fraction from 0 to 1. Respiration increases it and buffering/pH logic responds to it.
* `ph`: Current environmental pH from 0 to 14.
* `capacity`: Hard carrying capacity across all species. Excess population is culled if total cells exceed this value.
* `drop_frequency`: Expected food clusters per real second. The implementation uses an expected-rate spawning approximation:

  $$
  \mathbb{E}[N_{\Delta t}] = \lambda \Delta t
  $$

* `drop_range_min` and `drop_range_max`: Uniform random range for the amount in each new food cluster.

## 7. Data Export And Import

`simulation_stats.csv` stores three logical sections when exported from the app:

* Environment configuration.
* Species configuration.
* Global and per-species time-series statistics.

The Global Stats header is:

```csv
time,total,births,deaths,elapsed_seconds,dt_seconds,tickrate_hz,overall_births,overall_deaths
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
* `bacteria/ecm.py`: Extracellular matrix diffusion and decay for quorum-sensing molecules and bacteriocins.
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