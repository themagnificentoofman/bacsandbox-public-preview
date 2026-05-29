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

The extracellular matrix tracks quorum-sensing molecules, bacteriocins, and cross-feeding acetate (Section 4.15) on a spatial hash grid. ECM diffusion now processes every active grid cell each simulation tick. It uses an explicit Fickian 2D stencil with a stability clamp:
The extracellular matrix tracks quorum-sensing molecules, bacteriocins, and cross-feeding acetate (Section 4.15) on a spatial hash grid. ECM diffusion now processes every active grid cell each simulation tick. It uses an explicit Fickian 2D stencil with a stability clamp:

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

Temperature kinetics (Ratkowsky model):

* `RATKOWSKY_B_SQRT_PER_C = 0.0410`
* `RATKOWSKY_C_PER_C = 0.161`
* `DEFAULT_TEMP_MIN_C = 4.0`
* `DEFAULT_TEMP_MAX_C = 47.0`
* `DEFAULT_ENV_TEMPERATURE_C = 37.0`
* `ENV_TEMPERATURE_MIN_C = 0.0`
* `ENV_TEMPERATURE_MAX_C = 60.0`
* `RATKOWSKY_NORMALIZATION_SAMPLES = 200`

Cross-feeding (acetate overflow, see Section 4.15):

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

Lag phase (Baranyi & Roberts dynamic state, see Section 4.16):

* `DEFAULT_BARANYI_Q0 = 0.05` (fresh-inoculum physiological state; $h_0 = \ln(1 + 1/q_0) \approx 3.04$)
* `BARANYI_Q_STATE_MAX` $= 1.0 \times 10^{6}$ (numerical cap where $\alpha > 0.999999$)
* `BARANYI_Q0_MIN` $= 1.0 \times 10^{-4}$ (floor; guarantees a finite lag)
* `BARANYI_Q0_MAX` $= 1.0 \times 10^{3}$ (ceiling; effectively disables the lag)

Phenotypic switching (persister bet-hedging, see Section 4.17):

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

Phage predation and horizontal gene transfer (see Section 4.18):

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
derived from the four-parameter Ratkowsky square-root model (Ratkowsky, Lowry, McMeekin, Stokes
& Chandler, 1983), which is biologically validated across the entire biokinetic temperature range:

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
switch** (Wolfe, 2005, *Microbiol. Mol. Biol. Rev.* 69:12–50): during fast aerobic growth on
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

so $\Delta F = A_{\text{con}} \, Y_E$. The ATP yields ($\sim 32$ ATP/glucose, Rich 2003; $\sim 10$
ATP/acetate via acetyl-CoA → TCA cycle + electron transport chain) make acetate a slightly less
energy-dense carbon source per unit mass, as observed experimentally.

**Gas coupling.** Acetate oxidation follows $\mathrm{C_2H_4O_2 + 2\,O_2 \rightarrow 2\,CO_2 + 2\,H_2O}$.
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
the dynamic approach of Baranyi & Roberts (1994), which augments each cell with a dimensionless
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
antibiotic *tolerance* (as distinct from inherited resistance). BacSandbox models it as a two-state
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

References: Balaban, Merrin, Chait, Kowalik & Leibler (2004), "Bacterial Persistence as a Phenotypic
Switch", *Science* 305:1622–1625; Lewis (2010), "Persister Cells", *Annu. Rev. Microbiol.* 64:357–372.

### 4.18 Phage Predation and Horizontal Gene Transfer (SIR Lytic Cycle)

Lytic bacteriophages are the most abundant biological entities on Earth and the dominant top-down
control on bacterial populations. BacSandbox models them not as discrete agents but as a diffusible
**free-virion field** carried in the ECM (`phage` layer), coupled to a per-cell susceptible /
infected / resistant (SIR) state machine. This reuses the existing reaction–diffusion infrastructure
and keeps the model lightweight while reproducing the classic predator–prey boom–bust oscillation.

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
structure (a fixed latent period followed by a discrete burst) follows Ellis & Delbrück (1939).

**Free-virion field.** Released virions diffuse (`ECM_PHAGE_DIFFUSION_UM2_S`) and lose infectivity
(`ECM_PHAGE_DECAY_S`) via the same Fickian stencil used for the other ECM layers (Section 4.9),
spreading the epidemic spatially and decaying it where hosts are absent.

**Resistance, mutation, and rebound.** Each cell is born resistant with probability equal to its
species' `phage_resistance` trait (a heritable, chromosomal fraction). Resistant cells never adsorb
phage, so they survive the crash. On every division a susceptible lineage can spontaneously mutate to
resistant with probability `PHAGE_RESISTANCE_MUTATION_RATE` (a Luria–Delbrück-style fluctuation),
seeding the post-crash recovery. Daughters inherit the parent's resistance state and carriage cost.
The net result is the textbook sequence: phage outbreak $\rightarrow$ susceptible population crash
$\rightarrow$ resistant-only survivors $\rightarrow$ rebound.

**Horizontal gene transfer (HGT).** Resistance also spreads laterally. When a resistant donor and a
susceptible recipient are in physical contact (the same contact test used for sharing/T6SS), and both
species are HGT-competent (`hgt_competence` $> 0$), the trait transfers by a conjugation-like process
with probability

$$
P_{\text{HGT}} = 1 - e^{-\,k_H\, c_{\text{donor}}\, c_{\text{recip}}\,\Delta t},
$$

where $k_H$ is `HGT_TRANSFER_RATE_PER_SEC` and $c$ are the two cells' competences. The recipient
becomes resistant and inherits a **plasmid carriage cost**: its basal maintenance metabolism is
multiplied by $1 +$ `HGT_METABOLIC_BURDEN`, modelling the well-documented fitness penalty of
maintaining a resistance plasmid (San Millán & MacLean 2017). This cost is inherited by daughters, so
in the absence of phage pressure the burdened lineage is slightly outcompeted — capturing the
cost–benefit tension that governs the spread of mobile resistance elements.

**Outbreak control.** The global `phage_inflow` environment parameter seeds free virions via a
Poisson process ($\lambda\,\Delta t$ packets of `PHAGE_SEED_PACKET_VIRIONS` at random petri
locations), so a user can trigger or halt an epidemic live from the Globals tab. A value of $0$
(the default) means no outbreak and exactly reproduces the pre-M5 behaviour.

References: Ellis & Delbrück (1939), "The Growth of Bacteriophage", *J. Gen. Physiol.* 22:365–384;
Levin, Stewart & Chao (1977), "Resource-Limited Growth, Competition, and Predation: A Model and
Experimental Studies with Bacteria and Bacteriophage", *Am. Nat.* 111:3–24; Ochman, Lawrence &
Groisman (2000), "Lateral gene transfer and the nature of bacterial innovation", *Nature*
405:299–304; San Millán & MacLean (2017), "Fitness Costs of Plasmids", *Microbiol. Spectr.* 5(5).


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

## 6. Global Environment Parameters

Globals are stored in `env_config.json` and can be edited in the Globals tab.

* `o2`: Dissolved oxygen fraction from 0 to 1. Aerobes benefit from oxygen; anaerobes are harmed by it.
* `co2`: Dissolved carbon dioxide fraction from 0 to 1. Respiration increases it and buffering/pH logic responds to it.
* `ph`: Current environmental pH from 0 to 14.
* `temperature`: Environment temperature in degrees Celsius (slider range 0 to 60). It drives the Ratkowsky growth multiplier $f_{T}$ applied per species (see Section 4.14).
* `temperature`: Environment temperature in degrees Celsius (slider range 0 to 60). It drives the Ratkowsky growth multiplier $f_{T}$ applied per species (see Section 4.14).
* `capacity`: Hard carrying capacity across all species. Excess population is culled if total cells exceed this value.
* `drop_frequency`: Expected food clusters per real second. The implementation uses an expected-rate spawning approximation $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$.
* `drop_range_min` and `drop_range_max`: Uniform random range for the amount in each new food cluster.
* `phage_inflow`: Expected free-phage seed packets introduced per real second (see Section 4.18). Each packet deposits `PHAGE_SEED_PACKET_VIRIONS` virions at a random petri location via the same Poisson rate $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$. Set above 0 to trigger an epidemic and back to 0 to halt seeding. Default 0 (no outbreak).
* `phage_inflow`: Expected free-phage seed packets introduced per real second (see Section 4.18). Each packet deposits `PHAGE_SEED_PACKET_VIRIONS` virions at a random petri location via the same Poisson rate $\mathbb{E}[N_{\Delta t}] = \lambda \Delta t$. Set above 0 to trigger an epidemic and back to 0 to halt seeding. Default 0 (no outbreak).

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