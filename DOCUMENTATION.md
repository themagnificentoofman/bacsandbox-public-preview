# BacSandbox: Formal Documentation

## 1. Introduction
BacSandbox is a micro-scale bacterial environment simulator. It simulates a 2D petri-dish-like environment where independent bacterial species interact, consume resources, reproduce, and undergo modeled cellular senescence and death.

## 2. Mathematical and Biological Models
The simulation rigorously applies real-world biological concepts and physics models:

* **Volumetric and Mass Scaling**: Food capacity scales volumetrically $size^3$ based on simple geometric assumptions ($\sim 1 \text{ pg}$ mass per $1 \mu\text{m}^3$ volume).
* **Metabolic Rate (Kleiber's Law)**: Metabolic rate scales as $size^{2.25}$ (or mass$^{3/4}$), dictating the basal glucose consumption rate per tick.
* **Nutrient Absorption (Square-Cube Law)**: Consumption speed scales by surface area ($size^2$), limiting the exchange rate of nutrients.
* **Kinematics (Stokes' Law)**: Locomotion speed scales inversely with radius ($v \propto 1/r$) simulating viscous drag in low Reynolds number environments.
* **Run-and-Tumble Chemotaxis**: E. coli model of navigation. Cells run on average for 1s and tumble for 0.1s to reorient. Chemotactic receptor response scales logarithmically with nutrient concentration falling off via an inverse-square relation to distance.
* **Rotational Diffusion**: Brownian rotational noise is applied, scaling as $\sim 1/r^3$ to simulate thermal fluctuations at the microscale.
* **Reproduction (Modified Verhulst/Logistic Growth)**: Logistic growth governed by environmental carrying capacity and individual cellular fitness (incorporating pH and gas toxicity). Binary fission requires threshold biomass.
* **Cellular Senescence (Free Radical Theory & Gompertz Law)**: 
*   Age induces a linear increase in basal metabolic maintenance.
  * Extrinsic and intrinsic mortality risk increases over time following the Gompertz law of mortality exponentially $\mu(x) = A \cdot e^{\alpha x} + B$.
* **Gas Diffusion & pH Buffering (Fick's First Law / Le Chatelier's)**: Dissolved O2/CO2 equilibrate with atmospheric rates ($\sim 0.288 /\text{s}$). pH continuously drifts toward neutral 7.0 via a simulated bicarbonate buffer system ($\sim 0.144 /\text{s}$).

## 3. Species Characteristics
Species behaviors are determined by diverse attributes configured in `species_data.csv`:
* **Basic Traits**: Name, Initial Count, Reproduction Rate, Aging Rate.
* **Physical Traits**: Base Size ($\mu\text{m}$), Base Speed ($\mu\text{m/s}$), Sensor Radius ($\mu\text{m}$).
* **Sensitivities**: O2 Sensitivity ($-1$ to $1$; $<0$ anaerobes, $>0$ aerobes), Preferred pH ($0-14$), Food Sensitivity.
* **Metabolic Traits**: Food Capacity (pg), Metabolic Rate (pg/tick), Consumption Speed (pg/tick).
* **Behavioral Traits**: Greediness (preference for rich food), Aggression (toxin power), Defense (toxin resistance).
* **Thresholds**: Reproduction threshold, Share threshold, Predation threshold.

## 4. UI and Controls
The application features a real-time UI driven by PyGame, separated into five active tabs, augmented by live sliders:
* **Grid Tab**: The primary simulation view. Includes dynamic panning, zooming (using mouse or on-screen buttons), and a minimap.
* **Graph Tab**: Real-time line charts tracking historical total populations and per-species populations. Includes export and import commands to save or replay timeline stats dynamically.
* **Stats Tab**: Displays numerical readouts of live simulation events tracking active populations, global births, and total deaths globally and per-species.
* **Species Tab**: A deeply interactive data-grid for creating, modifying, duplicating, or deleting species parameters on the fly.
* **Globals Tab**: An embedded editor allowing tuning of starting environmental factors (O2, CO2, pH, Capacity) and food node drop frequencies.
* **Environment Sliders**: Live UI toolbars visible in the Grid, Graph, and Stats views to manually adjust the environment state during simulation.

## 5. Software Architecture
The application employs a strict Separation of Concerns pattern, modularized into the `bacteria` Python package:
* `app.py`: The entry-point class handling the PyGame rendering loop, state synchronization, and Tab dispatching.
* `core_sim.py`: Implements the headless mathematical and physical models (`Bacterium`, `BacteriaSpecies`, `Simulation`, `FoodNode`).
* `data.py`: Manages state serialization, loading, and saving mappings for `species_data.csv` and `env_config.json`.
* `ui_components.py`: Houses decoupled, reusable PyGame widgets (`Button`, `Slider`, `CellEditor`).
* `constants.py`: Defines immutable graphical primitives and simulation timing defaults.
* `__main__.py`: Allows the module to be cleanly booted via the `python -m bacteria` CLI command.

## 6. Advanced Biological & Physiological Engines

### 6.1. Binary Fission Synthesis (Stouthamer 1973)
Cellular division represents a massive metabolic investment rather than simple splitting. Based on Stouthamer's (1973) theoretical models on the ATP required for the synthesis of microbial cell material (chromosomes, ribosomes, and cell wall assembly), binary fission strictly requires stored biomass exceeding:
$$ \text{ATP}_{\text{cost}} = 40.0 \times \text{Volume}_{\text{child}} $$
where the cell volume is calculated as $size^3$. Parent cells must ensure that their glucose/biomass pool can satisfy this synthesis cost while leaving both the parent and child with survivable residual energy.

### 6.2. ECM Diffusion Grid (Fick's Second Law)
The Extracellular Matrix (ECM) employs a discrete spatial hash matrix to model the spread and dilution of biochemical agents such as Quorum Sensing autoinducers and Bacteriocin peptides.
The continuous spatial-temporal dynamics are evaluated via Fick's Second Law of Diffusion coupled with biological protease degradation:
$$ \frac{\partial C}{\partial t} = D \nabla^2 C - \lambda C $$
where $C$ is the molecule concentration, $D$ represents the diffusion scalar (modeling kinematic spread), and $\lambda$ represents the decay rate or environmental half-life.

### 6.3. Secretion Systems and Contact Warfare
The simulation strictly models both localized and contact-dependent competition mechanics:
* **Bacteriocin Production**: Cells excrete diffusible antimicrobial peptides into the ECM by expending metabolic ATP. Competing cells absorb these localized toxins, leading to continuous damage unless mitigated by matching immunity factors. Damage reduction scales as:
  $$ \text{Mitigation}_{\text{Bac}} = 1.0 - \min(0.99, I_{\text{factor}} \times 1.5) $$
* **Type VI Secretion System (T6SS)**: For immediate interspecies combat, cells utilize the T6SS—a spring-loaded needle delivering effector toxins directly into neighbors upon physical collision. The offensive damage injection scales with the attacker's physical surface area ($size^2$) and inherent T6SS activity modifier. The defender reduces this shock damage via competitive immunity scaling logarithmically against the size differential:
  $$ \text{Damage}_{\text{T6SS}} = \left( T_{\text{act}} \times \text{Size}_{\text{attacker}}^2 \times 5.0 \right) \times \left( 1.0 - \min\left(0.95, \frac{I_{\text{factor}} \times \text{Size}_{\text{target}}}{\text{Size}_{\text{attacker}} + 0.01}\right) \right) $$

## 7. Performance Architecture (Large $N$ Scale)

To sustain extremely large unit counts ($N > 5000$) without suffering from $O(N \times M)$ processing halts or garbage collection spikes, the engine incorporates several aggressive systemic overrides:

### 7.1. Object Pooling
To prevent memory fragmentation and continuous deallocation during binary fission and cell death, both `FoodNode` and `Bacterium` utilize a strict object-pooling paradigm (`_pool`, `acquire()`, `release()`). When an entity dies, its memory sector is scrubbed via `.init()` and placed in a suspended ring-buffer rather than destructing.

### 7.2. Spatial Hashing ($O(1)$ Lookup)
Iterating physical food searches across the entire map ($N \times M$ Manhattan distances) destroys CPU throughput. The engine generates a discretized `FoodSpatialHash` grid mapping clustered subsets of coordinates:
$$ \text{Cell}_{\text{index}} = \left( \lfloor \frac{x}{\text{cell\_size}} \rfloor, \lfloor \frac{y}{\text{cell\_size}} \rfloor \right) $$
Pathfinding algorithms locally query this hash matrix rather than validating global coordinate spans.

### 7.3. Hybrid ABM/ODE and Viewport Culling
Entities are visually culled to discard rendering sequences when entirely out of the active user viewport. The backend parallels this using **Hybrid Agent-Based / Ordinary Differential Equation** translation toggles. When entities are flagged as `is_out_of_bounds`, complex stochastic chemotaxis vectors and localized grid collisions are dynamically suppressed. Instead, simple rigid-body kinematics equations are applied for basic translation:
$$ \vec{x}_{t+\Delta t} = \vec{x}_t + \vec{v} \Delta t $$
Where $\vec{v}$ remains continuous from their exit vector, drastically suppressing the math computation limits.
