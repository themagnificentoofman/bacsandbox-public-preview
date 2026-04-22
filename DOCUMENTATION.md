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
The application features a real-time UI driven by PyGame, separated into tabs:
* **Grid**: The primary simulation view. Includes dynamic panning, zooming (using mouse or on-screen buttons), and a minimap.
* **Species Editor**: A deeply interactive data-grid for creating, modifying, duplicating, or deleting species parameters on the fly.
* **Globals Editor**: Allows tuning of environmental factors (O2, CO2, pH, Capacity) and food node drop frequencies.
* **Stats**: Tracks populations, births, and deaths globally and per-species.
* **Graph & Data Management**: Real-time graphing and the ability to export (`simulation_stats.csv`) and import simulation data dynamically.
