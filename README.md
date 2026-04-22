# Bacteria Sandbox

Welcome to **Bacteria Sandbox**, a robust 2D simulation built with Python and Pygame!

## Project Overview
Bacteria Sandbox is an interactive grid-based simulation where various bacterial species vie for resources, multiply, and perish according to configurable environmental pressures and species-specific traits.

For a detailed breakdown of the implemented mathematical biology models, tunable species traits, physics kinematics, and UI features, please see the formally curated [DOCUMENTATION.md](DOCUMENTATION.md).

## Setup & Installation
1. Ensure you have Python 3 installed.
2. Install the required packages via `requirements.txt`:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the application:
   ```bash
   python -m bacteria
   ```

## User Interface
The application features a real-time UI driven by PyGame, separated into tabs:
- **Grid Tab**: A huge 65536x65536 scrolling space for the primary simulation view. Use your mouse to drag the map or the bottom/right scrollbars to navigate, and the minimap for quick movement.
- **Graph Tab**: Real-time line charts tracking historical total populations and per-species populations. Includes export and import commands to save or replay timeline stats securely via CSV.
- **Stats Tab**: Displays numerical readouts of live simulation events tracking active populations, global births, and total deaths globally and per-species.
- **Species Tab**: An interactive data-grid for creating, modifying, duplicating, or deleting species parameters on the fly. Also features Undo/Redo tracking.
- **Globals Tab**: An embedded spreadsheet editor defining baseline environmental factors (O2, CO2, pH, Capacity) and food node drop frequencies.
- **Environment Sliders**: Use the UI sliders visible in the Grid, Graph, and Stats views to manually adjust live environmental features such as the amount of O2, CO2, Food Drop Frequencies and Range, as well as the pH scale in real-time.

## Configuration Files
The project maintains its state using a couple of configuration files:

### `species_data.csv`
Defines the parameters for every configured species. Can be edited externally or in-app.
* **name**: Identifier for the species.
* **init**: Initial population spawned.
* **rep**: Base reproduction rate.
* **aging**: Base aging rate multiplier.
* **so**: O2 sensitivity (-1 to 1; <0 needs CO₂, >0 needs O₂).
* **sf**: Food sensitivity (0 to 1).
* **speed**: Base movement speed.
* **ph**: Preferred pH level (0-14).
* **food_capacity**: Base max food (pg Glucose, scales w/ $\mu m^3$).
* **metabolic_rate**: Base food used (pg/tick, scales w/ Kleiber).
* **consumption_speed**: Base eat speed (pg/tick, scales w/ $\mu m^2$).
* **greediness**: (0.0 to 1.0) Likelihood to prioritize distant larger food piles over close smaller ones.
* **aggression**: Toxin power against other species (0-1).
* **defense**: Toxin resistance (0-1).
* **rep_thresh**: Min food capacity % to replicate (0-1).
* **share_thresh**: Min food capacity % to share (0-1).
* **predation_thresh**: Max food capacity % to predate (0-1).
* **base_size**: Base size ($\mu m$, default 1.0).
* **sensor_radius**: Food detection radius ($\mu m$, default 200.0).

### `env_config.json`
Stores the global simulation state:
* `o2`: Active O2 percent (0-1).
* `co2`: Active CO2 percent (0-1).
* `ph`: Current pH level (0-14).
* `capacity`: Hard carrying capacity for all species combined.
* `drop_frequency`: Amount of new food clusters to attempt to spawn per grouping tick.
* `drop_range_min` / `drop_range_max`: Min/Max amount of food in newly dropped clusters.

## Data Export
As the simulation runs, statistical points matching population changes and event counts are recorded in `simulation_stats.csv`. This includes a detailed log of births, deaths, and running populations over time, great for external data processing or spreadsheet imports.
