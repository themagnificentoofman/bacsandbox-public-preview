from pygame.locals import *
import pygame

# Frames per second for the simulation
FPS = 144

# World scale: simulation coordinates are micrometers; this controls scale-bar display only.
WORLD_UNIT_NAME = 'um'
# Scale indicator uses a 100-pixel reference bar because it is legible across common window sizes.
SCALE_BAR_REFERENCE_PIXELS = 100

# Default species values are biologically plausible interactive seeds, not hard limits.
DEFAULT_FOOD_CAPACITY_PG = 10.0
DEFAULT_METABOLIC_RATE_PG_S = 0.1
DEFAULT_CONSUMPTION_SPEED_PG_S = 1.0
# E. coli swimming speeds are commonly reported around 20-40 um/s in liquid media.
DEFAULT_BACTERIA_SPEED_UM_S = 40.0
DEFAULT_GREEDINESS = 0.5
DEFAULT_T6SS_ACTIVITY = 0.5
DEFAULT_BACTERIOCIN_PRODUCTION = 0.5
DEFAULT_IMMUNITY_FACTOR = 0.5
DEFAULT_REPRODUCTION_THRESHOLD = 0.8
DEFAULT_SHARE_THRESHOLD = 0.6
DEFAULT_PREDATION_THRESHOLD = 0.9
DEFAULT_BASE_SIZE_UM = 1.0
DEFAULT_SENSOR_RADIUS_UM = 200.0
# Geometry comparisons need a micron-scale epsilon to avoid division by zero in normalized vectors.
NUMERIC_DISTANCE_EPSILON_UM = 1.0e-3
# Gradient vectors below this magnitude are treated as no chemotactic signal.
CHEMOTAXIS_GRADIENT_EPSILON = 1.0e-3
# Lifespan variance is denominator-protected even if imported legacy data is pathological.
MIN_LIFESPAN_VARIANCE_DENOMINATOR = 1.0e-3

# E. coli run-and-tumble timing uses real elapsed seconds, not frame counts.
BACTERIA_RUN_MEAN_SEC = 1.0
# Standard deviation is 20% of the mean run duration, matching broad single-cell variability.
BACTERIA_RUN_STD_SEC = 0.2
# Bounds prevent nonphysical negative/very long runs after random sampling.
BACTERIA_RUN_MIN_SEC = 0.25
BACTERIA_RUN_MAX_SEC = 3.0
# E. coli tumbles are brief flagellar-bundle disruptions, typically near 0.1 s.
BACTERIA_TUMBLE_DURATION_SEC = 0.1
# Observed post-tumble heading changes average about 68 degrees, approximated as 1.19 radians.
BACTERIA_TUMBLE_MEAN_TURN_RAD = 1.19
# Tumbles vary substantially from cell to cell, so a broad angular spread keeps runs stochastic.
BACTERIA_TUMBLE_TURN_STD_RAD = 0.45
# During tumble, motors reorient rather than propel; translation is strongly reduced but remains forward.
BACTERIA_TUMBLE_TRANSLATION_FACTOR = 0.05
# Positive chemotactic gradients extend runs, matching temporal comparison behavior in E. coli.
BACTERIA_CHEMOTAXIS_RUN_EXTENSION_FACTOR = 2.5
# A 1 um sphere in water at room temperature has rotational diffusion around 0.16 rad^2/s.
BACTERIA_ROTATIONAL_DIFFUSION_COEFFICIENT_RAD2_S = 0.16
# Brownian angular displacement variance is 2 D_rot dt for one rotational degree of freedom.
BROWNIAN_ROTATIONAL_VARIANCE_FACTOR = 2.0
# Near-field grazing increases effective viscous drag; a 50% speed factor keeps cells moving while feeding.
BACTERIA_GRAZING_SPEED_FACTOR = 0.5
# Circular boundary correction nudges cells just inside the wall to avoid repeated floating-point collisions.
BACTERIA_BOUNDARY_OFFSET_UM = 0.1
# Wall impacts include small angular roughness to represent imperfect collisions with dish boundaries.
BACTERIA_BOUNDARY_REFLECTION_NOISE_RAD = 0.2
# Rendered starvation is a subtle visual cue; bacterial cell walls do not elastically shrink by large factors.
BACTERIA_STARVED_RENDER_SCALE = 0.9
BACTERIA_FULL_RENDER_SCALE_RANGE = 0.1
# Individual lifespan variability is modeled as a Gaussian multiplier with a floor to avoid extreme outliers.
LIFESPAN_VARIANCE_MIN = 0.4
LIFESPAN_VARIANCE_MEAN = 1.0
LIFESPAN_VARIANCE_STD = 0.2
# Senescence maintenance doubles over a few simulated minutes to keep aging visible in interactive runs.
SENESCENCE_BASE_LIFESPAN_SEC = 180.0
# Gompertz alpha sets mortality doubling near 17 seconds before individual lifespan scaling.
GOMPERTZ_MORTALITY_ALPHA = 0.04
# Stouthamer-style biomass synthesis cost in glucose-equivalent pg per 1 um^3 child volume.
FISSION_COST_PG_PER_UM3 = 40.0

# Softened receptor field avoids singular food signals at zero distance while preserving nearby gradients.
CHEMOTAXIS_CORE_RADIUS_UM = 5.0
# Reference concentration for Weber-Fechner log receptor compression in glucose-equivalent pg/um^2.
CHEMOTAXIS_REFERENCE_PG_PER_UM2 = 1.0e-3
# Low-Reynolds swimming reorients gradually toward sensed gradients during runs.
CHEMOTAXIS_STEER_RATE_RAD_S = 1.5

# Wet bacterial mass is approximated as 1 pg/um^3; dry structural biomass is about 30% of wet mass.
CELL_WET_MASS_PG_PER_UM3 = 1.0
CELL_DRY_MASS_FRACTION = 0.30
# Only labile necromass is immediately scavengable; the rest is refractory, respired, or diffuses away.
LABILE_NECROMASS_FRACTION = 0.35
# Remaining cytosolic nutrient can leak at death, but is not recovered perfectly.
CYTOSOL_RECOVERY_FRACTION = 0.75
# Food nodes below this amount are numerically and biologically negligible.
MIN_FOOD_NODE_AMOUNT_PG = 1.0e-3
# Food-node amount is treated as a spherical glucose-equivalent volume for contact geometry.
FOOD_NODE_MIN_RADIUS_UM = 0.5
FOOD_NODE_MIN_VOLUME_PG = 0.1
# Bacteria contact a nutrient patch at radius + cell size plus a near-field margin.
FOOD_EATING_DISTANCE_MARGIN_UM = 1.0
# Uptake can occur just outside contact radius because dissolved molecules diffuse locally.
FOOD_EATING_DISTANCE_GRACE_UM = 3.0

# pH fitness uses a Gaussian tolerance of 1.5 pH units around the species optimum.
PH_TOLERANCE_GAUSSIAN_SD = 1.5
# Aerobic fitness saturates by 10% O2; strict anaerobes are strongly inhibited around 2% O2.
AEROBE_O2_SATURATION_LEVEL = 0.10
ANAEROBE_O2_TOXICITY_THRESHOLD = 0.02
# Toxic death under poor environmental fitness remains a low per-second background hazard.
ENV_TOXIC_DEATH_RATE_PER_SEC = 0.05

# Ribosomally synthesized bacteriocins cost nutrient to produce and create larger local ECM amounts.
BACTERIOCIN_PRODUCTION_COST_RATE = 0.5
BACTERIOCIN_ECM_SECRETION_MULTIPLIER = 10.0
# Immunity acts as competitive mitigation but does not make cells absolutely invulnerable.
BACTERIOCIN_MAX_IMMUNITY_MITIGATION = 0.99
BACTERIOCIN_IMMUNITY_SCALING = 1.5
BACTERIOCIN_DAMAGE_COEFFICIENT = 0.1
# T6SS damage scales with attacker surface area and is capped by target immunity.
T6SS_DAMAGE_SCALING_FACTOR = 5.0
T6SS_MAX_IMMUNITY_DEFENSE = 0.95
# T6SS contact uses a micron-scale soft margin because cells are rendered as disks in 2D.
BACTERIA_CONTACT_MARGIN_UM = 2.0
# Epsilon prevents singular immunity ratios when attacker size approaches zero.
T6SS_IMMUNITY_SIZE_EPSILON_UM = 0.01
# Phagocytosis requires a 2x linear size advantage, equivalent to about 8x volume.
PHAGOCYTOSIS_SIZE_RATIO_THRESHOLD = 2.0
# Consumed prey are marked below starvation level so the species step removes them without necromass double-drop.
CONSUMED_PREY_SENTINEL_FOOD_LEVEL = -9999.0

# Environment is modeled as 1 mL of fluid; these are dissolved gas capacities in pg at model saturation.
ENV_100_O2_PG = 4.0e7
ENV_100_CO2_PG = 1.5e9
# Aerobic respiration stoichiometry: glucose + 6 O2 -> 6 CO2 + 6 H2O.
GLUCOSE_MOLAR_MASS_G = 180.16
O2_SIX_MOLES_MASS_G = 192.0
CO2_SIX_MOLES_MASS_G = 264.0
# Carbonic acid and lactic acid terms are simplified pH effects of CO2 and fermentation products.
CARBONIC_ACID_PH_EFFECT = 0.02
HETEROLACTIC_CO2_TRACE_FRACTION = 0.1
LACTIC_ACID_PH_EFFECT = 2.5
# Open shallow dishes equilibrate gas/pH toward atmosphere/neutral buffer on a seconds scale.
GAS_DIFFUSION_RATE_SEC = 0.288
BUFFER_RESTORATION_RATE_SEC = 0.144
# Shallow open dishes equilibrate toward atmospheric O2, near-zero excess CO2, and neutral pH.
ATMOSPHERIC_O2_FRACTION = 0.21
ATMOSPHERIC_CO2_FRACTION = 0.0
NEUTRAL_PH = 7.0

# Explicit Fickian diffusion stencil is nonnegative in 2D when alpha = D dt / dx^2 <= 0.25.
ECM_ALPHA_MAX = 0.24
# Quorum-sensing molecules diffuse faster than bacteriocins and decay on a minutes-scale half-life.
ECM_QS_DIFFUSION_UM2_S = 400.0
ECM_QS_DECAY_S = 0.00385
# Bacteriocins are larger peptides/proteins, so they diffuse slower and decay slightly faster.
ECM_BAC_DIFFUSION_UM2_S = 60.0
ECM_BAC_DECAY_S = 0.00578
# ECM cells below this concentration are pruned to prevent unbounded sparse-grid growth.
ECM_MIN_CONCENTRATION = 1.0e-3
# 20 um ECM cells keep explicit diffusion stable at interactive tick rates while preserving local gradients.
ECM_CELL_SIZE_UM = 20

# Spatial hash cell sizes balance local lookup speed with bacterial interaction radii.
FOOD_SPATIAL_HASH_CELL_SIZE_UM = 50
BACTERIA_INTERACTION_HASH_CELL_SIZE_UM = 20
BACTERIA_INTERACTION_AXIS_CUTOFF_UM = 15
# Petri dish padding keeps the initial circular boundary outside the visible window edge.
PETRI_RADIUS_PADDING_UM = 100

# Background color (RGB)
BACKGROUND = (30, 30, 30)
# Color for text and certain UI elements
WHITE = (255, 255, 255)
# Button base/hover colors
BUTTON_COLOR = (70, 70, 200)
BUTTON_HOVER = (90, 90, 220)
# Slider color
SLIDER_COLOR = (200, 200, 70)
# Filenames for saving/exporting data
EXPORT_FILE = 'simulation_stats.csv'
DATA_FILE = 'species_data.csv'
ENV_CONFIG_FILE = 'env_config.json'
# Font placeholder (initialized later)
FONT = None

