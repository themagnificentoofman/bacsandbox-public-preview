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
# ATP/biomass synthesis cost in glucose-equivalent pg per 1 um^3 child volume; see DOCUMENTATION.md §4.7 (Reproduction).
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

# --- Phase I: Diet-driven food substrates (malnutrition-research grounding) ---
# Three discrete dietary substrate classes let the environment reproduce specific diets
# (malnutrition / standard supplementary / MDCF) instead of one abstract glucose pool. The
# integer IDs index per-species affinity arrays and the per-substrate yield/colour tables.
SUBSTRATE_GENERIC_CALORIE = 0   # simple, universally fermentable sugar (e.g. rice/starch hydrolysate)
SUBSTRATE_COMPLEX_GLYCAN = 1    # dietary fibre / complex plant glycan (partially fermentable)
SUBSTRATE_MDCF_PREBIOTIC = 2    # Microbiota-Directed Complementary Food prebiotic (selective)
SUBSTRATE_IDS = (SUBSTRATE_GENERIC_CALORIE, SUBSTRATE_COMPLEX_GLYCAN, SUBSTRATE_MDCF_PREBIOTIC)
SUBSTRATE_NAMES = ('generic_calorie', 'complex_glycan', 'mdcf_prebiotic')
# Legacy, necromass, and decomposition food all default to the generic calorie class so prior
# behaviour (and old saves) is reproduced exactly when no diet is configured.
DEFAULT_SUBSTRATE_ID = SUBSTRATE_GENERIC_CALORIE

# Metabolic energy yield per pg ingested, expressed as a fraction of pure glucose (1.0).
# Simple sugars are fully oxidised (1.0); complex glycans lose energy to incomplete fermentation
# and host-indigestible residue (~0.7); MDCF prebiotics are targeted, efficiently fermentable
# substrates (~0.85). see DOCUMENTATION.md §4.19 (Diet-Driven Substrate Metabolism).
SUBSTRATE_ENERGY_YIELD = {
    SUBSTRATE_GENERIC_CALORIE: 1.0,
    SUBSTRATE_COMPLEX_GLYCAN: 0.7,
    SUBSTRATE_MDCF_PREBIOTIC: 0.85,
}
# Render colours (RGB) so the three substrate classes are visually distinguishable on the dish.
SUBSTRATE_COLORS = {
    SUBSTRATE_GENERIC_CALORIE: (200, 200, 70),   # yellow-green (legacy food colour)
    SUBSTRATE_COMPLEX_GLYCAN: (210, 140, 50),    # amber/orange (fibre)
    SUBSTRATE_MDCF_PREBIOTIC: (140, 90, 200),    # violet (engineered prebiotic)
}
# Default dimensionless per-species uptake-capacity multipliers. The legacy node path scales
# consumption_speed, while the reduced model scales each substrate's base Vmax; zero disables
# uptake in either path, and legacy species remain calorie-only; see DOCUMENTATION.md section 10.4.
DEFAULT_AFF_CALORIE = 1.0
DEFAULT_AFF_GLYCAN = 0.0
DEFAULT_AFF_MDCF = 0.0

# Diet profiles: categorical probability that each fresh nutrient drop is of a given substrate
# class. These ratios encode the three experimental feeding arms. Malnutrition-like diets are
# low-diversity and calorie-dominated; standard supplementary food adds fibre; the MDCF arm is
# enriched for the targeted prebiotic. Fractions per profile sum to 1.0. Values are calibration
# estimates informed by the MDCF trial design and are flagged for refinement against empirical
# dietary-intake data. see DOCUMENTATION.md §4.19 (Diet-Driven Substrate Metabolism).
DIET_PROFILES = {
    'malnutrition': {SUBSTRATE_GENERIC_CALORIE: 0.85, SUBSTRATE_COMPLEX_GLYCAN: 0.10, SUBSTRATE_MDCF_PREBIOTIC: 0.05},
    'standard':     {SUBSTRATE_GENERIC_CALORIE: 0.60, SUBSTRATE_COMPLEX_GLYCAN: 0.30, SUBSTRATE_MDCF_PREBIOTIC: 0.10},
    'mdcf':         {SUBSTRATE_GENERIC_CALORIE: 0.40, SUBSTRATE_COMPLEX_GLYCAN: 0.25, SUBSTRATE_MDCF_PREBIOTIC: 0.35},
}
DEFAULT_DIET_PROFILE = 'standard'

# --- Phase I: Spatial heterogeneity (gut lumen vs mucus layer) ---
# The petri dish is geometrically partitioned into a central, well-mixed LUMEN and a thin outer
# MUCUS annulus. A point is in the mucus layer when its radial position (as a fraction of the dish
# radius) is at or beyond this threshold; the inner 75% is lumen, the outer 25% is mucus. This
# mirrors the radial organisation of the intestinal cross-section (concentric lumen/mucus zones).
# see DOCUMENTATION.md §4.20 (Spatial Heterogeneity: Lumen vs Mucus Layer)
MUCUS_LAYER_RADIUS_FRACTION = 0.75
# Mucus-associated cells adhere to the glycoprotein matrix and move much more slowly than free
# luminal swimmers, creating an anchoring refuge. Motility is reduced to 30% in the mucus layer.
# see DOCUMENTATION.md §4.20 (Spatial Heterogeneity: Lumen vs Mucus Layer)
MUCUS_MOTILITY_FACTOR = 0.3
# Nutrient washout (peristaltic / luminal-flow clearance) of free food mass, per real second.
# The lumen is rapidly flushed; the mucus layer is a stabilized reservoir with near-zero
# clearance. Exponential decay: m(t+dt) = m(t) * (1 - rate*dt). Conservative default rates.
LUMEN_FOOD_WASHOUT_RATE_PER_SEC = 0.02
MUCUS_FOOD_WASHOUT_RATE_PER_SEC = 0.002
# Fraction of a nutrient drop's mass retained when it is deposited in the mucus layer; the lumen
# (factor 1.0) receives the bulk of fresh dietary substrate, the mucus annulus far less.
MUCUS_DROP_AMOUNT_FACTOR = 0.5
# Relative carrying-capacity cull weight for luminal cells versus mucus cells (weight 1.0).
# Stronger luminal competition makes well-mixed lumen cells more likely to be lost under
# population pressure, leaving the mucus a protective refuge.
LUMEN_COMPETITION_CULL_WEIGHT = 3.0

# --- Initial inoculum spawning (founder placement at t = 0) ---
# Real microbial communities are not seeded as a uniform sheet: an inoculum arrives as a small
# number of founder cells that establish localised microcolonies (clonal patches) which then grow
# clonally by binary fission (daughters placed at the parent, see attempt_fission). These constants
# parameterise that founder placement so the starting layout is biologically motivated, petri-dish
# aware, and fully reproducible under the RNG seed.
#
# Spawn mode selects the macroscopic geometry of the founder distribution:
#   'clustered'    -> a few founder nuclei, cells scattered in Gaussian microcolonies around them
#                     (default; mimics sparse colonisation and seeds spatial structure for
#                      competition / cross-feeding to act on).
#   'uniform_disc' -> area-uniform random placement across the dish (r = R*sqrt(U) sampling so
#                     density is constant per unit area, not per unit radius).
#   'central'      -> a single compact droplet near the dish centre (point-inoculation assay).
INITIAL_SPAWN_MODE = 'clustered'
# Fraction of founder cells placed inside the central lumen (radial fraction < the lumen boundary
# below) rather than the outer mucus annulus. Fresh dietary substrate is delivered chiefly to the
# well-mixed lumen, so a gut-realistic inoculum is lumen-biased; 1.0 = lumen only, 0.0 = no bias.
INITIAL_SPAWN_LUMEN_BIAS = 0.9
# Radial fraction (of the dish radius) marking the lumen/mucus boundary used by the lumen bias.
# Reuses the established lumen geometry so spawning and the spatial model stay consistent.
INITIAL_SPAWN_LUMEN_RADIUS_FRACTION = MUCUS_LAYER_RADIUS_FRACTION
# Number of founder nuclei (microcolony seed points) generated per species in 'clustered' mode.
# A small count yields a few distinct colonies; biologically this is the number of cells that
# survive transit and successfully attach to establish clonal patches.
INITIAL_SPAWN_CLUSTER_COUNT = 3
# Gaussian spread of a microcolony, expressed as a fraction of the dish radius; this is the
# standard deviation of the offset applied to each founder cell around its nucleus. Small spread
# keeps colonies tight and locally clonal.
INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION = 0.06
# Radial fraction (of the dish radius) of the compact droplet used in 'central' mode; the whole
# inoculum is confined within this radius of the dish centre (area-uniform within the droplet).
INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION = 0.15
# Inter-species arrangement controls how multiple species are positioned relative to one another:
#   'independent' -> each species draws its own independent founder nuclei (default; colonies of
#                    different species partially overlap by chance, the generic well-mixed case).
#   'mixed'       -> all species share the same founder nuclei, co-locating colonies to maximise
#                    early contact-dependent interactions (T6SS, bacteriocins, HGT, cross-feeding).
#   'segregated'  -> each species is confined to its own angular sector of the dish, separating
#                    founders to study invasion / range-expansion dynamics from distinct origins.
INITIAL_SPECIES_ARRANGEMENT = 'independent'

# Canonical option lists for the two categorical spawn settings exposed in the Globals tab. Defined
# here so the UI cycler (app.py) and the simulation validator (core_sim._apply_spawn_env) share one
# source of truth; the first entry of each is the compiled-in default above.
SPAWN_MODES = ('clustered', 'uniform_disc', 'central')
SPECIES_ARRANGEMENTS = ('independent', 'mixed', 'segregated')

# pH fitness uses a Gaussian tolerance of 1.5 pH units around the species optimum.
PH_TOLERANCE_GAUSSIAN_SD = 1.5
# Aerobic fitness saturates by 10% O2; strict anaerobes are strongly inhibited around 2% O2.
AEROBE_O2_SATURATION_LEVEL = 0.10
ANAEROBE_O2_TOXICITY_THRESHOLD = 0.02

# --- Temperature kinetics: Ratkowsky square-root model ---
# see DOCUMENTATION.md §4.14 (Temperature Kinetics)
# Full (four-parameter) form: sqrt(mu) = b*(T - T_min)*(1 - exp(c*(T - T_max))).
#   mu       = specific growth rate (here used as a dimensionless rate multiplier after normalization)
#   T_min    = theoretical minimum temperature below which no growth occurs
#   T_max    = theoretical maximum temperature above which no growth occurs
#   b        = slope parameter governing the low-temperature square-root rise
#   c        = Arrhenius-like steepness of high-temperature thermal inactivation
# Because the simulation uses the curve NORMALIZED to each species' own optimum, b cancels
# out of the multiplier; it is retained for documentation and possible absolute-rate use.
# Default b/c approximate mesophilic (E. coli-like) growth kinetics.
RATKOWSKY_B_SQRT_PER_C = 0.0410      # sqrt(rate) slope per degC; mesophile literature ~0.035-0.041
RATKOWSKY_C_PER_C = 0.161            # high-temperature inactivation steepness per degC (E. coli-like)
# Default cardinal temperatures: E. coli grows ~4-49 degC; defaults bracket a mesophile range.
DEFAULT_TEMP_MIN_C = 4.0             # theoretical minimum growth temperature (degC)
DEFAULT_TEMP_MAX_C = 47.0            # theoretical maximum growth temperature (degC)
# Default incubation temperature; 37 degC is standard human/host body temperature for E. coli.
DEFAULT_ENV_TEMPERATURE_C = 37.0
# Environment temperature slider bounds (degC); spans deep-freeze to autoclave-adjacent extremes.
ENV_TEMPERATURE_MIN_C = 0.0
ENV_TEMPERATURE_MAX_C = 60.0
# One-time per-species sample count used to locate the Ratkowsky growth optimum for normalization.
RATKOWSKY_NORMALIZATION_SAMPLES = 200

# --- Lag phase: Baranyi dynamic physiological-state model ---
# A freshly inoculated cell carries an internal "work to be done" before it can divide, captured by a
# dimensionless physiological state q. The adjustment function alpha(q) = q / (1 + q) multiplies the
# growth rate: alpha ~ 0 during lag (q small), rising to 1 once the cell has adapted (q large).
# see DOCUMENTATION.md §4.16 (Lag Phase)
DEFAULT_BARANYI_Q0 = 0.05  # initial physiological state of fresh inoculum; h0 = ln(1 + 1/q0) ~ 3.04 (moderate lag)
# q grows as dq/dt = mu_max * q (here mu_max is the species' temperature-scaled specific growth rate).
# Cap q to keep alpha numerically at ~1 without risking float overflow on long runs; adaptation is
# biologically complete (alpha > 0.999999) far below this bound.
BARANYI_Q_STATE_MAX = 1.0e6
BARANYI_Q0_MIN = 1.0e-4  # floor so lag is finite (h0 <= ln(1 + 1e4) ~ 9.2); zero q0 would never adapt
BARANYI_Q0_MAX = 1.0e3   # ceiling: high q0 means effectively no lag (alpha starts near 1)

# --- Phenotypic switching: persister / proliferative bet-hedging (M4) ---
# A subpopulation stochastically switches into a dormant, stress-tolerant "persister" state that
# halts division, minimizes maintenance metabolism, and resists antimicrobials and senescence. The
# switch is reversible: persisters resuscitate when nutrients return. Switching is governed by two
# fullness thresholds (a starvation trigger S1 and a recovery trigger S2) and per-second rates.
# Rates are sim-accelerated (the simulation runs in real seconds, not hours) so that the
# bet-hedging dynamics are observable within a session while preserving the relative ordering of the
# spontaneous, stress-induced, and resuscitation pathways. see DOCUMENTATION.md §4.17 (Phenotypic Switching)
DEFAULT_PERSISTER_SWITCH = 0.0       # per-species switching propensity in [0,1]; 0 disables persistence (pre-M4 behaviour)
DEFAULT_PERSISTER_ENTER_THRESH = 0.15  # S1: internal fullness below which starvation-induced switching engages
DEFAULT_PERSISTER_EXIT_THRESH = 0.50   # S2: internal fullness above which persisters resuscitate to proliferative
PERSISTER_SPONTANEOUS_SWITCH_RATE_PER_SEC = 0.002  # baseline stochastic normal->persister (bet-hedging) rate
PERSISTER_STRESS_SWITCH_RATE_PER_SEC = 0.2         # starvation-induced normal->persister rate (stringent/SOS response)
PERSISTER_WAKE_RATE_PER_SEC = 0.2                  # nutrient-triggered persister->normal resuscitation rate
PERSISTER_SPONTANEOUS_WAKE_RATE_PER_SEC = 0.005    # baseline stochastic persister->normal reversion rate
PERSISTER_METABOLIC_FACTOR = 0.1     # dormant maintenance burn as a fraction of basal metabolism (dormant, not dead)
PERSISTER_RESISTANCE_FACTOR = 0.1    # incoming bacteriocin/T6SS/senescence hazard scaled to 10% (~10x tolerance)
PERSISTER_THRESHOLD_MIN_GAP = 0.05   # minimum S2 - S1 separation so enter/exit bands cannot overlap and chatter
PERSISTER_RING_COLOR = (40, 40, 50)  # dark desaturated outline marking dormant persister cells in the renderer

# --- Lytic bacteriophage predation + horizontal gene transfer (M5) ---
# Free virions are tracked as a diffusible ECM field. Susceptible cells adsorb phage by mass action,
# enter a latent (eclipse) period, then lyse, releasing a burst of progeny virions. A heritable
# resistant subpopulation (arising by rare mutation or by horizontal transfer) survives the crash and
# rebounds, producing the classic predator-prey / SIR boom-bust dynamics. The model couples one-step
# growth (latent period + burst size), mass-action adsorption / predator-prey coupling, horizontal
# gene transfer, and a plasmid/resistance carriage cost. Rates are sim-accelerated (real seconds, not
# the ~20-40 min latent periods of nature) for observable dynamics.
# see DOCUMENTATION.md §4.18 (Phage Predation and Horizontal Gene Transfer)
ECM_PHAGE_DIFFUSION_UM2_S = 500.0      # virion diffusion coefficient; lower than small solutes (larger particles)
ECM_PHAGE_DECAY_S = 0.02               # virion inactivation (loss of infectivity) rate per second
PHAGE_ECM_UNITS_PER_VIRION = 1.0       # ECM grid units per virion; kept O(1) so a single adsorption is meaningful
DEFAULT_PHAGE_ADSORPTION_RATE = 0.02   # mass-action adsorption rate constant per virion-unit per second
DEFAULT_PHAGE_LATENT_PERIOD_SEC = 3.0  # eclipse + latent period before lysis (sim-accelerated)
DEFAULT_PHAGE_BURST_SIZE = 20.0        # progeny virions released per lysed cell (lytic burst, scaled from ~50-200)
DEFAULT_PHAGE_INFLOW = 0.0             # global free-phage introduction rate (seed packets/sec); 0 = no outbreak
PHAGE_SEED_PACKET_VIRIONS = 10.0       # virions per introduced seed packet (enough to exceed the adsorbable floor)
PHAGE_MIN_ADSORBABLE = 1.0e-3          # local phage below this (= ECM prune floor) is negligible for adsorption
PHAGE_RESISTANCE_MUTATION_RATE = 0.001 # spontaneous susceptible->resistant mutation probability per division
DEFAULT_PHAGE_RESISTANCE = 0.0         # per-species initial resistant fraction in [0,1] (0 = fully susceptible)
DEFAULT_HGT_COMPETENCE = 0.0           # per-species HGT donate/receive competence in [0,1]; 0 disables transfer
HGT_TRANSFER_RATE_PER_SEC = 0.1        # contact-dependent horizontal transfer rate constant (conjugation-like)
HGT_METABOLIC_BURDEN = 0.2             # fractional rise in basal metabolism from carrying the transferred plasmid
PHAGE_INFECTED_RING_COLOR = (255, 140, 0)  # orange outline marking phage-infected (lysing) cells
PHAGE_RESISTANT_RING_COLOR = (0, 200, 220)  # cyan outline marking phage-resistant cells
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

# --- Cross-feeding metabolite: acetate overflow (the "acetate switch") ---
# see DOCUMENTATION.md §4.15 (Metabolic Cross-Feeding)
# During fast aerobic growth on glucose, E. coli excretes overflow acetate; other taxa oxidize it as a
# secondary carbon source. We model acetate as a third diffusible ECM layer linking producers to consumers.
ACETATE_MOLAR_MASS_G = 60.05  # acetic acid / acetate, g/mol (CH3COOH)
# Overflow upper bound: 1 glucose can be split into 2 acetate by mass (glucose -> 2 acetate, ignoring CO2 loss).
ACETATE_PER_GLUCOSE_MAX_MASS_RATIO = (2.0 * ACETATE_MOLAR_MASS_G) / GLUCOSE_MOLAR_MASS_G  # ~0.6667 g acetate / g glucose
# ATP yields set the relative energy value of acetate vs glucose so cross-feeders gain less per pg.
GLUCOSE_AEROBIC_ATP_PER_MOL = 32.0  # ~30-32 ATP/glucose, full aerobic oxidation; see DOCUMENTATION.md §4.15 (Metabolic Cross-Feeding)
ACETATE_AEROBIC_ATP_PER_MOL = 10.0  # ~10 ATP/acetate via acetyl-CoA -> TCA + ETC
# Energy per unit mass relative to glucose: (ATP/g acetate) / (ATP/g glucose).
ACETATE_ENERGY_YIELD_GLUCOSE_EQUIV = (
    (ACETATE_AEROBIC_ATP_PER_MOL / ACETATE_MOLAR_MASS_G)
    / (GLUCOSE_AEROBIC_ATP_PER_MOL / GLUCOSE_MOLAR_MASS_G)
)  # ~0.9376 glucose-energy-equivalent per pg acetate consumed
# Acetate oxidation stoichiometry: C2H4O2 + 2 O2 -> 2 CO2 + 2 H2O.
O2_PER_ACETATE_MASS_RATIO = 64.0 / ACETATE_MOLAR_MASS_G   # 2 mol O2 (64.0 g) per mol acetate
CO2_PER_ACETATE_MASS_RATIO = 88.02 / ACETATE_MOLAR_MASS_G  # 2 mol CO2 (88.02 g) per mol acetate
# Acetate is a small ion that diffuses quickly in water and washes out slowly by dilution.
ECM_ACETATE_DIFFUSION_UM2_S = 1000.0  # ~1.09e3 um^2/s aqueous acetate self-diffusion; see DOCUMENTATION.md §4.15 (Metabolic Cross-Feeding)
ECM_ACETATE_DECAY_S = 0.001  # slow abiotic dilution/washout half-life
# Picogram acetate is tiny; scale into ECM grid units so it survives the ECM_MIN_CONCENTRATION prune floor.
# The same factor is inverted on uptake so secreted/consumed mass is conserved.
ACETATE_ECM_UNITS_PER_PG = 100.0
# Backward-compatible defaults: 0 secretion and 0 uptake reproduce pre-M2 behaviour exactly.
DEFAULT_ACETATE_SECRETION = 0.0
DEFAULT_ACETATE_UPTAKE = 0.0

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

