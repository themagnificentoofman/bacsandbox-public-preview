import copy
import math, random
import sys
import os
import csv

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'

from .constants import *
from .ecm import ECM
from .metabolic_integration import DietBolus, FbaPhysiology
from .physics import InteractionSummary, PhysicsOrchestrator

def _spawn_sector(species_index, species_total):
    """Angular wedge [a0, a1) a species may occupy. Under the 'segregated' arrangement each species
    is confined to an equal 2*pi / species_total slice anchored at its index, separating founder
    origins so invasion / range-expansion dynamics start from distinct sectors. All other
    arrangements return the full circle."""
    if INITIAL_SPECIES_ARRANGEMENT == 'segregated' and species_total > 1:
        width = 2.0 * math.pi / species_total
        a0 = species_index * width
        return a0, a0 + width
    return 0.0, 2.0 * math.pi

def _spawn_radial_band(radius):
    """Pick the [r_min, r_max] radial band a founder is drawn from, applying the lumen bias: with
    probability INITIAL_SPAWN_LUMEN_BIAS the cell goes in the inner lumen (radii below the lumen
    boundary), otherwise the outer mucus annulus. Biasing the inoculum toward the well-mixed lumen
    matches where fresh dietary substrate is delivered."""
    lumen_r = max(0.0, min(INITIAL_SPAWN_LUMEN_RADIUS_FRACTION, 1.0)) * radius
    if random.random() < max(0.0, min(INITIAL_SPAWN_LUMEN_BIAS, 1.0)):
        return 0.0, lumen_r
    return lumen_r, radius

def _sample_disc_point(cx, cy, r_min, r_max, a0, a1):
    """Sample an area-uniform point in the annular sector (r in [r_min, r_max], theta in [a0, a1]).
    Using r = sqrt(U*(r_max^2 - r_min^2) + r_min^2) makes the areal density constant instead of
    clumping points near the centre (which a naive uniform-in-radius draw would do)."""
    u = random.random()
    r = math.sqrt(u * (r_max * r_max - r_min * r_min) + r_min * r_min)
    th = random.uniform(a0, a1)
    return (cx + r * math.cos(th), cy + r * math.sin(th))

def generate_cluster_nuclei(count, center, radius, species_index=0, species_total=1):
    """Generate `count` microcolony seed points (founder nuclei) inside the dish, honouring the
    lumen bias and the species' angular sector. Biologically these are the few cells that survive
    transit and attach to nucleate clonal patches."""
    cx, cy = center
    a0, a1 = _spawn_sector(species_index, species_total)
    nuclei = []
    for _ in range(max(1, int(count))):
        r_min, r_max = _spawn_radial_band(radius)
        nuclei.append(_sample_disc_point(cx, cy, r_min, r_max, a0, a1))
    return nuclei

def generate_initial_positions(count, center, radius, species_index=0, species_total=1, nuclei=None):
    """Produce `count` founder (x, y) positions for one species under the configured spawn policy.

    All draws come from the module `random` stream, so seeding the RNG before calling this yields a
    reproducible layout. Modes (INITIAL_SPAWN_MODE):
      - 'central'      : compact area-uniform droplet within INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION
                         of the dish radius (point-inoculation assay).
      - 'uniform_disc' : area-uniform placement, lumen-biased band per cell.
      - 'clustered'    : Gaussian microcolonies around founder nuclei (the default). For the 'mixed'
                         arrangement the caller passes shared `nuclei` so species co-locate.
    Every position is constrained to lie within the petri-dish radius."""
    cx, cy = center
    radius = max(NUMERIC_DISTANCE_EPSILON_UM, float(radius))
    count = max(0, int(count))
    if count == 0:
        return []
    a0, a1 = _spawn_sector(species_index, species_total)

    if INITIAL_SPAWN_MODE == 'central':
        # Single compact droplet near the dish centre; well inside the lumen by construction.
        cr = max(0.0, min(INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION, 1.0)) * radius
        return [_sample_disc_point(cx, cy, 0.0, cr, a0, a1) for _ in range(count)]

    if INITIAL_SPAWN_MODE == 'uniform_disc':
        out = []
        for _ in range(count):
            r_min, r_max = _spawn_radial_band(radius)
            out.append(_sample_disc_point(cx, cy, r_min, r_max, a0, a1))
        return out

    # 'clustered' (default; also the fallback for any unrecognised mode).
    if not nuclei:
        nuclei = generate_cluster_nuclei(INITIAL_SPAWN_CLUSTER_COUNT, center, radius,
                                         species_index, species_total)
    spread = max(0.0, INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION) * radius
    out = []
    for i in range(count):
        # Round-robin founders across nuclei for an even, deterministic colony-size split.
        nx, ny = nuclei[i % len(nuclei)]
        x, y = nx, ny
        # Gaussian offset forms the microcolony; resample a few times so no cell escapes the dish.
        for _attempt in range(8):
            cxp = random.gauss(nx, spread)
            cyp = random.gauss(ny, spread)
            if math.hypot(cxp - cx, cyp - cy) <= radius:
                x, y = cxp, cyp
                break
        out.append((x, y))
    return out

class FoodNode:
    """
    Represents a food node on the grid.
    """
    _pool = []
    
    def __init__(self, x, y, amount, substrate_id=DEFAULT_SUBSTRATE_ID):
        self.x = x
        self.y = y
        self.amount = amount
        # Phase I: dietary substrate class (SUBSTRATE_* id) carried by this nutrient patch.
        self.substrate_id = substrate_id
    
    @classmethod
    def acquire(cls):
        """Acquire a FoodNode from pool or create new."""
        if cls._pool:
            node = cls._pool.pop()
            return node
        return cls.__new__(cls)
    
    @classmethod
    def release(cls, node):
        """Release FoodNode back to pool."""
        if node is not None and len(cls._pool) < 2000:
            cls._pool.append(node)
    
    def init(self, x, y, amount, substrate_id=DEFAULT_SUBSTRATE_ID):
        """Initialize/reinitialize pool object."""
        self.x = x
        self.y = y
        self.amount = amount
        # Phase I: dietary substrate class (defaults to generic calorie for necromass/legacy food).
        self.substrate_id = substrate_id
        return self

class Bacterium:
    """
    Individual bacterium with position and species reference.
    """
    _pool = []
    
    def __init__(self, species, x, y, size=1.0, initial_food=None, initial_q=None):
        self.species = species
        self.x = x
        self.y = y
        self.base_size = size  # Unit: micrometers (um)
        # Volumetric scaling (size^3) assuming constant cellular density context.
        # Unit mappings: 1 distance unit = 1 um. 1 food unit = 1 pg glucose-equivalent organic nutrient.
        # Simple geometric scaling for 3D objects; ~1 pg wet mass per 1 um^3 bacterial volume; see DOCUMENTATION.md §4.1 (Nutrient Storage and Cell Size).
        self.food_capacity = getattr(species, 'food_capacity', DEFAULT_FOOD_CAPACITY_PG) * (self.base_size ** 3)
        self.food_level = initial_food if initial_food is not None else self.food_capacity
        
        # Allometric metabolic scaling: metabolic rate scales as Mass^(3/4). If Mass ~ size^3, then Meta ~ size^(2.25).
        # see DOCUMENTATION.md §4.2 (Metabolism and Uptake)
        # Units: pg glucose-equivalent nutrient consumed per second, later scaled by elapsed dt.
        self.metabolic_rate = getattr(species, 'metabolic_rate', DEFAULT_METABOLIC_RATE_PG_S) * (self.base_size ** 2.25)
        
        # Surface area scaling (size^2) governing nutrient absorption/consumption rates.
        # Square-cube law; surface area limits exchange with environment; see DOCUMENTATION.md §4.2 (Metabolism and Uptake).
        # Units: pg glucose-equivalent nutrient absorbed per second, later scaled by elapsed dt.
        self.consumption_speed = getattr(species, 'consumption_speed', DEFAULT_CONSUMPTION_SPEED_PG_S) * (self.base_size ** 2)
        
        # Dimensionless per-substrate capacity multipliers use SUBSTRATE_* order. The legacy node
        # path computes consumption_speed * multiplier, while the reduced model applies the same
        # trait to base Vmax; zero exactly disables uptake; see DOCUMENTATION.md section 10.4.
        self.substrate_affinity = (
            max(0.0, getattr(species, 'aff_calorie', DEFAULT_AFF_CALORIE)),
            max(0.0, getattr(species, 'aff_glycan', DEFAULT_AFF_GLYCAN)),
            max(0.0, getattr(species, 'aff_mdcf', DEFAULT_AFF_MDCF)),
        )
        
        # Kinematic scaling in low Reynolds number environments (microscale).
        # Viscous drag at low Reynolds number (Fd = 6*pi*eta*r*v); V scales inversely with radius (base_size); see DOCUMENTATION.md §4.5 (Rotational Diffusion).
        # Units: um/second.
        self.speed = getattr(species, 'speed', DEFAULT_BACTERIA_SPEED_UM_S) / max(0.1, self.base_size)
        
        cap = self.food_capacity if self.food_capacity > 0 else 1.0
        self.size = self.base_size * (BACTERIA_STARVED_RENDER_SCALE + BACTERIA_FULL_RENDER_SCALE_RANGE * (self.food_level / cap))
        
        self.age = 0.0
        self.lifespan_variance = max(LIFESPAN_VARIANCE_MIN, random.gauss(LIFESPAN_VARIANCE_MEAN, LIFESPAN_VARIANCE_STD))
        
        # Rotational kinematics in low Reynolds number regime.
        # Torque tau = 8 * pi * eta * r^3 * omega.
        # Assuming flagellar motor output torque scales with cell surface area (r^2), max angular velocity (omega_max) scales as r^2 / r^3 = 1 / r.
        # Base maximum turn rate is modeled as pi radians/second for a 1 um bacterium.
        self.turn_rate = math.pi / max(0.1, self.base_size)
        self.heading = random.uniform(0, 2 * math.pi) if 'math' in globals() else random.random()
        
        # Bacterial run-and-tumble temporal mechanics (E. coli chemotaxis model).
        # E. coli naturally "run" for ~1 second and "tumble" for ~0.1 seconds to reorient; see DOCUMENTATION.md §4.4 (Run-and-Tumble Motion).
        self.base_run_duration = self._sample_run_duration()
        self.base_tumble_duration = BACTERIA_TUMBLE_DURATION_SEC
        self.state_timer = self.base_run_duration
        self.is_tumbling = False
        self._tumble_reoriented = False
        self.last_metabolic_consumption_pg = 0.0
        # M2 cross-feeding bookkeeping (pg). last_food_uptake_pg drives overflow acetate secretion;
        # last_acetate_consumption_pg feeds acetate-oxidation O2/CO2 coupling in the env update.
        self.last_food_uptake_pg = 0.0
        self.last_acetate_consumption_pg = 0.0
        # Phase I diet analytics: per-step substrate mass (pg) ingested this tick, indexed by
        # SUBSTRATE_* id (0 = generic calorie, 1 = complex glycan, 2 = MDCF prebiotic). One slot
        # per substrate class in SUBSTRATE_ENERGY_YIELD so the Simulation can aggregate diet uptake.
        self.last_uptake_by_substrate = [0.0] * len(SUBSTRATE_ENERGY_YIELD)
        # M3 Baranyi physiological state q (see DOCUMENTATION.md §4.16 (Lag Phase)). Fresh inoculum starts at
        # the species' q0 (low -> long lag); daughter cells inherit the parent's adapted q so they do not re-lag.
        self.q = initial_q if initial_q is not None else getattr(species, 'baranyi_q0', DEFAULT_BARANYI_Q0)
        # M4 phenotypic switching: cells begin proliferative and may stochastically enter the dormant,
        # stress-tolerant persister state under starvation (see DOCUMENTATION.md §4.17 (Phenotypic Switching)).
        # Daughters are always born proliferative since persisters do not divide.
        self.is_persister = False
        # M5 phage predation (SIR host state) + HGT. A cell is born resistant with the species' initial
        # resistant fraction (heritable, chromosomal); infection and lysis are tracked per cell. The
        # metabolic burden factor (>=1) models the fitness cost of a horizontally acquired resistance
        # plasmid (plasmid carriage cost); it is 1.0 for cells with no carriage cost. see DOCUMENTATION.md §4.18 (Phage Predation and Horizontal Gene Transfer).
        self.phage_resistant = random.random() < getattr(species, 'phage_resistance', 0.0)
        self.infected = False
        self.lysis_timer = 0.0
        self.metabolic_burden_factor = 1.0

    def _sample_run_duration(self):
        """Return one E. coli-style run duration sampled in real elapsed seconds."""
        # Equation: t_run ~ Normal(mean=1.0 s, sd=0.2 s), clamped to avoid nonphysical negative runs.
        return max(BACTERIA_RUN_MIN_SEC, min(BACTERIA_RUN_MAX_SEC, random.gauss(BACTERIA_RUN_MEAN_SEC, BACTERIA_RUN_STD_SEC)))
    
    @classmethod
    def acquire(cls):
        """Acquire a Bacterium from pool or create new."""
        if cls._pool:
            return cls._pool.pop()
        return cls.__new__(cls)
    
    @classmethod
    def release(cls, bacterium):
        """Release Bacterium back to pool."""
        if bacterium is not None and len(cls._pool) < 5000:
            cls._pool.append(bacterium)
    
    def init(self, species, x, y, size=1.0, initial_food=None, initial_q=None):
        self.__init__(species, x, y, size, initial_food, initial_q)
        return self

    def step(self, food_nodes, bounds, dt, viewport_bounds=None, food_hash=None, temp_factor=1.0):
        # temp_factor in [0, 1] is the Ratkowsky temperature multiplier supplied by the species.
        # It scales metabolism, nutrient uptake, and motility so that cold/hot extremes slow the
        # cell toward dormancy (factor -> 0) while the species optimum leaves rates unchanged (-> 1).
        is_out_of_bounds = False
        if viewport_bounds is not None:
            is_out_of_bounds = (
                self.x < viewport_bounds.get('left') or 
                self.x > viewport_bounds.get('right') or 
                self.y < viewport_bounds.get('top') or 
                self.y > viewport_bounds.get('bottom')
            )

        if self.food_level < 0:
            return
            
        self.age += dt
        # Senescence: Replicative senescence and ROS (Reactive Oxygen Species) accumulation.
        # As age increases, cellular damage increases metabolic maintenance overhead.
        # Free-radical theory of aging: models linear increase in basal energy. see DOCUMENTATION.md §4.11 (Senescence).
        # Average time to double metabolic cost (senescence base lifespan) approx 3 minutes (180s) in sim.
        age_factor = 1.0 + (self.age / (SENESCENCE_BASE_LIFESPAN_SEC * self.lifespan_variance))
        
        # Metabolism: lose nutrient over elapsed real time.
        # Temperature (Ratkowsky) scales enzymatic turnover: cold/hot extremes depress basal demand,
        # which is why chilled cultures persist far longer than cultures at their growth optimum.
        # M4: dormant persisters burn only a fraction of basal maintenance (see DOCUMENTATION.md §4.17 (Phenotypic Switching)).
        persister_metabolic_factor = PERSISTER_METABOLIC_FACTOR if self.is_persister else 1.0
        # M5: a horizontally acquired resistance plasmid raises basal maintenance (carriage cost).
        metabolic_demand_pg = (self.metabolic_rate * dt) * age_factor * temp_factor * persister_metabolic_factor * self.metabolic_burden_factor
        self.last_metabolic_consumption_pg = min(max(0.0, self.food_level), metabolic_demand_pg)
        self.food_level = max(0, self.food_level - metabolic_demand_pg)

        # M3 Baranyi lag: advance the physiological state q. dq/dt = mu_max * q, with mu_max taken as
        # the species' temperature-scaled specific growth rate. Cold/hot extremes (temp_factor -> 0)
        # slow adaptation and lengthen the lag, matching observed temperature-dependent lag times.
        # Solved analytically over the tick (q_{t+dt} = q_t * e^{mu_max*dt}) so the result is exact and
        # frame-rate independent rather than a forward-Euler approximation.
        mu_max = getattr(self.species, 'base_rep', 0.0)
        if mu_max > 0.0 and self.q < BARANYI_Q_STATE_MAX:
            self.q = min(BARANYI_Q_STATE_MAX, self.q * math.exp(mu_max * temp_factor * dt))
        
        # Update size based on fullness
        cap = self.food_capacity if self.food_capacity > 0 else 1.0
        fullness = max(0.0, min(1.0, self.food_level / cap))
        self.size = self.base_size * (BACTERIA_STARVED_RENDER_SCALE + BACTERIA_FULL_RENDER_SCALE_RANGE * fullness)

        # M4 phenotypic switching (see DOCUMENTATION.md §4.17 (Phenotypic Switching)). Reversible, stochastic
        # transition between the proliferative and dormant persister states, gated by internal fullness as a
        # faithful proxy for sustained LOCAL nutrient balance. Probabilities use the exponential CDF
        # 1 - e^{-rate*dt} so the per-tick chance is exact and frame-rate independent. Disabled when the
        # species' switching propensity is 0, reproducing pre-M4 behaviour exactly.
        persister_switch = getattr(self.species, 'persister_switch', 0.0)
        if persister_switch > 0.0:
            s1 = getattr(self.species, 'persister_enter_thresh', DEFAULT_PERSISTER_ENTER_THRESH)
            s2 = getattr(self.species, 'persister_exit_thresh', DEFAULT_PERSISTER_EXIT_THRESH)
            if not self.is_persister:
                # Baseline bet-hedging plus a starvation-induced term that strengthens as reserves empty.
                enter_rate = PERSISTER_SPONTANEOUS_SWITCH_RATE_PER_SEC
                if fullness < s1:
                    enter_rate += PERSISTER_STRESS_SWITCH_RATE_PER_SEC * (1.0 - fullness / max(s1, NUMERIC_DISTANCE_EPSILON_UM))
                if random.random() < 1.0 - math.exp(-enter_rate * persister_switch * dt):
                    self.is_persister = True
            else:
                # Spontaneous reversion plus a resuscitation term that strengthens as reserves refill.
                wake_rate = PERSISTER_SPONTANEOUS_WAKE_RATE_PER_SEC
                if fullness > s2:
                    wake_rate += PERSISTER_WAKE_RATE_PER_SEC * ((fullness - s2) / max(1.0 - s2, NUMERIC_DISTANCE_EPSILON_UM))
                if random.random() < 1.0 - math.exp(-wake_rate * dt):
                    self.is_persister = False
        
        # Compute effective greediness that rises toward 1 as stomach empties
        base_gx = getattr(self.species, 'greediness', 0.5)
        hunger = 1.0 - fullness
        gx_eff = max(0.0, min(1.0, base_gx + (1.0 - base_gx) * hunger))
        
        # Pathfind with chemotaxis: score food nodes by physical diffusion gradients.
        # Fickian nutrient fields and Weber-Fechner receptor compression; see DOCUMENTATION.md §4.3 (Food Sensing and Chemotaxis).
        target = None
        best_score = -float('inf')
        closest_node = None
        min_dist = float('inf')
        gradient_x = 0.0
        gradient_y = 0.0
        max_range = self.species.sensor_radius
        
        candidate_nodes = food_hash.query(self.x, self.y, max_range) if food_hash else food_nodes
        if is_out_of_bounds:
            candidate_nodes = [] # Hybrid ABM/ODE Fast-forward. ODE translation applies instead.
            
        for node in candidate_nodes:
            if node.amount <= MIN_FOOD_NODE_AMOUNT_PG:
                continue
            # Phase I diet model: ignore substrate classes this cell cannot metabolise (affinity 0)
            # so chemotaxis and targeting are selective — MDCF prebiotics attract only competent taxa.
            node_aff = self.substrate_affinity[getattr(node, 'substrate_id', DEFAULT_SUBSTRATE_ID)]
            if node_aff <= 0.0:
                continue
            dx = node.x - self.x
            dy = node.y - self.y
            dist = math.hypot(dx, dy)
            if dist < min_dist:
                min_dist = dist
                closest_node = node
            if dist > max_range:
                continue
            
            # Equation: C = amount / (pi * (r^2 + a^2)) models a softened 2D nutrient patch without a zero-distance singularity.
            concentration = node.amount / (math.pi * ((dist ** 2) + (CHEMOTAXIS_CORE_RADIUS_UM ** 2)))
            # Equation: score = greed * affinity * log1p(C / C_ref); affinity weights attraction toward usable substrates.
            score = gx_eff * node_aff * math.log1p(concentration / CHEMOTAXIS_REFERENCE_PG_PER_UM2)
            if dist > NUMERIC_DISTANCE_EPSILON_UM:
                gradient_x += (dx / dist) * score
                gradient_y += (dy / dist) * score
            
            if score > best_score:
                best_score = score
                target = node
        current_score = best_score if target is not None else 0.0
        has_gradient = math.hypot(gradient_x, gradient_y) > CHEMOTAXIS_GRADIENT_EPSILON
                
        # State management: Run vs Tumble
        self.last_score = getattr(self, 'last_score', 0.0)
        self.state_timer -= dt
        if self.state_timer <= 0.0:
            self.is_tumbling = not self.is_tumbling
            if self.is_tumbling:
                self.state_timer = self.base_tumble_duration
                self._tumble_reoriented = False
            else:
                if current_score > self.last_score:
                    # Gradient is improving - extend run duration by factor of 2.5
                    self.state_timer = self._sample_run_duration() * BACTERIA_CHEMOTAXIS_RUN_EXTENSION_FACTOR
                else:
                    self.state_timer = self._sample_run_duration()
        
        self.last_score = current_score

        # Base movement speed; flagellar motor output is enzyme-driven, so it scales with temperature.
        spd = self.speed * dt * temp_factor

        if target is not None and has_gradient and not self.is_tumbling:
            desired_heading = math.atan2(gradient_y, gradient_x)
            # Equation: signed angular error uses atan2(sin, cos) to wrap to [-pi, pi] without discontinuity.
            heading_error = math.atan2(math.sin(desired_heading - self.heading), math.cos(desired_heading - self.heading))
            # Equation: dtheta <= omega * dt gives smooth low-Reynolds reorientation while translation remains forward-only.
            max_turn = CHEMOTAXIS_STEER_RATE_RAD_S * gx_eff * dt
            self.heading += max(-max_turn, min(max_turn, heading_error))

        # Multi-node consumption: feed from all candidate nodes within eat distance.
        # Membrane transporter kinetics are temperature-limited, so uptake scales with temp_factor.
        # max_eat_amt is the shared transporter throughput budget (mass, pg) available this tick.
        max_eat_amt = self.consumption_speed * dt * temp_factor
        amt_eaten_total = 0.0
        # Reset this tick's per-substrate uptake ledger before the consumption loop so each frame
        # reports only the mass ingested this step (the Simulation accumulates the running totals).
        for _si in range(len(self.last_uptake_by_substrate)):
            self.last_uptake_by_substrate[_si] = 0.0
        is_grazing = False
        
        for node in candidate_nodes:
            if node.amount <= MIN_FOOD_NODE_AMOUNT_PG:
                continue
            # Phase I diet model: per-substrate affinity gates and scales uptake; substrates with
            # zero affinity cannot be metabolised and are skipped entirely.
            sid = getattr(node, 'substrate_id', DEFAULT_SUBSTRATE_ID)
            node_aff = self.substrate_affinity[sid]
            if node_aff <= 0.0:
                continue
            # Geometric volumetric scaling: V = 4/3 * pi * r^3 -> r = cbrt(3V / 4pi)
            target_radius = max(FOOD_NODE_MIN_RADIUS_UM, ((3.0 * max(FOOD_NODE_MIN_VOLUME_PG, node.amount)) / (4.0 * math.pi)) ** (1.0 / 3.0))
            eat_dist = target_radius + self.size + FOOD_EATING_DISTANCE_MARGIN_UM
            
            d = math.hypot(node.x - self.x, node.y - self.y)
            
            if d <= eat_dist:
                is_grazing = True
                
            # If close enough to feed (small margin)
            remaining_budget = max_eat_amt - amt_eaten_total
            if d <= eat_dist + FOOD_EATING_DISTANCE_GRACE_UM and remaining_budget > 0.0:
                # Energy yield (glucose-equivalent fraction) gained per pg of this substrate ingested.
                yield_f = SUBSTRATE_ENERGY_YIELD.get(sid, 1.0)
                # Internal storage headroom expressed back in ingestible substrate mass (pg).
                room_mass = max(0.0, self.food_capacity - self.food_level) / max(yield_f, MIN_FOOD_NODE_AMOUNT_PG)
                # Affinity scales the per-substrate uptake rate; the draw can never exceed the
                # remaining shared transporter budget, the node's mass, or storage headroom.
                eat_amt = min(remaining_budget * node_aff, remaining_budget, node.amount, room_mass)
                if eat_amt > 0.0:
                    # Mass is conserved (node loses eat_amt); stored energy gains eat_amt * yield.
                    self.food_level += eat_amt * yield_f
                    node.amount -= eat_amt
                    amt_eaten_total += eat_amt
                    # Attribute the ingested mass to its substrate class for diet analytics.
                    if 0 <= sid < len(self.last_uptake_by_substrate):
                        self.last_uptake_by_substrate[sid] += eat_amt

        # Record this tick's glucose-equivalent uptake mass so the species step can excrete overflow acetate (M2).
        self.last_food_uptake_pg = amt_eaten_total

        if is_grazing:
            spd *= BACTERIA_GRAZING_SPEED_FACTOR # Grazing increases drag/friction but does not create a full stop.
                
        if self.is_tumbling:
            # During tumble, motility is spent almost entirely on reorientation.
            if not getattr(self, '_tumble_reoriented', False):
                if has_gradient:
                    desired_heading = math.atan2(gradient_y, gradient_x)
                    # Equation: theta_new = theta_gradient + Normal(0, sigma), a local-gradient approximation of biased post-tumble runs.
                    self.heading = desired_heading + random.gauss(0.0, BACTERIA_TUMBLE_TURN_STD_RAD)
                else:
                    # Equation: |delta theta| ~ Normal(68 deg, sigma), with random sign for unbiased tumbles.
                    turn_mag = max(0.0, random.gauss(BACTERIA_TUMBLE_MEAN_TURN_RAD, BACTERIA_TUMBLE_TURN_STD_RAD))
                    self.heading += random.choice((-1.0, 1.0)) * turn_mag
                self._tumble_reoriented = True
            spd *= BACTERIA_TUMBLE_TRANSLATION_FACTOR # Very slow forward translation during tumble
        
        # Apply Brownian rotational noise (rotational diffusion).
        # Rotational diffusion coefficient scales roughly ~ 1 / r^3; see DOCUMENTATION.md §4.5 (Rotational Diffusion).
        # D_rot = kT / (8 pi eta r^3). A standard value for 1 um sphere in water at room temp is ~0.16 rad^2/s.
        rotational_diffusion = BACTERIA_ROTATIONAL_DIFFUSION_COEFFICIENT_RAD2_S / max(0.1, self.base_size ** 3)
        # Brownian angular variance is sigma^2 = 2 D_rot dt, so the sampled angle uses sqrt(2 D_rot dt).
        rotational_sigma_rad = math.sqrt(BROWNIAN_ROTATIONAL_VARIANCE_FACTOR * rotational_diffusion * dt)
        self.heading += random.gauss(0, rotational_sigma_rad)
        
        # Spatial heterogeneity: cells in the outer mucus annulus adhere to the matrix and slow
        # down, forming an anchoring refuge from luminal washout (see DOCUMENTATION.md §4.20 (Spatial
        # Heterogeneity: Lumen vs Mucus Layer)). `bounds` carries the petri centre/radius also used
        # for the circular boundary reflection below.
        cx_b, cy_b, r_b = bounds
        if r_b > 0.0 and math.hypot(self.x - cx_b, self.y - cy_b) >= r_b * MUCUS_LAYER_RADIUS_FRACTION:
            spd *= MUCUS_MOTILITY_FACTOR
        
        # Forward translation
        self.x += math.cos(self.heading) * spd
        self.y += math.sin(self.heading) * spd
            
        # Circular petri dish boundary constraints
        cx, cy, r_bound = bounds
        dx = self.x - cx
        dy = self.y - cy
        dist = math.hypot(dx, dy)
        
        if dist > r_bound and dist > NUMERIC_DISTANCE_EPSILON_UM:
            nx = dx / dist
            ny = dy / dist
            
            # Reposition inside boundary slightly to prevent boundary sticking
            self.x = cx + nx * (r_bound - BACTERIA_BOUNDARY_OFFSET_UM)
            self.y = cy + ny * (r_bound - BACTERIA_BOUNDARY_OFFSET_UM)
            
            # Proper geometric reflection logic off circle normal
            # velocity vector
            vx = math.cos(self.heading)
            vy = math.sin(self.heading)
            
            # v' = v - 2(v.n)n
            dot = vx * nx + vy * ny
            
            # Mathematical Rationale: Only apply specular reflection if the velocity vector
            # is pointing out of the circle bounds (dot product > 0). If it is already pointing
            # inwards, the particle shouldn't be reflected back out.
            if dot > 0:
                vx_new = vx - 2 * dot * nx
                vy_new = vy - 2 * dot * ny
                self.heading = math.atan2(vy_new, vx_new) + random.uniform(-BACTERIA_BOUNDARY_REFLECTION_NOISE_RAD, BACTERIA_BOUNDARY_REFLECTION_NOISE_RAD)

    def step_fba_kinematics(
        self,
        ecm,
        bounds,
        dt_seconds,
        physiology,
        viewport_bounds=None,
        temp_factor=1.0,
    ):
        """Move from continuous substrate gradients without legacy pg metabolism or uptake."""
        if not isinstance(physiology, FbaPhysiology):
            raise TypeError('physiology must be FbaPhysiology')
        if dt_seconds < 0.0 or not math.isfinite(dt_seconds):
            raise ValueError('dt_seconds must be finite and nonnegative')

        self.age += dt_seconds
        self.size = self.base_size * physiology.render_scale
        coordinate = ecm.world_to_cell(self.x, self.y)
        cx, cy = coordinate
        neighbor_coordinates = (
            coordinate,
            (cx - 1, cy),
            (cx + 1, cy),
            (cx, cy - 1),
            (cx, cy + 1),
        )
        concentrations = ecm.sample_substrate_batch_m(neighbor_coordinates)
        scores = tuple(
            sum(float(row[index]) * self.substrate_affinity[index] for index in range(3))
            for row in concentrations
        )
        current_score = scores[0] * physiology.chemotactic_drive
        gradient_x = (scores[2] - scores[1]) * physiology.chemotactic_drive
        gradient_y = (scores[4] - scores[3]) * physiology.chemotactic_drive
        has_gradient = math.hypot(gradient_x, gradient_y) > CHEMOTAXIS_GRADIENT_EPSILON

        self.last_score = getattr(self, 'last_score', 0.0)
        self.state_timer -= dt_seconds
        if self.state_timer <= 0.0:
            self.is_tumbling = not self.is_tumbling
            if self.is_tumbling:
                self.state_timer = self.base_tumble_duration
                self._tumble_reoriented = False
            elif current_score > self.last_score:
                self.state_timer = (
                    self._sample_run_duration() * BACTERIA_CHEMOTAXIS_RUN_EXTENSION_FACTOR
                )
            else:
                self.state_timer = self._sample_run_duration()
        self.last_score = current_score

        if has_gradient and not self.is_tumbling:
            desired_heading = math.atan2(gradient_y, gradient_x)
            heading_error = math.atan2(
                math.sin(desired_heading - self.heading),
                math.cos(desired_heading - self.heading),
            )
            max_turn = (
                CHEMOTAXIS_STEER_RATE_RAD_S
                * physiology.chemotactic_drive
                * dt_seconds
            )
            self.heading += max(-max_turn, min(max_turn, heading_error))

        speed = self.speed * dt_seconds * temp_factor * physiology.motility_scale
        if self.is_tumbling:
            if not getattr(self, '_tumble_reoriented', False):
                if has_gradient:
                    self.heading = math.atan2(gradient_y, gradient_x) + random.gauss(
                        0.0,
                        BACTERIA_TUMBLE_TURN_STD_RAD,
                    )
                else:
                    turn_magnitude = max(
                        0.0,
                        random.gauss(
                            BACTERIA_TUMBLE_MEAN_TURN_RAD,
                            BACTERIA_TUMBLE_TURN_STD_RAD,
                        ),
                    )
                    self.heading += random.choice((-1.0, 1.0)) * turn_magnitude
                self._tumble_reoriented = True
            speed *= BACTERIA_TUMBLE_TRANSLATION_FACTOR

        rotational_diffusion = (
            BACTERIA_ROTATIONAL_DIFFUSION_COEFFICIENT_RAD2_S
            / max(0.1, self.base_size ** 3)
        )
        rotational_sigma_rad = math.sqrt(
            BROWNIAN_ROTATIONAL_VARIANCE_FACTOR * rotational_diffusion * dt_seconds
        )
        self.heading += random.gauss(0.0, rotational_sigma_rad)

        center_x, center_y, radius = bounds
        if radius > 0.0 and math.hypot(
            self.x - center_x,
            self.y - center_y,
        ) >= radius * MUCUS_LAYER_RADIUS_FRACTION:
            speed *= MUCUS_MOTILITY_FACTOR
        self.x += math.cos(self.heading) * speed
        self.y += math.sin(self.heading) * speed

        delta_x = self.x - center_x
        delta_y = self.y - center_y
        distance = math.hypot(delta_x, delta_y)
        if distance > radius and distance > NUMERIC_DISTANCE_EPSILON_UM:
            normal_x = delta_x / distance
            normal_y = delta_y / distance
            self.x = center_x + normal_x * (radius - BACTERIA_BOUNDARY_OFFSET_UM)
            self.y = center_y + normal_y * (radius - BACTERIA_BOUNDARY_OFFSET_UM)
            velocity_x = math.cos(self.heading)
            velocity_y = math.sin(self.heading)
            outward_component = velocity_x * normal_x + velocity_y * normal_y
            if outward_component > 0.0:
                reflected_x = velocity_x - 2.0 * outward_component * normal_x
                reflected_y = velocity_y - 2.0 * outward_component * normal_y
                self.heading = math.atan2(reflected_y, reflected_x) + random.uniform(
                    -BACTERIA_BOUNDARY_REFLECTION_NOISE_RAD,
                    BACTERIA_BOUNDARY_REFLECTION_NOISE_RAD,
                )

    def attempt_fission(self, envf, rep_prob, rep_thresh, dt):
        """
        Evaluate and execute cell division (binary fission).
        ATP/biomass synthesis cost: chromosome, ribosome, and cell-wall assembly equate to ~40 pg
        glucose-equivalent nutrient per 1 um^3 of new cell volume. Requires sufficient accumulated biomass.
        see DOCUMENTATION.md §4.7 (Reproduction).
        """
        if rep_prob > 0 and random.random() < rep_prob and self.food_level > self.food_capacity * rep_thresh:
            child_volume = self.base_size ** 3
            atp_cost = FISSION_COST_PG_PER_UM3 * child_volume
            
            if self.food_level > atp_cost:
                remaining_food = self.food_level - atp_cost
                parent_food = remaining_food * 0.5
                child_food = remaining_food * 0.5
                
                # Ensure both cells have enough nutrient to survive the next elapsed-time update.
                min_survival_food = self.metabolic_rate * dt
                if parent_food > min_survival_food and child_food > min_survival_food:
                    self.food_level = parent_food
                    # Daughter inherits the parent's adapted physiological state q (no re-lag): an
                    # already-dividing cell's progeny are pre-adapted (see DOCUMENTATION.md §4.16 (Lag Phase)).
                    child = Bacterium(self.species, self.x, self.y, size=max(0.1, random.gauss(self.base_size, self.base_size * 0.2)), initial_food=child_food, initial_q=self.q)
                    # M5: phage resistance is heritable. A susceptible lineage gains resistance by rare
                    # spontaneous mutation per division (see DOCUMENTATION.md §4.18 (Phage Predation and
                    # Horizontal Gene Transfer)), the seed of the post-crash rebound. The plasmid carriage
                    # cost is inherited with the trait.
                    child.phage_resistant = self.phage_resistant
                    if not child.phage_resistant and random.random() < PHAGE_RESISTANCE_MUTATION_RATE:
                        child.phage_resistant = True
                    child.metabolic_burden_factor = self.metabolic_burden_factor
                    return child
        return None

class BacteriaSpecies:
    """
    Represents one species, tracking its parameters,
    individual members, and population history.
    """
    def __init__(self, name, init, rep, aging, so, sf, speed, ph, color,
                 food_capacity=DEFAULT_FOOD_CAPACITY_PG, metabolic_rate=DEFAULT_METABOLIC_RATE_PG_S,
                 consumption_speed=DEFAULT_CONSUMPTION_SPEED_PG_S, greediness=DEFAULT_GREEDINESS,
                 t6ss_activity=DEFAULT_T6SS_ACTIVITY, bacteriocin_production=DEFAULT_BACTERIOCIN_PRODUCTION,
                 immunity_factor=DEFAULT_IMMUNITY_FACTOR, rep_thresh=DEFAULT_REPRODUCTION_THRESHOLD,
                 share_thresh=DEFAULT_SHARE_THRESHOLD, predation_thresh=DEFAULT_PREDATION_THRESHOLD,
                 base_size=DEFAULT_BASE_SIZE_UM, sensor_radius=DEFAULT_SENSOR_RADIUS_UM,
                 temp_min=DEFAULT_TEMP_MIN_C, temp_max=DEFAULT_TEMP_MAX_C,
                 acetate_secretion=DEFAULT_ACETATE_SECRETION, acetate_uptake=DEFAULT_ACETATE_UPTAKE,
                 baranyi_q0=DEFAULT_BARANYI_Q0,
                 persister_switch=DEFAULT_PERSISTER_SWITCH,
                 persister_enter_thresh=DEFAULT_PERSISTER_ENTER_THRESH,
                 persister_exit_thresh=DEFAULT_PERSISTER_EXIT_THRESH,
                 phage_resistance=DEFAULT_PHAGE_RESISTANCE,
                 hgt_competence=DEFAULT_HGT_COMPETENCE,
                 aff_calorie=DEFAULT_AFF_CALORIE, aff_glycan=DEFAULT_AFF_GLYCAN,
                 aff_mdcf=DEFAULT_AFF_MDCF):
        # Store parameters
        self.name = name
        self.base_rep = rep       # base reproduction probability
        self.base_aging = aging   # base aging probability
        self.o2_sens = so  # -1 to 1
        self.food_sens = sf
        self.speed = speed  # pixels/sec
        self.preferred_ph = ph
        self.color = color        # draw color for this species
        self.food_capacity = food_capacity
        self.metabolic_rate = metabolic_rate
        self.consumption_speed = consumption_speed
        # New attributes
        self.greediness = max(0.0, min(float(greediness), 1.0))
        self.t6ss_activity = max(0.0, min(float(t6ss_activity), 1.0))
        self.bacteriocin_production = max(0.0, min(float(bacteriocin_production), 1.0))
        self.immunity_factor = max(0.0, min(float(immunity_factor), 1.0))
        self.rep_thresh = max(0.0, min(float(rep_thresh), 1.0))
        self.share_thresh = max(0.0, min(float(share_thresh), 1.0))
        self.predation_thresh = max(0.0, min(float(predation_thresh), 1.0))
        self.base_size = max(0.1, float(base_size))
        self.sensor_radius = max(0.0, float(sensor_radius))
        # Cardinal growth temperatures (Ratkowsky square-root model). T_max is kept strictly above
        # T_min so the curve is well defined; the optimum is located numerically once at construction.
        self.temp_min = float(temp_min)
        self.temp_max = max(self.temp_min + 1.0, float(temp_max))
        # Precompute the peak of sqrt(mu) over [T_min, T_max] so temperature_factor returns a value
        # in [0, 1] that equals 1.0 at this species' own thermal optimum.
        self._temp_sqrt_mu_peak = self._compute_temp_sqrt_mu_peak()
        # M2 cross-feeding traits, both in [0, 1]. acetate_secretion is the fraction of fresh glucose
        # uptake excreted as overflow acetate; acetate_uptake scales the cell's acetate-oxidation rate.
        self.acetate_secretion = max(0.0, min(float(acetate_secretion), 1.0))
        self.acetate_uptake = max(0.0, min(float(acetate_uptake), 1.0))
        # M3 Baranyi lag: initial physiological state q0 of fresh inoculum. Small q0 -> long lag
        # (h0 = ln(1 + 1/q0)); clamped to a finite, non-zero range so adaptation always completes.
        self.baranyi_q0 = max(BARANYI_Q0_MIN, min(float(baranyi_q0), BARANYI_Q0_MAX))
        # M4 phenotypic switching traits. persister_switch in [0,1] is the switching propensity (0
        # disables persistence). The two thresholds are internal-fullness levels: enter (S1) is the
        # starvation trigger, exit (S2) is the recovery trigger; S2 is kept strictly above S1 so the
        # enter/exit bands do not overlap and chatter.
        self.persister_switch = max(0.0, min(float(persister_switch), 1.0))
        self.persister_enter_thresh = max(0.0, min(float(persister_enter_thresh), 1.0))
        self.persister_exit_thresh = max(self.persister_enter_thresh + PERSISTER_THRESHOLD_MIN_GAP,
                                         min(float(persister_exit_thresh), 1.0))
        # M5 phage/HGT traits, both in [0, 1]. phage_resistance is the heritable initial fraction of
        # cells born resistant to adsorption; hgt_competence gates contact-dependent horizontal
        # transfer of resistance (0 disables both donating and receiving).
        self.phage_resistance = max(0.0, min(float(phage_resistance), 1.0))
        self.hgt_competence = max(0.0, min(float(hgt_competence), 1.0))
        # Dimensionless per-substrate uptake-capacity multipliers are nonnegative. The legacy path
        # scales consumption_speed and the reduced model scales substrate-specific base Vmax;
        # zero disables uptake exactly; see DOCUMENTATION.md section 10.4.
        self.aff_calorie = max(0.0, float(aff_calorie))
        self.aff_glycan = max(0.0, float(aff_glycan))
        self.aff_mdcf = max(0.0, float(aff_mdcf))
        # Requested founder-population size. Members are NOT placed here: the inoculum is positioned
        # by the owning Simulation (see Simulation._seed_initial_population) once the petri-dish
        # geometry is known, so founder placement is dish-aware, lumen-biased, and seed-reproducible
        # instead of the legacy fixed top-left rectangle.
        self.init_count = max(0, int(init))
        self.members = []
        self.history = []         # track population over time
        self.total_births = 0     # overall births
        self.total_deaths = 0     # overall deaths
        # Phase I reproducibility: elapsed sim-seconds at which this species first reached zero
        # population (irreversible from 0 here). None until/unless the species goes extinct.
        self.extinction_time = None

    def _ratkowsky_sqrt_mu(self, temp_c):
        """Raw square-root growth term sqrt(mu) of the Ratkowsky model (unnormalized).

        Equation: sqrt(mu) = b*(T - T_min)*(1 - exp(c*(T - T_max))).
        Returns 0.0 outside the (T_min, T_max) biokinetic window where no growth occurs.
        see DOCUMENTATION.md §4.14 (Temperature Kinetics).
        """
        if temp_c <= self.temp_min or temp_c >= self.temp_max:
            return 0.0
        root = (RATKOWSKY_B_SQRT_PER_C * (temp_c - self.temp_min)
                * (1.0 - math.exp(RATKOWSKY_C_PER_C * (temp_c - self.temp_max))))
        return root if root > 0.0 else 0.0

    def _compute_temp_sqrt_mu_peak(self):
        """Locate the maximum of sqrt(mu) across [T_min, T_max] by uniform sampling (one-time)."""
        peak = 0.0
        n = RATKOWSKY_NORMALIZATION_SAMPLES
        span = self.temp_max - self.temp_min
        for i in range(n + 1):
            t = self.temp_min + span * (i / n)
            peak = max(peak, self._ratkowsky_sqrt_mu(t))
        return peak if peak > 0.0 else 1.0

    def temperature_factor(self, temp_c):
        """Ratkowsky growth-rate multiplier in [0, 1], normalized to this species' optimum.

        mu(T) ~ sqrt(mu)^2, so the multiplier is (sqrt(mu)/sqrt(mu)_peak)^2. It is 1.0 at the
        species optimum and 0.0 at/beyond its cardinal T_min and T_max temperatures.
        """
        root = self._ratkowsky_sqrt_mu(temp_c)
        if root <= 0.0:
            return 0.0
        ratio = root / self._temp_sqrt_mu_peak
        return max(0.0, min(1.0, ratio * ratio))

    def step(self, env, cap, w, h, food_nodes, dt, viewport_bounds=None, food_hash=None):
        """
        Process one simulation step for this species:
        - Move each bacterium
        - Decide births and deaths based on rates,
          environment, and carrying capacity
        Returns (birth_count, death_count).
        """
        births, deaths = [], []
        # Ensure species_list is set
        if not hasattr(self, 'species_list') or not self.species_list:
            self.species_list = [self]
        total = sum(len(sp.members) for sp in self.species_list)
        safe_cap = cap if cap > 0 else 1
        capf = max(0, 1 - total / safe_cap) if safe_cap else 0

        # pH suitability: Gaussian curve centered on preferred pH
        curr_ph = env.get('ph', NEUTRAL_PH)
        ph_diff = curr_ph - self.preferred_ph if self.preferred_ph is not None else 0
        ph_tolerance = PH_TOLERANCE_GAUSSIAN_SD
        ph_suit = math.exp(-0.5 * (ph_diff / ph_tolerance)**2)

        # O2/CO2 suitability
        o2_level = env.get('o2', ATMOSPHERIC_O2_FRACTION)
        if self.o2_sens > 0:
            # Aerobes: Need O2 to survive. Scale up to the documented saturation threshold.
            o2_suit = min(1.0, o2_level / AEROBE_O2_SATURATION_LEVEL) * self.o2_sens + (1 - self.o2_sens)
        elif self.o2_sens < 0:
            # Anaerobes: Poisoned by O2. Toxicity kicks in heavily over the documented threshold.
            toxicity = min(1.0, o2_level / ANAEROBE_O2_TOXICITY_THRESHOLD) * abs(self.o2_sens)
            o2_suit = 1.0 - toxicity
        else:
            o2_suit = 1.0  # Indifferent

        # Overall environmental fitness [0, 1]
        envf = ph_suit * max(0.0, o2_suit)

        # Temperature kinetics (Ratkowsky): a per-species multiplier in [0, 1] driving growth,
        # metabolism, uptake, and motility. Kept separate from envf so it can scale per-cell rates.
        curr_temp = env.get('temperature', DEFAULT_ENV_TEMPERATURE_C)
        temp_factor = self.temperature_factor(curr_temp)

        ds_cx, ds_cy = getattr(self.sim, 'petri_center', (w/2, h/2))[0] if hasattr(self, 'sim') else w/2, getattr(self.sim, 'petri_center', (w/2, h/2))[1] if hasattr(self, 'sim') else h/2
        ds_r = getattr(self.sim, 'petri_radius', math.hypot(w, h)/2 + PETRI_RADIUS_PADDING_UM) if hasattr(self, 'sim') else math.hypot(w, h)/2 + PETRI_RADIUS_PADDING_UM
        # M5: free-phage field for the lytic cycle. Adsorption, latent countdown, and lysis are all
        # resolved here per cell so a single tick handles the full SIR transition (S->I->lysed).
        ecm = getattr(getattr(self, 'sim', None), 'ecm', None)

        for b in list(self.members):
            b.step(food_nodes, bounds=(ds_cx, ds_cy, ds_r), dt=dt, viewport_bounds=viewport_bounds, food_hash=food_hash, temp_factor=temp_factor)

            # --- M5 Phage predation (lytic SIR cycle) ---
            phage_lysed = False
            if ecm is not None:
                if b.infected:
                    # Latent (eclipse) period countdown; on expiry the host lyses (see DOCUMENTATION.md §4.18 (Phage Predation and Horizontal Gene Transfer)).
                    b.lysis_timer -= dt
                    if b.lysis_timer <= 0.0:
                        deaths.append(b)
                        phage_lysed = True
                        # Burst release: progeny virions deposited into the local phage field.
                        ecm.add(b.x, b.y, DEFAULT_PHAGE_BURST_SIZE * PHAGE_ECM_UNITS_PER_VIRION, layer='phage')
                        # Lysed-cell debris returns as labile necromass (no cytosol recovery: hijacked by phage).
                        cell_volume_um3 = b.base_size ** 3
                        necromass_pg = cell_volume_um3 * CELL_WET_MASS_PG_PER_UM3 * CELL_DRY_MASS_FRACTION * LABILE_NECROMASS_FRACTION
                        if necromass_pg > MIN_FOOD_NODE_AMOUNT_PG:
                            food_nodes.append(FoodNode.acquire().init(b.x, b.y, necromass_pg))
                elif not b.phage_resistant:
                    # Mass-action adsorption: P(adsorb) = 1 - exp(-k * [phage] * dt); see DOCUMENTATION.md §4.18 (Phage Predation and Horizontal Gene Transfer).
                    local_phage = ecm.sample(b.x, b.y, layer='phage') / PHAGE_ECM_UNITS_PER_VIRION
                    if local_phage > PHAGE_MIN_ADSORBABLE and \
                            random.random() < 1.0 - math.exp(-DEFAULT_PHAGE_ADSORPTION_RATE * local_phage * dt):
                        b.infected = True
                        b.lysis_timer = DEFAULT_PHAGE_LATENT_PERIOD_SEC
                        # One free virion is removed from the field upon successful adsorption.
                        ecm.add(b.x, b.y, -PHAGE_ECM_UNITS_PER_VIRION, layer='phage')
            if phage_lysed:
                # The cell is destroyed this tick; skip its reproduction and senescence trials.
                continue

            # Reproduction trial based on the logistic growth model modified by environmental fitness.
            # dp/dt = r*p*(1 - p/K); `capf` represents the limiting (1 - p/K) term. see DOCUMENTATION.md §4.7 (Reproduction).
            # The Ratkowsky temp_factor multiplies the intrinsic rate r, gating growth to the biokinetic range.
            # The Baranyi adjustment function alpha(q) = q/(1+q) gates growth during lag: ~0 while the cell
            # is still adapting (q small) and -> 1 once adapted, reproducing the flat-then-exponential curve.
            baranyi_alpha = b.q / (1.0 + b.q)
            rep_prob = (self.base_rep * dt) * capf * envf * temp_factor * baranyi_alpha
            rep_prob = max(0, min(rep_prob, 1))  # Clamp to [0,1]
            
            # M4: dormant persisters halt division entirely until they resuscitate (see DOCUMENTATION.md §4.17 (Phenotypic Switching)).
            # M5: phage-hijacked (infected) cells cannot divide while the virus commandeers replication.
            if not b.is_persister and not b.infected:
                new_bacterium = b.attempt_fission(envf, rep_prob, self.rep_thresh, dt)
                if new_bacterium:
                    births.append(new_bacterium)
            
            # Senescence modeled via Gompertz law of mortality: mu(x) = A * e^(alpha * x) + B.
            # Exponential increase in death rate due to intrinsic cellular senescence; alpha governs the
            # steepness of the decay curve. see DOCUMENTATION.md §4.11 (Senescence).
            gompertz_alpha = GOMPERTZ_MORTALITY_ALPHA / max(MIN_LIFESPAN_VARIANCE_DENOMINATOR, b.lifespan_variance)
            age_death_penalty = (self.base_aging * dt) * math.exp(gompertz_alpha * b.age)
            # Cap the penalty to avoid overflow probabilities per simulation tick
            age_death_penalty = min(1.0, age_death_penalty)
            
            # Death trial: extrinsic (environmental toxicity) + intrinsic (Gompertz senescence) + starvation.
            # M4: persisters are stress-tolerant, so both the environmental and senescence hazards are
            # attenuated by the persister resistance factor (see DOCUMENTATION.md §4.17 (Phenotypic Switching)).
            persister_resistance = PERSISTER_RESISTANCE_FACTOR if b.is_persister else 1.0
            toxic_death_prob = (((1.0 - envf) * ENV_TOXIC_DEATH_RATE_PER_SEC) * dt + age_death_penalty) * persister_resistance
            if b.food_level <= 0 or random.random() < toxic_death_prob:
                deaths.append(b)
                if b.food_level >= 0:
                    # Decomposition equation: reusable = cytosolic leak + labile dry necromass, measured as glucose-equivalent nutrient.
                    cell_volume_um3 = b.base_size ** 3
                    necromass_pg = cell_volume_um3 * CELL_WET_MASS_PG_PER_UM3 * CELL_DRY_MASS_FRACTION * LABILE_NECROMASS_FRACTION
                    food_value = max(0.0, b.food_level) * CYTOSOL_RECOVERY_FRACTION + necromass_pg
                    if food_value > MIN_FOOD_NODE_AMOUNT_PG:
                        food_nodes.append(FoodNode.acquire().init(b.x, b.y, food_value))
        # apply deaths and births
        for d in deaths:
            self.members.remove(d)
        self.members.extend(births)
        # record new population size
        self.history.append(len(self.members))
        # Track overall births/deaths
        self.total_births += len(births)
        self.total_deaths += len(deaths)
        return len(births), len(deaths)


class FoodSpatialHash:
    """
    Spatial hash table for efficient food node lookup by grid cell.
    Supports O(1) average case insertion/lookup/removal.
    """
    def __init__(self, cell_size=50):
        self.cell_size = cell_size
        self.grid = {}
    
    def insert(self, food_node):
        grid_x = int(food_node.x // self.cell_size)
        grid_y = int(food_node.y // self.cell_size)
        key = (grid_x, grid_y)
        if key not in self.grid:
            self.grid[key] = []
        self.grid[key].append(food_node)
    
    def query(self, x, y, radius):
        cell_x = int(x // self.cell_size)
        cell_y = int(y // self.cell_size)
        search_radius = int(radius / self.cell_size) + 1
        results = []
        for dx in range(-search_radius, search_radius + 1):
            for dy in range(-search_radius, search_radius + 1):
                key = (cell_x + dx, cell_y + dy)
                if key in self.grid:
                    results.extend(self.grid[key])
        return results

class Simulation:
    """
    Manages multiple species, environment, and global stats.
    """
    def __init__(
        self,
        spp,
        env,
        cap,
        steps=500,
        seed=None,
        physics_backend='auto',
        *,
        ecm_backend='auto',
        metabolic_adapter=None,
    ):
        # Phase I reproducibility: a non-None seed makes every stochastic draw deterministic so
        # extinction/abundance baselines are repeatable. Seed BEFORE any other RNG use so runs that
        # construct species after this point (programmatic tests) share the same stream.
        self.seed = seed
        if seed is not None:
            random.seed(seed)
        # give each species a reference to the full list
        for sp in spp:
            setattr(sp, 'species_list', spp)
            setattr(sp, 'sim', self)
        self.species = spp
        self.env = env           # environment dict: oxygen/CO2/food
        self.cap = cap           # carrying capacity
        self.steps = steps       # legacy max update count (not used by the real-time loop)
        self.ecm = ECM(
            cell_size=ECM_CELL_SIZE_UM,
            backend=ecm_backend,
        ) # Inject Extracellular Matrix for molecular diffusion
        self.physics = PhysicsOrchestrator(backend=physics_backend)
        self.metabolic_adapter = metabolic_adapter
        self._pending_diet_boluses = []
        self.last_fba_result = None
        self._last_interaction_summary = InteractionSummary()
        # stats collected each step
        self.stats = self._empty_stats()
        self.species_step_stats = []  # List of dicts: [{'births': [...], 'deaths': [...]}]
        self.step_idx = 0        # current time step
        self.running = False     # is simulation active?
        self.food_nodes = []
        self.drop_frequency = env.get('drop_frequency', 5)
        self.drop_range = env.get('drop_range', (2, 8))
        # Phase I: active diet profile controls the substrate-class mix of nutrient drops
        # (malnutrition / standard / mdcf). Falls back to the standard supplementary diet.
        self.diet_profile = env.get('diet_profile', DEFAULT_DIET_PROFILE)
        # Lazily-populated flag: the founder inoculum is placed by _seed_initial_population on the
        # first update (or reset), once the petri-dish geometry is established from the window size.
        self._populated = False

    def _empty_stats(self):
        return {
            'time': [],
            'elapsed_seconds': [],
            'dt_seconds': [],
            'tickrate_hz': [],
            'total': [],
            'births': [],
            'deaths': [],
            'overall_births': [],
            'overall_deaths': [],
            't6ss_deaths': [],
            'hgt_transfers': [],
            'overall_t6ss_deaths': [],
            'overall_hgt_transfers': [],
            # Phase I spatial heterogeneity: per-step population split between the central lumen
            # and the outer mucus annulus (membership by radial fraction vs MUCUS_LAYER_RADIUS_FRACTION).
            'lumen': [],
            'mucus': [],
            # Phase I diet analytics: cumulative substrate mass (pg) ingested by the whole community,
            # one running total per SUBSTRATE_* class (0 = generic calorie, 1 = complex glycan,
            # 2 = MDCF prebiotic). Monotonic non-decreasing; supports diet-efficacy comparisons.
            'consumed_calorie': [],
            'consumed_glycan': [],
            'consumed_mdcf': [],
            'fba_biomass_gdw': [],
            'fba_acetate_mol': [],
            'fba_consumed_calorie_mol': [],
            'fba_consumed_glycan_mol': [],
            'fba_consumed_mdcf_mol': [],
        }

    def queue_diet_bolus(self, bolus):
        """Queue one validated continuous-substrate bolus for the next FBA tick."""
        if self.metabolic_adapter is None:
            raise RuntimeError('diet boluses require an explicit live FBA adapter')
        if not isinstance(bolus, DietBolus):
            raise TypeError('bolus must be a DietBolus')
        self._pending_diet_boluses.append(bolus)

    def _ensure_petri_geometry(self, w, h):
        """Establish the circular petri-dish geometry (centre + radius) from the window size on the
        first call. The dish is anchored at the WORLD ORIGIN (0, 0) so it coincides with the grid
        camera's home pose (camera default + Home button both target (0, 0) at 1x zoom); this keeps
        the freshly seeded inoculum framed without any camera offset. The radius pads the window
        diagonal so the dish encloses the visible area. Idempotent so later resizes do not move an
        established dish."""
        if not hasattr(self, 'petri_center'):
            self.petri_center = (0.0, 0.0)
            self.petri_radius = math.hypot(w, h) / 2 + PETRI_RADIUS_PADDING_UM

    def _apply_spawn_env(self):
        """Push any spawn-configuration overrides carried in the environment dict (set from the
        Globals tab) onto the module-level spawn parameters that the founder-placement helpers read.
        Keeping the helpers reading module globals (rather than threading a config object through
        every call) preserves their simple, directly-testable signatures. Absent keys leave the
        compiled-in constant defaults untouched, and every value is validated/clamped to its valid
        domain so a malformed UI entry can never crash placement."""
        env = self.env
        g = globals()
        # Categorical geometry selectors: accept only recognised options, else keep the current value.
        mode = env.get('spawn_mode')
        if mode in SPAWN_MODES:
            g['INITIAL_SPAWN_MODE'] = mode
        arrangement = env.get('species_arrangement')
        if arrangement in SPECIES_ARRANGEMENTS:
            g['INITIAL_SPECIES_ARRANGEMENT'] = arrangement
        # Numeric controls: clamp to valid domains. Fractions are dimensionless on [0, 1]; the
        # cluster count is a positive integer number of founder microcolonies. Each conversion is
        # guarded independently so one bad field cannot discard the others.
        def _clamp01(v):
            return max(0.0, min(1.0, float(v)))
        for key, name, conv in (
            ('spawn_lumen_bias', 'INITIAL_SPAWN_LUMEN_BIAS', _clamp01),
            ('spawn_cluster_count', 'INITIAL_SPAWN_CLUSTER_COUNT', lambda v: max(1, int(float(v)))),
            ('spawn_cluster_spread', 'INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION', _clamp01),
            ('spawn_central_radius', 'INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION', _clamp01),
        ):
            raw = env.get(key)
            if raw is None or raw == '':
                continue
            try:
                g[name] = conv(raw)
            except (TypeError, ValueError):
                pass  # ignore unparseable entries, retaining the prior value

    def _seed_initial_population(self):
        """Place every species' founder inoculum using the configured spawning policy
        (INITIAL_SPAWN_MODE / INITIAL_SPECIES_ARRANGEMENT). Requires the petri geometry to exist.
        Idempotent within a run via the self._populated guard, and reproducible under the RNG seed,
        which is applied before this runs (construction or reset)."""
        # Apply any Globals-tab spawn overrides before drawing positions so live edits take effect
        # on the next reset/re-seed.
        self._apply_spawn_env()
        center, radius = self.petri_center, self.petri_radius
        total = len(self.species)
        shared_nuclei = None
        # 'mixed' arrangement co-locates species by sharing one set of founder nuclei so colonies of
        # different taxa overlap, maximising early contact-dependent interactions (T6SS, HGT, etc.).
        if (INITIAL_SPECIES_ARRANGEMENT == 'mixed'
                and INITIAL_SPAWN_MODE not in ('central', 'uniform_disc')):
            shared_nuclei = generate_cluster_nuclei(INITIAL_SPAWN_CLUSTER_COUNT, center, radius)
        for idx, sp in enumerate(self.species):
            positions = generate_initial_positions(
                sp.init_count, center, radius,
                species_index=idx, species_total=total, nuclei=shared_nuclei)
            sp.members = [Bacterium(sp, x, y,
                                    size=max(0.1, random.gauss(sp.base_size, sp.base_size * 0.2)))
                          for (x, y) in positions]
        self._populated = True

    def reset(self):
        """
        Reset populations to initial values and clear stats.
        """
        # Phase I reproducibility: re-seeding on reset makes repeated runs from the same seed
        # regenerate identical initial layouts and stochastic trajectories.
        if self.seed is not None:
            random.seed(self.seed)
        if hasattr(self, 'petri_center'):
            # Geometry already known: re-place the inoculum immediately via the spawning policy.
            self._seed_initial_population()
        else:
            # Geometry not yet established (reset before the first frame): defer placement to the
            # next update, which establishes the dish then seeds the founders.
            for sp in self.species:
                sp.members = []
            self._populated = False
        for sp in self.species:
            sp.history = []
            sp.total_births = 0
            sp.total_deaths = 0
            sp.extinction_time = None
        self.stats = self._empty_stats()
        self.species_step_stats = []
        self._last_interaction_summary = InteractionSummary()
        self._pending_diet_boluses.clear()
        self.last_fba_result = None
        if self.metabolic_adapter is not None:
            self.ecm.clear_substrate_layers()
            self.metabolic_adapter.state_store.clear()
            self.metabolic_adapter.orchestrator.clear_cache()
        self.step_idx = 0
        self.running = False

    def _process_interactions(self, spatial_hash=None, CELL_SIZE=None, dt=None):
        """
        Process ECM biochemistry, then dispatch ordered contact interactions.

        The legacy spatial_hash and CELL_SIZE arguments remain accepted for direct callers, but the
        authoritative contact pass builds its signed hash from the flattened arrays.
        """
        if dt is None:
            raise TypeError('dt is required')

        # This stable biochemical-then-contact sequence is the deterministic two-phase contract;
        # it intentionally does not promise historical full-tick trajectory parity. See future
        # DOCUMENTATION.md §10.2.
        for sp in self.species:
            # Bacteriocin production cost is expression times the per-second nutrient cost and dt;
            # see DOCUMENTATION.md §4.12.
            bac_prod_rate = (sp.bacteriocin_production * BACTERIOCIN_PRODUCTION_COST_RATE) * dt
            # Species immunity attenuates local bacteriocin dose up to the documented ceiling;
            # see DOCUMENTATION.md §4.12.
            immunity_mitigation = min(BACTERIOCIN_MAX_IMMUNITY_MITIGATION, sp.immunity_factor * BACTERIOCIN_IMMUNITY_SCALING)

            # M2 cross-feeding gates (computed once per species per tick).
            # Acetate oxidation is an aerobic, temperature-limited process, so secretion of overflow
            # acetate (a hallmark of fast aerobic glucose growth) and acetate uptake both scale with
            # the species' Ratkowsky temperature factor; uptake additionally requires dissolved O2.
            cf_temp_factor = sp.temperature_factor(self.env.get('temperature', DEFAULT_ENV_TEMPERATURE_C))
            cf_o2_gate = min(1.0, max(0.0, self.env.get('o2', ATMOSPHERIC_O2_FRACTION)) / AEROBE_O2_SATURATION_LEVEL)
            cf_secretes = sp.acetate_secretion > 0.0 and cf_temp_factor > 0.0
            cf_consumes = sp.acetate_uptake > 0.0 and cf_temp_factor > 0.0 and cf_o2_gate > 0.0

            for b in sp.members:
                if b.food_level < 0:
                    continue

                # Reset per-tick acetate consumption before any cross-feeding so the env update only
                # couples O2/CO2 to acetate actually oxidized this tick.
                b.last_acetate_consumption_pg = 0.0

                # --- M2 overflow acetate secretion ---
                # Excrete a fraction of the glucose taken up last tick as acetate (capped by the
                # 1 glucose -> 2 acetate mass ceiling and by the carbon the cell currently holds).
                if cf_secretes and b.last_food_uptake_pg > 0.0:
                    secreted_pg = b.last_food_uptake_pg * sp.acetate_secretion * ACETATE_PER_GLUCOSE_MAX_MASS_RATIO
                    secreted_pg = min(secreted_pg, max(0.0, b.food_level))
                    if secreted_pg > 0.0:
                        b.food_level -= secreted_pg
                        # Scale pg into ECM grid units so the trace survives the prune floor; inverted on uptake.
                        self.ecm.add(b.x, b.y, secreted_pg * ACETATE_ECM_UNITS_PER_PG, layer='acetate')

                # --- M2 acetate uptake (cross-feeding) ---
                # Oxidize locally available acetate as a secondary carbon source, limited by transporter
                # budget (temperature- and O2-gated), local availability, and remaining storage room.
                if cf_consumes and b.food_level < b.food_capacity:
                    available_pg = self.ecm.sample(b.x, b.y, layer='acetate') / ACETATE_ECM_UNITS_PER_PG
                    if available_pg > 0.0:
                        budget_pg = b.consumption_speed * sp.acetate_uptake * dt * cf_temp_factor * cf_o2_gate
                        # Room is expressed in acetate mass: each pg acetate stores ACETATE_ENERGY_YIELD pg glucose-equiv.
                        room_pg = (b.food_capacity - b.food_level) / ACETATE_ENERGY_YIELD_GLUCOSE_EQUIV
                        consumed_pg = max(0.0, min(available_pg, budget_pg, room_pg))
                        if consumed_pg > 0.0:
                            b.food_level += consumed_pg * ACETATE_ENERGY_YIELD_GLUCOSE_EQUIV
                            b.last_acetate_consumption_pg = consumed_pg
                            self.ecm.add(b.x, b.y, -consumed_pg * ACETATE_ECM_UNITS_PER_PG, layer='acetate')
                
                # Bacteriocin Emission (metabolic cost):
                # Emission occurs consistently, scaled by bacteriocin_production.
                # Ribosomally synthesized antimicrobial peptides; production costs ATP/amino acids. see DOCUMENTATION.md §4.12 (Competition, Bacteriocins, and T6SS).
                if bac_prod_rate > 0 and b.food_level > bac_prod_rate:
                    b.food_level -= bac_prod_rate
                    self.ecm.add(b.x, b.y, bac_prod_rate * BACTERIOCIN_ECM_SECRETION_MULTIPLIER, layer='bac')
                
                # Environmental bacteriocin exposure
                local_bac = self.ecm.sample(b.x, b.y, layer='bac')
                if local_bac > 0:
                    # Damage scaled by immunity factor acting as competitive inhibitor.
                    # M4: dormant persisters tolerate the antimicrobial dose (see DOCUMENTATION.md §4.17 (Phenotypic Switching)).
                    persister_resistance = PERSISTER_RESISTANCE_FACTOR if b.is_persister else 1.0
                    bac_dmg = (local_bac * BACTERIOCIN_DAMAGE_COEFFICIENT) * (1.0 - immunity_mitigation) * persister_resistance * dt
                    b.food_level -= bac_dmg

        summary = self.physics.process(self.species, dt)
        self._last_interaction_summary = summary
        return summary

    def _sample_substrate(self):
        """Draw a substrate class for one nutrient drop from the active diet profile.

        The diet profile (malnutrition / standard / mdcf) defines the categorical probability of
        each substrate class, reproducing distinct experimental feeding arms.
        Inverse-CDF sampling over the profile's ratios; unknown profiles fall back to generic calorie.
        see DOCUMENTATION.md §4.19 (Diet-Driven Substrate Metabolism).
        """
        ratios = DIET_PROFILES.get(self.diet_profile, DIET_PROFILES[DEFAULT_DIET_PROFILE])
        r = random.random()
        cumulative = 0.0
        for sid in SUBSTRATE_IDS:
            cumulative += ratios.get(sid, 0.0)
            if r < cumulative:
                return sid
        return DEFAULT_SUBSTRATE_ID

    def update(self, w, h, dt, viewport_bounds=None):
        """
        Advance one time step if running.
        viewport_bounds: dict with keys 'left', 'right', 'top', 'bottom' for visible area
        """
        self._ensure_petri_geometry(w, h)
        # Lazily place the founder inoculum once the dish geometry exists; this runs even while
        # paused so the starting layout is visible before the simulation is started.
        if not self._populated:
            self._seed_initial_population()

        if not self.running:
            return

        if self.metabolic_adapter is not None:
            return self._update_fba(w, h, dt, viewport_bounds)
            
        # Poisson-rate approximation: expected drops = lambda * dt, so spawning is real-time rather than FPS-bound.
        expected_drops = max(0.0, self.drop_frequency * dt)
        drop_count = int(expected_drops)
        if random.random() < (expected_drops - drop_count):
            drop_count += 1
        for _ in range(drop_count):
            r_f = self.petri_radius * math.sqrt(random.random())
            th = random.random() * 2 * math.pi
            x = self.petri_center[0] + r_f * math.cos(th)
            y = self.petri_center[1] + r_f * math.sin(th)
            amount = random.uniform(*self.drop_range)
            # Spatial heterogeneity: the mucus annulus receives less direct nutrient deposition than
            # the well-mixed lumen, which luminal flow supplies with most fresh dietary substrate.
            if r_f >= self.petri_radius * MUCUS_LAYER_RADIUS_FRACTION:
                amount *= MUCUS_DROP_AMOUNT_FACTOR
            # Phase I: assign a substrate class to this drop per the active diet profile.
            substrate_id = self._sample_substrate()
            node = FoodNode.acquire().init(x, y, amount, substrate_id)
            self.food_nodes.append(node)

        # --- M5 Phage seeding: Poisson introduction of free virions (epidemic trigger) ---
        # phage_inflow is the expected number of seed packets per second; each packet deposits a bolus
        # of virions at a random petri location. Default 0 means no outbreak; a slider can trigger or
        # halt epidemics live. Poisson rate (lambda * dt) keeps seeding real-time rather than FPS-bound.
        phage_inflow = max(0.0, float(self.env.get('phage_inflow', DEFAULT_PHAGE_INFLOW)))
        if phage_inflow > 0.0:
            expected_packets = phage_inflow * dt
            packet_count = int(expected_packets)
            if random.random() < (expected_packets - packet_count):
                packet_count += 1
            for _ in range(packet_count):
                r_p = self.petri_radius * math.sqrt(random.random())
                th_p = random.random() * 2 * math.pi
                px = self.petri_center[0] + r_p * math.cos(th_p)
                py = self.petri_center[1] + r_p * math.sin(th_p)
                self.ecm.add(px, py, PHAGE_SEED_PACKET_VIRIONS * PHAGE_ECM_UNITS_PER_VIRION, layer='phage')

        # Spatial washout: luminal flow / peristalsis clears free nutrient mass quickly, while the
        # mucus layer is a stabilized refuge. Exponential mass loss per real second, by zone.
        cx0, cy0 = self.petri_center
        inv_r = 1.0 / max(self.petri_radius, NUMERIC_DISTANCE_EPSILON_UM)
        for node in self.food_nodes:
            frac = math.hypot(node.x - cx0, node.y - cy0) * inv_r
            rate = MUCUS_FOOD_WASHOUT_RATE_PER_SEC if frac >= MUCUS_LAYER_RADIUS_FRACTION else LUMEN_FOOD_WASHOUT_RATE_PER_SEC
            node.amount -= node.amount * rate * dt

        to_remove = [n for n in self.food_nodes if n.amount < MIN_FOOD_NODE_AMOUNT_PG]
        for node in to_remove:
            self.food_nodes.remove(node)
            FoodNode.release(node)
            
        MAX_FOOD_NODES = 1000
        if len(self.food_nodes) > MAX_FOOD_NODES:
            excess = self.food_nodes[:-MAX_FOOD_NODES]
            for node in excess:
                FoodNode.release(node)
            self.food_nodes = self.food_nodes[-MAX_FOOD_NODES:]
            
        # Spatial Hash mapping
        food_hash = FoodSpatialHash(cell_size=FOOD_SPATIAL_HASH_CELL_SIZE_UM)
        for node in self.food_nodes:
            food_hash.insert(node)
        
        # Process Inter/Intra species interactions
        interaction_summary = self._process_interactions(dt=dt)

        # Tick ECM molecular diffusion mechanics
        self.ecm.tick(dt)

        bsum, dsum = 0, 0
        births_list, deaths_list = [], []
        for sp in self.species:
            b, d = sp.step(self.env, self.cap, w, h, self.food_nodes, dt, viewport_bounds, food_hash)
            bsum += b
            dsum += d
            births_list.append(b)
            deaths_list.append(d)
        # --- Enforce hard carrying capacity ---
        total_after = sum(len(sp.members) for sp in self.species)
        excess = total_after - self.cap
        forced_deaths_per_species = {sp: 0 for sp in self.species}
        if excess > 0:
            all_bacteria = []
            for sp in self.species:
                for b in sp.members:
                    all_bacteria.append((sp, b))
            if len(all_bacteria) > 0:
                n_kill = min(excess, len(all_bacteria))
                # Spatial competition: the well-mixed lumen imposes stronger resource competition
                # than the stabilized mucus refuge, so luminal cells are preferentially lost under
                # carrying-capacity pressure. Weighted sampling without replacement: key = U^(1/w);
                # the n_kill largest keys are removed, so lumen cells (up-weighted by
                # LUMEN_COMPETITION_CULL_WEIGHT) are more likely to be culled. With equal weights this
                # reduces exactly to uniform random sampling. see DOCUMENTATION.md §4.20 (Spatial Heterogeneity: Lumen vs Mucus Layer).
                cx0, cy0 = self.petri_center
                inv_r = 1.0 / max(self.petri_radius, NUMERIC_DISTANCE_EPSILON_UM)
                keyed = []
                for sp, b in all_bacteria:
                    frac = math.hypot(b.x - cx0, b.y - cy0) * inv_r
                    w_comp = 1.0 if frac >= MUCUS_LAYER_RADIUS_FRACTION else LUMEN_COMPETITION_CULL_WEIGHT
                    keyed.append((random.random() ** (1.0 / w_comp), sp, b))
                keyed.sort(key=lambda t: t[0], reverse=True)
                to_kill = [(sp, b) for _, sp, b in keyed[:n_kill]]
                for sp, b in to_kill:
                    sp.members.remove(b)
                    sp.total_deaths += 1
                    forced_deaths_per_species[sp] += 1
        for idx, sp in enumerate(self.species):
            deaths_list[idx] += forced_deaths_per_species[sp]
        sum_forced = sum(forced_deaths_per_species.values())
        elapsed_seconds = (self.stats['elapsed_seconds'][-1] if self.stats['elapsed_seconds'] else 0.0) + dt
        # Phase I extinction baseline: record the elapsed sim-time at which each species first hits
        # zero population. A species cannot regrow from zero in this model, so the first crossing is
        # the time-to-extinction; survivors keep extinction_time = None.
        for sp in self.species:
            if sp.extinction_time is None and len(sp.members) == 0:
                sp.extinction_time = elapsed_seconds
        step_births = bsum
        step_deaths = dsum + sum_forced
        overall_births = (self.stats['overall_births'][-1] if self.stats['overall_births'] else 0) + step_births
        overall_deaths = (self.stats['overall_deaths'][-1] if self.stats['overall_deaths'] else 0) + step_deaths
        # Contact-event totals add this ordered pass's transition counts to the previous tick's
        # cumulative values; see DOCUMENTATION.md §10.2.
        overall_t6ss_deaths = (self.stats['overall_t6ss_deaths'][-1] if self.stats['overall_t6ss_deaths'] else 0) + interaction_summary.t6ss_deaths
        overall_hgt_transfers = (self.stats['overall_hgt_transfers'][-1] if self.stats['overall_hgt_transfers'] else 0) + interaction_summary.hgt_transfers
        self.stats['time'].append(self.step_idx)
        self.stats['elapsed_seconds'].append(elapsed_seconds)
        self.stats['dt_seconds'].append(dt)
        self.stats['tickrate_hz'].append((1.0 / dt) if dt > 0 else 0.0)
        self.stats['total'].append(sum(len(sp.members) for sp in self.species))
        self.stats['births'].append(step_births)
        self.stats['deaths'].append(step_deaths)
        self.stats['overall_births'].append(overall_births)
        self.stats['overall_deaths'].append(overall_deaths)
        self.stats['t6ss_deaths'].append(interaction_summary.t6ss_deaths)
        self.stats['hgt_transfers'].append(interaction_summary.hgt_transfers)
        self.stats['overall_t6ss_deaths'].append(overall_t6ss_deaths)
        self.stats['overall_hgt_transfers'].append(overall_hgt_transfers)
        # Phase I spatial heterogeneity: classify each live cell by radial fraction so the Stats tab
        # and export can report lumen vs mucus occupancy. mucus = outer annulus (frac >= threshold).
        cx_z, cy_z = self.petri_center
        inv_r_z = 1.0 / max(self.petri_radius, NUMERIC_DISTANCE_EPSILON_UM)
        mucus_count = 0
        total_count = 0
        # Phase I diet analytics: sum this step's per-cell substrate uptake, then carry forward the
        # previous cumulative total so the series is monotonic (total mass ever ingested per class).
        step_uptake = [0.0, 0.0, 0.0]  # indices map to SUBSTRATE_* ids 0/1/2
        for sp in self.species:
            for b in sp.members:
                total_count += 1
                if math.hypot(b.x - cx_z, b.y - cy_z) * inv_r_z >= MUCUS_LAYER_RADIUS_FRACTION:
                    mucus_count += 1
                ub = getattr(b, 'last_uptake_by_substrate', None)
                if ub:
                    for _si in range(min(len(ub), len(step_uptake))):
                        step_uptake[_si] += ub[_si]
        self.stats['mucus'].append(mucus_count)
        self.stats['lumen'].append(total_count - mucus_count)
        prev_cal = self.stats['consumed_calorie'][-1] if self.stats['consumed_calorie'] else 0.0
        prev_gly = self.stats['consumed_glycan'][-1] if self.stats['consumed_glycan'] else 0.0
        prev_mdcf = self.stats['consumed_mdcf'][-1] if self.stats['consumed_mdcf'] else 0.0
        self.stats['consumed_calorie'].append(prev_cal + step_uptake[0])
        self.stats['consumed_glycan'].append(prev_gly + step_uptake[1])
        self.stats['consumed_mdcf'].append(prev_mdcf + step_uptake[2])
        self.species_step_stats.append({'births': births_list, 'deaths': deaths_list})
        for sp in self.species:
            if sp.history:
                sp.history[-1] = len(sp.members)
            else:
                sp.history.append(len(sp.members))
                
        # Real-world Environmental Exchange & Metabolism
        # Unit Definitions Context: 1 food unit maps to 1 pg glucose-equivalent organic nutrient.
        # Environment is modeled as 1 mL of fluid (10^12 um^3).
        # Atmospheric O2 dissolved in 1 mL water represents about 8.4 ug of O2 (8.4e6 pg).
        # `env['o2']` of 1.0 and `env['co2']` of 1.0 use the documented dissolved gas capacities.
        
        o2_consumption_fraction = 0.0
        co2_production_fraction = 0.0
        acid_production = 0.0
        
        for sp in self.species:
            pop = len(sp.members)
            if pop == 0: continue
            
            # total_glucose_pg: summed glucose-equivalent nutrient actually consumed by size- and age-scaled cells this update.
            total_glucose_pg = sum(getattr(b, 'last_metabolic_consumption_pg', 0.0) for b in sp.members)
            
            if sp.o2_sens > 0:
                # Aerobic Respiration: C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O
                # Stoichiometry mapping cell metabolism to gas usage exactly; see DOCUMENTATION.md §4.10 (Gas Exchange, Fermentation, and pH).
                # 1 mole Glucose (180.16 g) consumes 6 moles O2 (192.0 g) and produces 6 moles CO2 (264.0 g).
                o2_consumed_pg = total_glucose_pg * (O2_SIX_MOLES_MASS_G / GLUCOSE_MOLAR_MASS_G) * sp.o2_sens
                co2_produced_pg = total_glucose_pg * (CO2_SIX_MOLES_MASS_G / GLUCOSE_MOLAR_MASS_G) * sp.o2_sens
                
                o2_consumption_fraction += o2_consumed_pg / ENV_100_O2_PG
                co2_production_fraction += co2_produced_pg / ENV_100_CO2_PG
                
                # Carbonic acid (H2CO3) weakly dissociates (pKa 6.35). Moderate effect on environmental pH.
                acid_production += (co2_produced_pg / ENV_100_CO2_PG) * CARBONIC_ACID_PH_EFFECT
                
            elif sp.o2_sens < 0:
                # Anaerobic Fermentation (Lactic Acid Fermentation): C6H12O6 -> 2 Lactic Acid
                # Stoichiometry mapping homolactic fermentation for anaerobic life; see DOCUMENTATION.md §4.10 (Gas Exchange, Fermentation, and pH).
                # 1 mole Glucose (180.16 g) produces 2 moles Lactic Acid (180.16 g, conservation of mass).
                lactic_acid_pg = total_glucose_pg * abs(sp.o2_sens)
                
                # Small trace amounts of CO2 depending on fermentation pathway (e.g. heterolactic models).
                co2_produced_pg = total_glucose_pg * HETEROLACTIC_CO2_TRACE_FRACTION * abs(sp.o2_sens)
                co2_production_fraction += co2_produced_pg / ENV_100_CO2_PG
                
                # Lactic acid strongly drops pH due to highly acidic dissociation equilibrium (pKa 3.86).
                # Henderson-Hasselbalch representation inside small buffered environments; see DOCUMENTATION.md §4.10 (Gas Exchange, Fermentation, and pH).
                acid_production += (lactic_acid_pg / ENV_100_O2_PG) * LACTIC_ACID_PH_EFFECT

            # M2 acetate oxidation gas coupling: C2H4O2 + 2 O2 -> 2 CO2 + 2 H2O.
            # Acetate cross-feeders consume extra O2 and produce CO2 (and carbonic acid) independent of
            # the glucose-driven respiration above, so this is summed for every species that ate acetate.
            total_acetate_pg = sum(getattr(b, 'last_acetate_consumption_pg', 0.0) for b in sp.members)
            if total_acetate_pg > 0.0:
                acetate_o2_pg = total_acetate_pg * O2_PER_ACETATE_MASS_RATIO
                acetate_co2_pg = total_acetate_pg * CO2_PER_ACETATE_MASS_RATIO
                o2_consumption_fraction += acetate_o2_pg / ENV_100_O2_PG
                co2_production_fraction += acetate_co2_pg / ENV_100_CO2_PG
                acid_production += (acetate_co2_pg / ENV_100_CO2_PG) * CARBONIC_ACID_PH_EFFECT


        self.env['o2'] -= o2_consumption_fraction
        self.env['co2'] += co2_production_fraction
        self.env['ph'] -= acid_production

        # Atmospheric gas dissolution and physiological pH buffering mechanics.
        # Fick's first law of diffusion (J = -D * dC/dx): concentration gradient drives gas toward atmospheric equilibrium.
        # Assuming a shallow 1 mL volume open to sea-level atmospheric pressure, equilibrium typically reaches ~30% per second. see DOCUMENTATION.md §4.10 (Gas Exchange, Fermentation, and pH).
        diffusion_rate = GAS_DIFFUSION_RATE_SEC * dt
        self.env['o2'] += (ATMOSPHERIC_O2_FRACTION - self.env['o2']) * diffusion_rate
        self.env['co2'] += (ATMOSPHERIC_CO2_FRACTION - self.env['co2']) * diffusion_rate
        
        # Buffer system (e.g. bicarbonate buffer system naturally driving pH towards neutral 7 in open systems).
        # Henderson-Hasselbalch equation and Le Chatelier's equilibrium principle in aqueous environments. see DOCUMENTATION.md §4.10 (Gas Exchange, Fermentation, and pH).
        buffer_rate = BUFFER_RESTORATION_RATE_SEC * dt
        self.env['ph'] += (NEUTRAL_PH - self.env['ph']) * buffer_rate

        # Cap bounds
        self.env['o2'] = max(0.0, min(1.0, self.env['o2']))
        self.env['co2'] = max(0.0, min(1.0, self.env['co2']))
        self.env['ph'] = max(0.0, min(14.0, self.env['ph']))

        self.step_idx += 1

    def _snapshot_fba_transaction(self):
        agents = tuple(
            bacterium
            for species in self.species
            for bacterium in species.members
        )
        identity_objects = (
            self,
            self.ecm,
            self.metabolic_adapter,
            self.metabolic_adapter.state_store,
            *self.species,
            *agents,
        )
        memo = {id(value): value for value in identity_objects}
        return {
            'substrate_snapshot': self.ecm.snapshot_substrate_layers(),
            'signal_layers': tuple(
                (name, dict(self.ecm._layers[name]))
                for name in self.ecm._layers
            ),
            'ecm_backend_used': self.ecm._backend_used,
            'state_store': self.metabolic_adapter.state_store.snapshot(),
            'solver_cache': self.metabolic_adapter.orchestrator.snapshot_cache(),
            'species': tuple(
                (species, copy.deepcopy(species.__dict__, memo))
                for species in self.species
            ),
            'agents': tuple(
                (agent, copy.deepcopy(agent.__dict__, memo))
                for agent in agents
            ),
            'stats': copy.deepcopy(self.stats),
            'species_step_stats': copy.deepcopy(self.species_step_stats),
            'environment': copy.deepcopy(self.env),
            'pending_boluses': tuple(self._pending_diet_boluses),
            'last_fba_result': self.last_fba_result,
            'step_idx': self.step_idx,
            'random_state': random.getstate(),
        }

    def _restore_fba_transaction(self, snapshot):
        self.ecm.restore_substrate_layers(snapshot['substrate_snapshot'])
        for name, values in snapshot['signal_layers']:
            layer = self.ecm._layers[name]
            layer.clear()
            layer.update(values)
        self.ecm._backend_used = snapshot['ecm_backend_used']
        self.metabolic_adapter.state_store.restore(snapshot['state_store'])
        self.metabolic_adapter.orchestrator.restore_cache(snapshot['solver_cache'])
        for species, attributes in snapshot['species']:
            species.__dict__.clear()
            species.__dict__.update(attributes)
        for agent, attributes in snapshot['agents']:
            agent.__dict__.clear()
            agent.__dict__.update(attributes)
        self.stats.clear()
        self.stats.update(snapshot['stats'])
        self.species_step_stats[:] = snapshot['species_step_stats']
        self.env.clear()
        self.env.update(snapshot['environment'])
        self._pending_diet_boluses[:] = snapshot['pending_boluses']
        self.last_fba_result = snapshot['last_fba_result']
        self.step_idx = snapshot['step_idx']
        random.setstate(snapshot['random_state'])

    def _update_fba(self, w, h, dt, viewport_bounds=None):
        snapshot = self._snapshot_fba_transaction()
        try:
            result = self._advance_fba_transaction(w, h, dt, viewport_bounds)
            self.metabolic_adapter.emit_success(result)
            return result
        except BaseException as error:
            self._restore_fba_transaction(snapshot)
            self.metabolic_adapter.emit_failure(
                snapshot['step_idx'],
                'simulation',
                error,
            )
            raise

    def _advance_fba_transaction(self, w, h, dt, viewport_bounds=None):
        """Advance one opt-in molar FBA tick without invoking legacy pg nutrition."""
        self.ecm.tick(dt)
        result = self.metabolic_adapter.advance_tick(
            species=self.species,
            ecm=self.ecm,
            tick_index=self.step_idx,
            dt_seconds=dt,
            boluses=tuple(self._pending_diet_boluses),
            defer_events=True,
        )
        self._pending_diet_boluses.clear()
        self.last_fba_result = result

        total_before_lifecycle = sum(len(species.members) for species in self.species)
        safe_capacity = self.cap if self.cap > 0 else 1
        carrying_capacity_factor = max(
            0.0,
            1.0 - total_before_lifecycle / safe_capacity,
        )
        births_list = []
        deaths_list = []
        for species in self.species:
            current_ph = self.env.get('ph', NEUTRAL_PH)
            ph_difference = current_ph - species.preferred_ph
            ph_suitability = math.exp(
                -0.5 * (ph_difference / PH_TOLERANCE_GAUSSIAN_SD) ** 2
            )
            oxygen = self.env.get('o2', ATMOSPHERIC_O2_FRACTION)
            if species.o2_sens > 0.0:
                oxygen_suitability = (
                    min(1.0, oxygen / AEROBE_O2_SATURATION_LEVEL) * species.o2_sens
                    + (1.0 - species.o2_sens)
                )
            elif species.o2_sens < 0.0:
                oxygen_suitability = 1.0 - (
                    min(1.0, oxygen / ANAEROBE_O2_TOXICITY_THRESHOLD)
                    * abs(species.o2_sens)
                )
            else:
                oxygen_suitability = 1.0
            environmental_fitness = ph_suitability * max(0.0, oxygen_suitability)
            temperature_factor = species.temperature_factor(
                self.env.get('temperature', DEFAULT_ENV_TEMPERATURE_C)
            )

            births = []
            deaths = []
            for bacterium in tuple(species.members):
                state = self.metabolic_adapter.state_store.state_for(bacterium)
                physiology = self.metabolic_adapter.lifecycle_policy.evaluate(
                    bacterium,
                    state,
                    dt_seconds=dt,
                    environmental_fitness=environmental_fitness,
                    carrying_capacity_factor=carrying_capacity_factor,
                )
                if not isinstance(physiology, FbaPhysiology):
                    raise TypeError('FBA lifecycle evaluate must return FbaPhysiology')
                bacterium.step_fba_kinematics(
                    self.ecm,
                    (self.petri_center[0], self.petri_center[1], self.petri_radius),
                    dt,
                    physiology,
                    viewport_bounds=viewport_bounds,
                    temp_factor=temperature_factor,
                )
                if physiology.starved:
                    deaths.append(bacterium)
                    continue
                if physiology.division_probability > 0.0 and random.random() < min(
                    1.0,
                    physiology.division_probability,
                ):
                    child = Bacterium(
                        species,
                        bacterium.x,
                        bacterium.y,
                        size=max(
                            0.1,
                            random.gauss(
                                bacterium.base_size,
                                bacterium.base_size * 0.2,
                            ),
                        ),
                        initial_q=bacterium.q,
                    )
                    self.metabolic_adapter.state_store.partition_for_fission(
                        bacterium,
                        child,
                        self.metabolic_adapter.lifecycle_policy,
                    )
                    births.append(child)
            for bacterium in deaths:
                species.members.remove(bacterium)
                self.metabolic_adapter.state_store.remove(bacterium)
            species.members.extend(births)
            species.total_births += len(births)
            species.total_deaths += len(deaths)
            species.history.append(len(species.members))
            births_list.append(len(births))
            deaths_list.append(len(deaths))

        excess = sum(len(species.members) for species in self.species) - self.cap
        if excess > 0:
            candidates = [
                (species, bacterium)
                for species in self.species
                for bacterium in species.members
            ]
            for species, bacterium in random.sample(candidates, min(excess, len(candidates))):
                species.members.remove(bacterium)
                species.total_deaths += 1
                self.metabolic_adapter.state_store.remove(bacterium)
                species_index = self.species.index(species)
                deaths_list[species_index] += 1

        elapsed_seconds = (
            self.stats['elapsed_seconds'][-1]
            if self.stats['elapsed_seconds']
            else 0.0
        ) + dt
        for species in self.species:
            if species.extinction_time is None and not species.members:
                species.extinction_time = elapsed_seconds
            species.history[-1] = len(species.members)

        step_births = sum(births_list)
        step_deaths = sum(deaths_list)
        self.stats['time'].append(self.step_idx)
        self.stats['elapsed_seconds'].append(elapsed_seconds)
        self.stats['dt_seconds'].append(dt)
        self.stats['tickrate_hz'].append((1.0 / dt) if dt > 0.0 else 0.0)
        self.stats['total'].append(sum(len(species.members) for species in self.species))
        self.stats['births'].append(step_births)
        self.stats['deaths'].append(step_deaths)
        self.stats['overall_births'].append(
            (self.stats['overall_births'][-1] if self.stats['overall_births'] else 0)
            + step_births
        )
        self.stats['overall_deaths'].append(
            (self.stats['overall_deaths'][-1] if self.stats['overall_deaths'] else 0)
            + step_deaths
        )
        for name in (
            't6ss_deaths',
            'hgt_transfers',
            'overall_t6ss_deaths',
            'overall_hgt_transfers',
        ):
            previous = self.stats[name][-1] if self.stats[name] else 0
            self.stats[name].append(previous if name.startswith('overall_') else 0)

        mucus_count = 0
        total_count = 0
        inverse_radius = 1.0 / max(self.petri_radius, NUMERIC_DISTANCE_EPSILON_UM)
        for species in self.species:
            for bacterium in species.members:
                total_count += 1
                if math.hypot(
                    bacterium.x - self.petri_center[0],
                    bacterium.y - self.petri_center[1],
                ) * inverse_radius >= MUCUS_LAYER_RADIUS_FRACTION:
                    mucus_count += 1
        self.stats['mucus'].append(mucus_count)
        self.stats['lumen'].append(total_count - mucus_count)
        for name in ('consumed_calorie', 'consumed_glycan', 'consumed_mdcf'):
            self.stats[name].append(self.stats[name][-1] if self.stats[name] else 0.0)

        self.stats['fba_biomass_gdw'].append(
            sum(
                self.metabolic_adapter.state_store.state_for(bacterium).biomass_gdw
                for species in self.species
                for bacterium in species.members
            )
        )
        self.stats['fba_acetate_mol'].append(result.acetate_secreted_mol)
        self.stats['fba_consumed_calorie_mol'].append(result.substrate_consumed_mol[0])
        self.stats['fba_consumed_glycan_mol'].append(result.substrate_consumed_mol[1])
        self.stats['fba_consumed_mdcf_mol'].append(result.substrate_consumed_mol[2])
        self.species_step_stats.append({'births': births_list, 'deaths': deaths_list})

        gas_diffusion = GAS_DIFFUSION_RATE_SEC * dt
        self.env['o2'] += (ATMOSPHERIC_O2_FRACTION - self.env['o2']) * gas_diffusion
        self.env['co2'] += (ATMOSPHERIC_CO2_FRACTION - self.env['co2']) * gas_diffusion
        buffer_rate = BUFFER_RESTORATION_RATE_SEC * dt
        self.env['ph'] += (NEUTRAL_PH - self.env['ph']) * buffer_rate
        self.env['o2'] = max(0.0, min(1.0, self.env['o2']))
        self.env['co2'] = max(0.0, min(1.0, self.env['co2']))
        self.env['ph'] = max(0.0, min(14.0, self.env['ph']))
        self.step_idx += 1
        return result

    def export(self, dm_species=None, default_env=None):
        """
        Write collected stats and per-species data to CSV file for analysis.
        """
        if not self.stats['time']:
            print("No data to export.")
            return
        with open(EXPORT_FILE, 'w', newline='') as f:
            writer = csv.writer(f)
            if default_env is not None:
                writer.writerow(['--- Environment Configuration ---'])
                for k, v in default_env.items():
                    writer.writerow([k, v])
                # Phase I reproducibility: persist the RNG seed used for this run so the
                # baseline can be exactly reproduced from the exported configuration.
                writer.writerow(['seed', self.seed if self.seed is not None else ''])
                writer.writerow([])
            if dm_species is not None and len(dm_species) > 0:
                writer.writerow(['--- Species Configuration ---'])
                headers = list(dm_species[0].keys())
                writer.writerow(headers)
                for sp in dm_species:
                    writer.writerow([sp.get(h, '') for h in headers])
                writer.writerow([])
            # Write header for global stats
            writer.writerow(['--- Global Stats ---'])
            writer.writerow([
                'time', 'total', 'births', 'deaths', 'elapsed_seconds', 'dt_seconds',
                'tickrate_hz', 'overall_births', 'overall_deaths', 't6ss_deaths',
                'hgt_transfers', 'overall_t6ss_deaths', 'overall_hgt_transfers',
                'lumen', 'mucus', 'consumed_calorie', 'consumed_glycan',
                'consumed_mdcf', 'fba_biomass_gdw', 'fba_acetate_mol',
                'fba_consumed_calorie_mol', 'fba_consumed_glycan_mol',
                'fba_consumed_mdcf_mol',
            ])
            # Zone and substrate-consumption columns may be absent on stats built by legacy imports
            # or manual construction; fall back to blank so the row width stays consistent without
            # raising IndexError.
            lumen_series = self.stats.get('lumen', [])
            mucus_series = self.stats.get('mucus', [])
            cal_series = self.stats.get('consumed_calorie', [])
            gly_series = self.stats.get('consumed_glycan', [])
            mdcf_series = self.stats.get('consumed_mdcf', [])
            t6ss_deaths_series = self.stats.get('t6ss_deaths', [])
            hgt_transfers_series = self.stats.get('hgt_transfers', [])
            overall_t6ss_deaths_series = self.stats.get('overall_t6ss_deaths', [])
            overall_hgt_transfers_series = self.stats.get('overall_hgt_transfers', [])
            fba_biomass_series = self.stats.get('fba_biomass_gdw', [])
            fba_acetate_series = self.stats.get('fba_acetate_mol', [])
            fba_calorie_series = self.stats.get('fba_consumed_calorie_mol', [])
            fba_glycan_series = self.stats.get('fba_consumed_glycan_mol', [])
            fba_mdcf_series = self.stats.get('fba_consumed_mdcf_mol', [])
            for i in range(len(self.stats['time'])):
                writer.writerow([
                    self.stats['time'][i],
                    self.stats['total'][i],
                    self.stats['births'][i],
                    self.stats['deaths'][i],
                    self.stats['elapsed_seconds'][i],
                    self.stats['dt_seconds'][i],
                    self.stats['tickrate_hz'][i],
                    self.stats['overall_births'][i],
                    self.stats['overall_deaths'][i],
                    t6ss_deaths_series[i] if i < len(t6ss_deaths_series) else 0,
                    hgt_transfers_series[i] if i < len(hgt_transfers_series) else 0,
                    overall_t6ss_deaths_series[i] if i < len(overall_t6ss_deaths_series) else 0,
                    overall_hgt_transfers_series[i] if i < len(overall_hgt_transfers_series) else 0,
                    lumen_series[i] if i < len(lumen_series) else '',
                    mucus_series[i] if i < len(mucus_series) else '',
                    cal_series[i] if i < len(cal_series) else '',
                    gly_series[i] if i < len(gly_series) else '',
                    mdcf_series[i] if i < len(mdcf_series) else '',
                    fba_biomass_series[i] if i < len(fba_biomass_series) else '',
                    fba_acetate_series[i] if i < len(fba_acetate_series) else '',
                    fba_calorie_series[i] if i < len(fba_calorie_series) else '',
                    fba_glycan_series[i] if i < len(fba_glycan_series) else '',
                    fba_mdcf_series[i] if i < len(fba_mdcf_series) else '',
                ])
            # Write per-species stats
            writer.writerow([])
            writer.writerow(['Index', 'Species', 'Step', 'Count', 'Births', 'Deaths', 'Overall Births', 'Overall Deaths'])
            for idx, sp in enumerate(self.species):
                cum_births = 0
                cum_deaths = 0
                for step in range(len(sp.history)):
                    births = self.species_step_stats[step]['births'][idx] if step < len(self.species_step_stats) else 0
                    deaths = self.species_step_stats[step]['deaths'][idx] if step < len(self.species_step_stats) else 0
                    cum_births += births
                    cum_deaths += deaths
                    writer.writerow([
                        idx,
                        sp.name,
                        step,
                        sp.history[step] if step < len(sp.history) else 0,
                        births,
                        deaths,
                        cum_births,
                        cum_deaths
                    ])
            # Phase I extinction baseline: one row per species with the time-to-extinction (elapsed
            # sim-seconds at first zero-population crossing) or blank if the species survived.
            writer.writerow([])
            writer.writerow(['--- Extinction Events ---'])
            writer.writerow(['Species', 'Time_To_Extinction_s', 'Extinct', 'Final_Count'])
            for sp in self.species:
                extinct = sp.extinction_time is not None
                writer.writerow([
                    sp.name,
                    sp.extinction_time if extinct else '',
                    int(extinct),
                    len(sp.members)
                ])
        print(f"Exported data to {EXPORT_FILE}")

    def import_stats(self, filename=EXPORT_FILE):
        """
        Import simulation stats from CSV and update simulation state.
        Only updates global stats and per-species history if format matches.
        """
        if not os.path.exists(filename):
            print(f"File {filename} not found.")
            return
        with open(filename, 'r', newline='') as f:
            reader = csv.reader(f)
            lines = list(reader)
        # Parse global stats using headers for forward/backward CSV compatibility.
        stats_header = []
        stats_start = 1
        for i, row in enumerate(lines):
            if row and row[0] == 'time' and len(row) >= 4 and row[1] == 'total':
                stats_header = row
                stats_start = i + 1
                break
        
        stats_end = stats_start
        while stats_end < len(lines) and lines[stats_end]:
            stats_end += 1
        self.stats = self._empty_stats()
        header_index = {name: idx for idx, name in enumerate(stats_header)}
        cumulative_births = 0
        cumulative_deaths = 0
        cumulative_t6ss_deaths = 0
        cumulative_hgt_transfers = 0
        for row in lines[stats_start:stats_end]:
            if len(row) >= 4:
                step_time = int(row[header_index.get('time', 0)])
                births = int(row[header_index.get('births', 2)])
                deaths = int(row[header_index.get('deaths', 3)])
                t6ss_deaths = int(row[header_index['t6ss_deaths']]) if 't6ss_deaths' in header_index and row[header_index['t6ss_deaths']] else 0
                hgt_transfers = int(row[header_index['hgt_transfers']]) if 'hgt_transfers' in header_index and row[header_index['hgt_transfers']] else 0
                cumulative_births += births
                cumulative_deaths += deaths
                cumulative_t6ss_deaths += t6ss_deaths
                cumulative_hgt_transfers += hgt_transfers
                elapsed = float(row[header_index['elapsed_seconds']]) if 'elapsed_seconds' in header_index and row[header_index['elapsed_seconds']] else (step_time + 1) / FPS
                dt_value = float(row[header_index['dt_seconds']]) if 'dt_seconds' in header_index and row[header_index['dt_seconds']] else 1.0 / FPS
                tickrate = float(row[header_index['tickrate_hz']]) if 'tickrate_hz' in header_index and row[header_index['tickrate_hz']] else FPS
                self.stats['time'].append(step_time)
                self.stats['total'].append(int(row[header_index.get('total', 1)]))
                self.stats['births'].append(births)
                self.stats['deaths'].append(deaths)
                self.stats['elapsed_seconds'].append(elapsed)
                self.stats['dt_seconds'].append(dt_value)
                self.stats['tickrate_hz'].append(tickrate)
                self.stats['overall_births'].append(int(row[header_index['overall_births']]) if 'overall_births' in header_index and row[header_index['overall_births']] else cumulative_births)
                self.stats['overall_deaths'].append(int(row[header_index['overall_deaths']]) if 'overall_deaths' in header_index and row[header_index['overall_deaths']] else cumulative_deaths)
                self.stats['t6ss_deaths'].append(t6ss_deaths)
                self.stats['hgt_transfers'].append(hgt_transfers)
                self.stats['overall_t6ss_deaths'].append(int(row[header_index['overall_t6ss_deaths']]) if 'overall_t6ss_deaths' in header_index and row[header_index['overall_t6ss_deaths']] else cumulative_t6ss_deaths)
                self.stats['overall_hgt_transfers'].append(int(row[header_index['overall_hgt_transfers']]) if 'overall_hgt_transfers' in header_index and row[header_index['overall_hgt_transfers']] else cumulative_hgt_transfers)
                for field_name in (
                    'fba_biomass_gdw',
                    'fba_acetate_mol',
                    'fba_consumed_calorie_mol',
                    'fba_consumed_glycan_mol',
                    'fba_consumed_mdcf_mol',
                ):
                    field_index = header_index.get(field_name)
                    value = (
                        float(row[field_index])
                        if field_index is not None and row[field_index]
                        else 0.0
                    )
                    self.stats[field_name].append(value)
        # Parse per-species stats
        species_start = stats_end + 2
        self.species_step_stats = []
        for sp in self.species:
            sp.history = []
            sp.total_births = 0
            sp.total_deaths = 0

        for row in lines[species_start:]:
            if len(row) >= 8 and row[0].isdigit():
                sp_idx, name, step, count, births, deaths, overall_births, overall_deaths = row
                # Find species by index
                sp_idx = int(sp_idx)
                if sp_idx < len(self.species):
                    sp = self.species[sp_idx]
                    step_idx = int(step)
                    while len(sp.history) <= step_idx:
                        sp.history.append(0)
                    sp.history[step_idx] = int(count) if count else 0
                    sp.total_births = int(overall_births) if overall_births else 0
                    sp.total_deaths = int(overall_deaths) if overall_deaths else 0
                    # Update species_step_stats
                    while len(self.species_step_stats) <= step_idx:
                        self.species_step_stats.append({'births': [0]*len(self.species), 'deaths': [0]*len(self.species)})
                    self.species_step_stats[step_idx]['births'][sp_idx] = int(births) if births else 0
                    self.species_step_stats[step_idx]['deaths'][sp_idx] = int(deaths) if deaths else 0

            elif len(row) >= 7 and not row[0].isdigit():
                name, step, count, births, deaths, overall_births, overall_deaths = row
                # Find species by name
                for sp in self.species:
                    if sp.name == name:
                        step_idx = int(step)
                        while len(sp.history) <= step_idx:
                            sp.history.append(0)
                        sp.history[step_idx] = int(count) if count else 0
                        sp.total_births = int(overall_births) if overall_births else 0
                        sp.total_deaths = int(overall_deaths) if overall_deaths else 0
                        # Update species_step_stats
                        while len(self.species_step_stats) <= step_idx:
                            self.species_step_stats.append({'births': [0]*len(self.species), 'deaths': [0]*len(self.species)})
                        idx = self.species.index(sp)
                        self.species_step_stats[step_idx]['births'][idx] = int(births) if births else 0
                        self.species_step_stats[step_idx]['deaths'][idx] = int(deaths) if deaths else 0

        if self.stats['time']:
            self.step_idx = self.stats['time'][-1] + 1
        else:
            self.step_idx = 0

        print(f"Imported stats from {filename}.")

