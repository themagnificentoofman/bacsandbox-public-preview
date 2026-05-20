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

class FoodNode:
    """
    Represents a food node on the grid.
    """
    _pool = []
    
    def __init__(self, x, y, amount):
        self.x = x
        self.y = y
        self.amount = amount
    
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
    
    def init(self, x, y, amount):
        """Initialize/reinitialize pool object."""
        self.x = x
        self.y = y
        self.amount = amount
        return self

class Bacterium:
    """
    Individual bacterium with position and species reference.
    """
    _pool = []
    
    def __init__(self, species, x, y, size=1.0, initial_food=None):
        self.species = species
        self.x = x
        self.y = y
        self.base_size = size  # Unit: micrometers (um)
        # Volumetric scaling (size^3) assuming constant cellular density context.
        # Unit mappings: 1 distance unit = 1 um. 1 food unit = 1 pg glucose-equivalent organic nutrient.
        # Reference: Simple geometric scaling for 3D objects; ~1 pg wet mass per 1 um^3 bacterial volume.
        self.food_capacity = getattr(species, 'food_capacity', DEFAULT_FOOD_CAPACITY_PG) * (self.base_size ** 3)
        self.food_level = initial_food if initial_food is not None else self.food_capacity
        
        # Kleiber's Law (Allometric scaling): Metabolic rate scales to Mass^(3/4)
        # Reference: Max Kleiber (1932), "Body size and metabolism". If Mass ~ size^3, then Meta ~ size^(2.25)
        # Units: pg glucose-equivalent nutrient consumed per second, later scaled by elapsed dt.
        self.metabolic_rate = getattr(species, 'metabolic_rate', DEFAULT_METABOLIC_RATE_PG_S) * (self.base_size ** 2.25)
        
        # Surface area scaling (size^2) governing nutrient absorption/consumption rates
        # Reference: Square-cube law; surface area limits exchange with environment
        # Units: pg glucose-equivalent nutrient absorbed per second, later scaled by elapsed dt.
        self.consumption_speed = getattr(species, 'consumption_speed', DEFAULT_CONSUMPTION_SPEED_PG_S) * (self.base_size ** 2)
        
        # Kinematic scaling in low Reynolds number environments (microscale)
        # Reference: Stokes' Law for viscous drag (Fd = 6*pi*eta*r*v). V scales inversely with radius (base_size)
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
        
        # Bacterial run-and-tumble temporal mechanics (E. coli chemotaxis model)
        # Reference: E. coli naturally "run" for ~1 second and "tumble" for ~0.1 seconds to reorient.
        self.base_run_duration = self._sample_run_duration()
        self.base_tumble_duration = BACTERIA_TUMBLE_DURATION_SEC
        self.state_timer = self.base_run_duration
        self.is_tumbling = False
        self._tumble_reoriented = False
        self.last_metabolic_consumption_pg = 0.0

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
    
    def init(self, species, x, y, size=1.0, initial_food=None):
        self.__init__(species, x, y, size, initial_food)
        return self

    def step(self, food_nodes, bounds, dt, viewport_bounds=None, food_hash=None):
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
        # Senescence: Replicative senescence and ROS (Reactive Oxygen Species) accumulation
        # As age increases, cellular damage increases metabolic maintenance overhead.
        # Reference: Free radical theory of aging (Harman, 1956). Models linear increase in basal energy.
        # Average time to double metabolic cost (senescence base lifespan) approx 3 minutes (180s) in sim.
        age_factor = 1.0 + (self.age / (SENESCENCE_BASE_LIFESPAN_SEC * self.lifespan_variance))
        
        # Metabolism: lose nutrient over elapsed real time.
        metabolic_demand_pg = (self.metabolic_rate * dt) * age_factor
        self.last_metabolic_consumption_pg = min(max(0.0, self.food_level), metabolic_demand_pg)
        self.food_level = max(0, self.food_level - metabolic_demand_pg)
        
        # Update size based on fullness
        cap = self.food_capacity if self.food_capacity > 0 else 1.0
        fullness = max(0.0, min(1.0, self.food_level / cap))
        self.size = self.base_size * (BACTERIA_STARVED_RENDER_SCALE + BACTERIA_FULL_RENDER_SCALE_RANGE * fullness)
        
        # Compute effective greediness that rises toward 1 as stomach empties
        base_gx = getattr(self.species, 'greediness', 0.5)
        hunger = 1.0 - fullness
        gx_eff = max(0.0, min(1.0, base_gx + (1.0 - base_gx) * hunger))
        
        # Pathfind with chemotaxis: score food nodes by physical diffusion gradients
        # Reference: Fickian nutrient fields and Weber-Fechner receptor compression.
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
            # Equation: score = greed * log1p(C / C_ref), matching receptor saturation/compression of strong signals.
            score = gx_eff * math.log1p(concentration / CHEMOTAXIS_REFERENCE_PG_PER_UM2)
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

        # Base movement speed
        spd = self.speed * dt

        if target is not None and has_gradient and not self.is_tumbling:
            desired_heading = math.atan2(gradient_y, gradient_x)
            # Equation: signed angular error uses atan2(sin, cos) to wrap to [-pi, pi] without discontinuity.
            heading_error = math.atan2(math.sin(desired_heading - self.heading), math.cos(desired_heading - self.heading))
            # Equation: dtheta <= omega * dt gives smooth low-Reynolds reorientation while translation remains forward-only.
            max_turn = CHEMOTAXIS_STEER_RATE_RAD_S * gx_eff * dt
            self.heading += max(-max_turn, min(max_turn, heading_error))

        # Multi-node consumption: feed from all candidate nodes within eat distance
        max_eat_amt = self.consumption_speed * dt
        amt_eaten_total = 0.0
        is_grazing = False
        
        for node in candidate_nodes:
            if node.amount <= MIN_FOOD_NODE_AMOUNT_PG:
                continue
            # Geometric volumetric scaling: V = 4/3 * pi * r^3 -> r = cbrt(3V / 4pi)
            target_radius = max(FOOD_NODE_MIN_RADIUS_UM, ((3.0 * max(FOOD_NODE_MIN_VOLUME_PG, node.amount)) / (4.0 * math.pi)) ** (1.0 / 3.0))
            eat_dist = target_radius + self.size + FOOD_EATING_DISTANCE_MARGIN_UM
            
            d = math.hypot(node.x - self.x, node.y - self.y)
            
            if d <= eat_dist:
                is_grazing = True
                
            # If close enough to feed (small margin)
            if d <= eat_dist + FOOD_EATING_DISTANCE_GRACE_UM and amt_eaten_total < max_eat_amt:
                eat_amt = min(max_eat_amt - amt_eaten_total, node.amount, max(0.0, self.food_capacity - self.food_level))
                self.food_level += eat_amt
                node.amount -= eat_amt
                amt_eaten_total += eat_amt

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
        
        # Apply Brownian rotational noise (rotational diffusion)
        # Reference: Rotational diffusion coefficient scales roughly ~ 1 / r^3
        # D_rot = kT / (8 pi eta r^3). A standard value for 1 um sphere in water at room temp is ~0.16 rad^2/s.
        rotational_diffusion = BACTERIA_ROTATIONAL_DIFFUSION_COEFFICIENT_RAD2_S / max(0.1, self.base_size ** 3)
        # Brownian angular variance is sigma^2 = 2 D_rot dt, so the sampled angle uses sqrt(2 D_rot dt).
        rotational_sigma_rad = math.sqrt(BROWNIAN_ROTATIONAL_VARIANCE_FACTOR * rotational_diffusion * dt)
        self.heading += random.gauss(0, rotational_sigma_rad)
        
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

    def attempt_fission(self, envf, rep_prob, rep_thresh, dt):
        """
        Evaluate and execute cell division (binary fission).
        Reference: Stouthamer (1973), A theoretical study on the amount of ATP required for synthesis of microbial cell material.
        Cost equivalent of chromosome, ribosomes, and cell wall assembly equates ~40 pg glucose-equivalent nutrient per 1 um^3 of new cell volume.
        Requires sufficient accumulated biomass.
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
                    return Bacterium(self.species, self.x, self.y, size=max(0.1, random.gauss(self.base_size, self.base_size * 0.2)), initial_food=child_food)
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
                 base_size=DEFAULT_BASE_SIZE_UM, sensor_radius=DEFAULT_SENSOR_RADIUS_UM):
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
        # Each bacterium gets a random size
        self.members = [Bacterium(self, random.randint(0,800), random.randint(0,600),
                                  size=max(0.1, random.gauss(self.base_size, self.base_size * 0.2)))
                        for _ in range(init)]
        self.history = []         # track population over time
        self.total_births = 0     # overall births
        self.total_deaths = 0     # overall deaths

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

        ds_cx, ds_cy = getattr(self.sim, 'petri_center', (w/2, h/2))[0] if hasattr(self, 'sim') else w/2, getattr(self.sim, 'petri_center', (w/2, h/2))[1] if hasattr(self, 'sim') else h/2
        ds_r = getattr(self.sim, 'petri_radius', math.hypot(w, h)/2 + PETRI_RADIUS_PADDING_UM) if hasattr(self, 'sim') else math.hypot(w, h)/2 + PETRI_RADIUS_PADDING_UM
        
        for b in list(self.members):
            b.step(food_nodes, bounds=(ds_cx, ds_cy, ds_r), dt=dt, viewport_bounds=viewport_bounds, food_hash=food_hash)
            # Reproduction trial based on the Verhulst (Logistic) growth model modified by environmental fitness.
            # Reference: Verhulst (1838). dp/dt = r*p*(1 - p/K). Here, `capf` represents the limiting (1 - p/K) term.
            rep_prob = (self.base_rep * dt) * capf * envf
            rep_prob = max(0, min(rep_prob, 1))  # Clamp to [0,1]
            
            new_bacterium = b.attempt_fission(envf, rep_prob, self.rep_thresh, dt)
            if new_bacterium:
                births.append(new_bacterium)
            
            # Senescence modeled via Gompertz law of mortality: mu(x) = A * e^(alpha * x) + B
            # Reference: Gompertz (1825). Exponential increase in death rate due to intrinsic cellular senescence.
            # Alpha governs the steepness of the decay curve.
            gompertz_alpha = GOMPERTZ_MORTALITY_ALPHA / max(MIN_LIFESPAN_VARIANCE_DENOMINATOR, b.lifespan_variance)
            age_death_penalty = (self.base_aging * dt) * math.exp(gompertz_alpha * b.age)
            # Cap the penalty to avoid overflow probabilities per simulation tick
            age_death_penalty = min(1.0, age_death_penalty)
            
            # Death trial: extrinsic (environmental toxicity) + intrinsic (Gompertz senescence) + starvation
            toxic_death_prob = ((1.0 - envf) * ENV_TOXIC_DEATH_RATE_PER_SEC) * dt + age_death_penalty
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
    def __init__(self, spp, env, cap, steps=500):
        # give each species a reference to the full list
        for sp in spp:
            setattr(sp, 'species_list', spp)
            setattr(sp, 'sim', self)
        self.species = spp
        self.env = env           # environment dict: oxygen/CO2/food
        self.cap = cap           # carrying capacity
        self.steps = steps       # legacy max update count (not used by the real-time loop)
        self.ecm = ECM(cell_size=ECM_CELL_SIZE_UM) # Inject Extracellular Matrix for molecular diffusion
        # stats collected each step
        self.stats = self._empty_stats()
        self.species_step_stats = []  # List of dicts: [{'births': [...], 'deaths': [...]}]
        self.step_idx = 0        # current time step
        self.running = False     # is simulation active?
        self.food_nodes = []
        self.drop_frequency = env.get('drop_frequency', 5)
        self.drop_range = env.get('drop_range', (2, 8))

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
            'overall_deaths': []
        }

    def reset(self):
        """
        Reset populations to initial values and clear stats.
        """
        for sp in self.species:
            initial = int(sp.history[0] if sp.history and sp.history[0] else len(sp.members))
            sp.members = []
            for _ in range(initial):
                if hasattr(self, 'petri_radius'):
                    r_b = random.random() * self.petri_radius * 0.8
                    th = random.random() * 2 * math.pi
                    x = self.petri_center[0] + r_b * math.cos(th)
                    y = self.petri_center[1] + r_b * math.sin(th)
                else:
                    x = random.randint(0, 800)
                    y = random.randint(0, 600)
                sp.members.append(Bacterium(sp, x, y, size=max(0.5, random.gauss(1.0, 0.2))))
            sp.history = []
            sp.total_births = 0
            sp.total_deaths = 0
        self.stats = self._empty_stats()
        self.species_step_stats = []
        self.step_idx = 0
        self.running = False

    def _process_interactions(self, spatial_hash, CELL_SIZE, dt):
        """
        Process Inter/Intra species interactions and ECM biochemical exposure.
        Reference: Contact-dependent inhibition (T6SS) and diffusible antimicrobial peptides (bacteriocins).
        """
        for sp in self.species:
            bac_prod_rate = (sp.bacteriocin_production * BACTERIOCIN_PRODUCTION_COST_RATE) * dt
            immunity_mitigation = min(BACTERIOCIN_MAX_IMMUNITY_MITIGATION, sp.immunity_factor * BACTERIOCIN_IMMUNITY_SCALING)
            
            for b in sp.members:
                if b.food_level < 0:
                    continue
                cx = int(b.x // CELL_SIZE)
                cy = int(b.y // CELL_SIZE)
                
                b_size_sq_dmg = (sp.t6ss_activity * (b.size ** 2) * T6SS_DAMAGE_SCALING_FACTOR) * dt
                
                # Bacteriocin Emission (metabolic cost):
                # Emission occurs consistently, scaled by bacteriocin_production
                # Reference: Ribosomally synthesized antimicrobial peptides. Production costs ATP/amino acids.
                if bac_prod_rate > 0 and b.food_level > bac_prod_rate:
                    b.food_level -= bac_prod_rate
                    self.ecm.add(b.x, b.y, bac_prod_rate * BACTERIOCIN_ECM_SECRETION_MULTIPLIER, layer='bac')
                
                # Environmental bacteriocin exposure
                local_bac = self.ecm.sample(b.x, b.y, layer='bac')
                if local_bac > 0:
                    # Damage scaled by immunity factor acting as competitive inhibitor.
                    bac_dmg = (local_bac * BACTERIOCIN_DAMAGE_COEFFICIENT) * (1.0 - immunity_mitigation) * dt
                    b.food_level -= bac_dmg

                # Check 9 surrounding cells for physical contact (T6SS / Predation / Sharing)
                for nx in (cx - 1, cx, cx + 1):
                    for ny in (cy - 1, cy, cy + 1):
                        neighbors = spatial_hash.get((nx, ny), [])
                        for other in neighbors:
                            if other is b or other.food_level < 0:
                                continue
                                
                            dx = b.x - other.x
                            dy = b.y - other.y
                            if abs(dx) > BACTERIA_INTERACTION_AXIS_CUTOFF_UM or abs(dy) > BACTERIA_INTERACTION_AXIS_CUTOFF_UM:
                                continue
                                
                            dist_sq = dx**2 + dy**2
                            interact_dist = b.size + other.size + BACTERIA_CONTACT_MARGIN_UM
                            if dist_sq < interact_dist**2:
                                if b.species == other.species:
                                    # Intraspecies cooperation: share food slowly to help starving kin
                                    if b.food_level > b.food_capacity * sp.share_thresh:
                                        diff = b.food_level - other.food_level
                                        if diff > 0:
                                            # Mathematical Rationale: Intraspecies sharing ensures that the recipient's food level
                                            # never exceeds its defined biological food capacity.
                                            transfer_amount = diff * dt
                                            max_transfer = other.food_capacity - other.food_level
                                            actual_transfer = min(transfer_amount, max_transfer)
                                            if actual_transfer > 0:
                                                b.food_level -= actual_transfer
                                                other.food_level += actual_transfer
                                else:
                                    # Interspecies competition
                                    # Type VI Secretion System (T6SS) Contact Dueling:
                                    # Reference: Spring-loaded contact weapon delivering effector toxins directly into neighbor.
                                    # Base damage scales with T6SS activity expression and physical area.
                                    immunity_defense = min(T6SS_MAX_IMMUNITY_DEFENSE, (other.species.immunity_factor * other.size) / (b.size + T6SS_IMMUNITY_SIZE_EPSILON_UM))
                                    actual_dmg = b_size_sq_dmg * (1.0 - immunity_defense)
                                    other.food_level -= actual_dmg

                                    # Phagocytosis (Engulfment): if much larger, instantly consume
                                    # The named size-ratio threshold gives roughly an eightfold volume advantage.
                                    if b.food_level < b.food_capacity * sp.predation_thresh:
                                        if b.size > other.size * PHAGOCYTOSIS_SIZE_RATIO_THRESHOLD:
                                            # Predation equation: digestible gain = cytosolic nutrient + labile dry necromass.
                                            prey_volume_um3 = other.base_size ** 3
                                            prey_necromass_pg = prey_volume_um3 * CELL_WET_MASS_PG_PER_UM3 * CELL_DRY_MASS_FRACTION * LABILE_NECROMASS_FRACTION
                                            digestible_gain = max(0.0, other.food_level) * CYTOSOL_RECOVERY_FRACTION + prey_necromass_pg
                                            b.food_level = min(b.food_capacity, b.food_level + digestible_gain)
                                            other.food_level = CONSUMED_PREY_SENTINEL_FOOD_LEVEL

    def update(self, w, h, dt, viewport_bounds=None):
        """
        Advance one time step if running.
        viewport_bounds: dict with keys 'left', 'right', 'top', 'bottom' for visible area
        """
        if not hasattr(self, 'petri_center'):
            self.petri_center = (w/2, h/2)
            self.petri_radius = math.hypot(w, h) / 2 + PETRI_RADIUS_PADDING_UM

        if not self.running:
            return
            
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
            node = FoodNode.acquire().init(x, y, amount)
            self.food_nodes.append(node)
            
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
        
        # --- Spatial Hash for Bacterial Interactions ---
        CELL_SIZE = BACTERIA_INTERACTION_HASH_CELL_SIZE_UM
        spatial_hash = {}
        for sp in self.species:
            for b in sp.members:
                cx = int(b.x // CELL_SIZE)
                cy = int(b.y // CELL_SIZE)
                if (cx, cy) not in spatial_hash:
                    spatial_hash[(cx, cy)] = []
                spatial_hash[(cx, cy)].append(b)

        # Process Inter/Intra species interactions
        self._process_interactions(spatial_hash, CELL_SIZE, dt)

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
                to_kill = random.sample(all_bacteria, min(excess, len(all_bacteria)))
                for sp, b in to_kill:
                    sp.members.remove(b)
                    sp.total_deaths += 1
                    forced_deaths_per_species[sp] += 1
        for idx, sp in enumerate(self.species):
            deaths_list[idx] += forced_deaths_per_species[sp]
        sum_forced = sum(forced_deaths_per_species.values())
        elapsed_seconds = (self.stats['elapsed_seconds'][-1] if self.stats['elapsed_seconds'] else 0.0) + dt
        step_births = bsum
        step_deaths = dsum + sum_forced
        overall_births = (self.stats['overall_births'][-1] if self.stats['overall_births'] else 0) + step_births
        overall_deaths = (self.stats['overall_deaths'][-1] if self.stats['overall_deaths'] else 0) + step_deaths
        self.stats['time'].append(self.step_idx)
        self.stats['elapsed_seconds'].append(elapsed_seconds)
        self.stats['dt_seconds'].append(dt)
        self.stats['tickrate_hz'].append((1.0 / dt) if dt > 0 else 0.0)
        self.stats['total'].append(sum(len(sp.members) for sp in self.species))
        self.stats['births'].append(step_births)
        self.stats['deaths'].append(step_deaths)
        self.stats['overall_births'].append(overall_births)
        self.stats['overall_deaths'].append(overall_deaths)
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
                # Reference: Stoichiometry mapping cell metabolism to gas usage exactly. 
                # 1 mole Glucose (180.16 g) consumes 6 moles O2 (192.0 g) and produces 6 moles CO2 (264.0 g).
                o2_consumed_pg = total_glucose_pg * (O2_SIX_MOLES_MASS_G / GLUCOSE_MOLAR_MASS_G) * sp.o2_sens
                co2_produced_pg = total_glucose_pg * (CO2_SIX_MOLES_MASS_G / GLUCOSE_MOLAR_MASS_G) * sp.o2_sens
                
                o2_consumption_fraction += o2_consumed_pg / ENV_100_O2_PG
                co2_production_fraction += co2_produced_pg / ENV_100_CO2_PG
                
                # Carbonic acid (H2CO3) weakly dissociates (pKa 6.35). Moderate effect on environmental pH.
                acid_production += (co2_produced_pg / ENV_100_CO2_PG) * CARBONIC_ACID_PH_EFFECT
                
            elif sp.o2_sens < 0:
                # Anaerobic Fermentation (Lactic Acid Fermentation): C6H12O6 -> 2 Lactic Acid
                # Reference: Stoichiometry mapping homolactic fermentation for anaerobic life.
                # 1 mole Glucose (180.16 g) produces 2 moles Lactic Acid (180.16 g, conservation of mass).
                lactic_acid_pg = total_glucose_pg * abs(sp.o2_sens)
                
                # Small trace amounts of CO2 depending on fermentation pathway (e.g. heterolactic models).
                co2_produced_pg = total_glucose_pg * HETEROLACTIC_CO2_TRACE_FRACTION * abs(sp.o2_sens)
                co2_production_fraction += co2_produced_pg / ENV_100_CO2_PG
                
                # Lactic acid strongly drops pH due to highly acidic dissociation equilibrium (pKa 3.86).
                # Reference: Henderson-Hasselbalch representation inside small buffered environments.
                acid_production += (lactic_acid_pg / ENV_100_O2_PG) * LACTIC_ACID_PH_EFFECT

        self.env['o2'] -= o2_consumption_fraction
        self.env['co2'] += co2_production_fraction
        self.env['ph'] -= acid_production

        # Atmospheric gas dissolution and physiological pH buffering mechanics
        # Reference: Fick's first law of diffusion (J = -D * dC/dx). Concentration gradient drives gas toward atmospheric equilibrium.
        # Assuming a shallow 1 mL volume open to sea-level atmospheric pressure, equilibrium typically reaches ~30% per second.
        diffusion_rate = GAS_DIFFUSION_RATE_SEC * dt
        self.env['o2'] += (ATMOSPHERIC_O2_FRACTION - self.env['o2']) * diffusion_rate
        self.env['co2'] += (ATMOSPHERIC_CO2_FRACTION - self.env['co2']) * diffusion_rate
        
        # Buffer system (e.g. bicarbonate buffer system naturally driving pH towards neutral 7 in open systems)
        # Reference: Henderson-Hasselbalch equation and Le Chatelier's equilibrium principle in aqueous environments.
        buffer_rate = BUFFER_RESTORATION_RATE_SEC * dt
        self.env['ph'] += (NEUTRAL_PH - self.env['ph']) * buffer_rate

        # Cap bounds
        self.env['o2'] = max(0.0, min(1.0, self.env['o2']))
        self.env['co2'] = max(0.0, min(1.0, self.env['co2']))
        self.env['ph'] = max(0.0, min(14.0, self.env['ph']))

        self.step_idx += 1

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
            writer.writerow(['time', 'total', 'births', 'deaths', 'elapsed_seconds', 'dt_seconds', 'tickrate_hz', 'overall_births', 'overall_deaths'])
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
                    self.stats['overall_deaths'][i]
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
        for row in lines[stats_start:stats_end]:
            if len(row) >= 4:
                step_time = int(row[header_index.get('time', 0)])
                births = int(row[header_index.get('births', 2)])
                deaths = int(row[header_index.get('deaths', 3)])
                cumulative_births += births
                cumulative_deaths += deaths
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

