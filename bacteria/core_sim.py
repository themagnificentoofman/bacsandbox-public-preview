import math, random
import sys
import os
import csv

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'

from .constants import *

class FoodNode:
    """
    Represents a food node on the grid.
    """
    def __init__(self, x, y, amount):
        self.x = x
        self.y = y
        self.amount = amount

class Bacterium:
    """
    Individual bacterium with position and species reference.
    """
    def __init__(self, species, x, y, size=1.0, initial_food=None):
        self.species = species
        self.x = x
        self.y = y
        self.base_size = size  # Unit: micrometers (um)
        # Volumetric scaling (size^3) assuming constant cellular density context.
        # Unit mappings: 1 distance unit = 1 um. 1 food unit = 1 picogram (pg) of Glucose.
        # Reference: Simple geometric scaling for 3D objects; ~1 pg mass per 1 um^3 volume.
        self.food_capacity = getattr(species, 'food_capacity', 10.0) * (self.base_size ** 3)
        self.food_level = initial_food if initial_food is not None else self.food_capacity
        
        # Kleiber's Law (Allometric scaling): Metabolic rate scales to Mass^(3/4)
        # Reference: Max Kleiber (1932), "Body size and metabolism". If Mass ~ size^3, then Meta ~ size^(2.25)
        # Units: pg of Glucose consumed per tick.
        self.metabolic_rate = getattr(species, 'metabolic_rate', 0.1) * (self.base_size ** 2.25)
        
        # Surface area scaling (size^2) governing nutrient absorption/consumption rates
        # Reference: Square-cube law; surface area limits exchange with environment
        # Units: pg of Glucose absorbed per tick.
        self.consumption_speed = getattr(species, 'consumption_speed', 1.0) * (self.base_size ** 2)
        
        # Kinematic scaling in low Reynolds number environments (microscale)
        # Reference: Stokes' Law for viscous drag (Fd = 6*pi*eta*r*v). V scales inversely with radius (base_size)
        # Units: um/second.
        self.speed = getattr(species, 'speed', 40.0) / max(0.1, self.base_size)
        
        cap = self.food_capacity if self.food_capacity > 0 else 1.0
        self.size = self.base_size * (0.5 + 0.5 * (self.food_level / cap))
        
        self.age = 0.0
        self.lifespan_variance = max(0.4, random.gauss(1.0, 0.2))
        
        # Rotational kinematics in low Reynolds number regime.
        # Torque tau = 8 * pi * eta * r^3 * omega.
        # Assuming flagellar motor output torque scales with cell surface area (r^2), max angular velocity (omega_max) scales as r^2 / r^3 = 1 / r.
        # Base maximum turn rate is modeled as pi radians/second for a 1 um bacterium.
        self.turn_rate = math.pi / max(0.1, self.base_size)
        self.heading = random.uniform(0, 2 * math.pi) if 'math' in globals() else random.random()
        
        # Bacterial run-and-tumble temporal mechanics (E. coli chemotaxis model)
        # Reference: E. coli naturally "run" for ~1 second and "tumble" for ~0.1 seconds to reorient.
        self.base_run_duration = random.gauss(1.0, 0.2)
        self.base_tumble_duration = 0.1
        self.state_timer = self.base_run_duration
        self.is_tumbling = False

    def step(self, food_nodes, grid_min, grid_max):
        self.age += 1.0 / FPS
        # Senescence: Replicative senescence and ROS (Reactive Oxygen Species) accumulation
        # As age increases, cellular damage increases metabolic maintenance overhead.
        # Reference: Free radical theory of aging (Harman, 1956). Models linear increase in basal energy.
        # Average time to double metabolic cost (senescence base lifespan) approx 3 minutes (180s) in sim.
        BASE_LIFESPAN_SEC = 180.0
        age_factor = 1.0 + (self.age / (BASE_LIFESPAN_SEC * self.lifespan_variance))
        
        # Metabolism: lose food each tick (scaled by FPS)
        self.food_level = max(0, self.food_level - (self.metabolic_rate / FPS) * age_factor)
        
        # Update size based on fullness
        cap = self.food_capacity if self.food_capacity > 0 else 1.0
        fullness = max(0.0, min(1.0, self.food_level / cap))
        self.size = self.base_size * (0.5 + 0.5 * fullness)
        
        # Compute effective greediness that rises toward 1 as stomach empties
        base_gx = getattr(self.species, 'greediness', 0.5)
        hunger = 1.0 - fullness
        gx_eff = max(0.0, min(1.0, base_gx + (1.0 - base_gx) * hunger))
        
        # Pathfind with chemotaxis: score food nodes by physical diffusion gradients
        # Reference: Fick's laws of diffusion model concentration dropping off inversely with distance or distance squared.
        # Chemotactic receptor response scales logarithmically with this concentration.
        target = None
        best_score = -float('inf')
        closest_node = None
        min_dist = float('inf')
        max_range = self.species.sensor_radius
        for node in food_nodes:
            dx = node.x - self.x
            dy = node.y - self.y
            dist = (dx*dx + dy*dy) ** 0.5
            if dist < min_dist:
                min_dist = dist
                closest_node = node
            if dist > max_range:
                continue
            
            # Chemotactic score: Inverse square law for concentration model at a distance.
            # Base concentration is proportional to node amount.
            concentration = node.amount / max(1.0, dist**2)
            # Receptor activation incorporates greediness (ignores weak signals if not greedy)
            score = concentration * gx_eff
            
            if score > best_score:
                best_score = score
                target = node
                
        # State management: Run vs Tumble
        self.state_timer -= 1.0 / FPS
        if self.state_timer <= 0.0:
            self.is_tumbling = not self.is_tumbling
            if self.is_tumbling:
                self.state_timer = self.base_tumble_duration
            else:
                self.state_timer = self.base_run_duration * random.gauss(1.0, 0.2)

        # Base movement speed
        spd = self.speed / FPS

        # Movement: steer toward target with physical rotation limits
        if target:
            # Geometric volumetric scaling: V = 4/3 * pi * r^3 -> r = cbrt(3V / 4pi)
            # where V = target.amount (pg = um^3)
            target_radius = max(0.5, ((3.0 * max(0.1, target.amount)) / (4.0 * math.pi)) ** (1.0 / 3.0))
            eat_dist = target_radius + self.size + 1.0 # Stop at the edge
            
            d = ((self.x - target.x)**2 + (self.y - target.y)**2) ** 0.5
            
            if d > eat_dist:
                target_angle = math.atan2(target.y - self.y, target.x - self.x)
                angle_diff = (target_angle - self.heading + math.pi) % (2 * math.pi) - math.pi
                max_turn = self.turn_rate / FPS
                self.heading += max(-max_turn, min(angle_diff, max_turn))
            else:
                # Motility inhibition down-regulation upon reaching concentrated nutrient sources.
                # Reference: Bacterial cells often halt flagellar rotation or utilize pili to attach to surfaces (biofilm initiation) when grazing.
                target_angle = math.atan2(target.y - self.y, target.x - self.x)
                self.heading = target_angle
                spd = 0.0
                
            # If close enough to feed (small margin)
            if d <= eat_dist + 3.0:
                eat_amt = min(self.consumption_speed / FPS, target.amount, max(0.0, self.food_capacity - self.food_level))
                self.food_level += eat_amt
                target.amount -= eat_amt
        else:
            if self.is_tumbling:
                # During tumble, motility is spent entirely on reorientation
                # Model tumble as applying random torque, capped by maximum turn rate
                self.heading += random.uniform(-self.turn_rate, self.turn_rate) / FPS
                spd *= 0.1 # Very slow translation during tumble
            else:
                # Run phase: try to follow closest node if exists
                if closest_node:
                    target_angle = math.atan2(closest_node.y - self.y, closest_node.x - self.x)
                    angle_diff = (target_angle - self.heading + math.pi) % (2 * math.pi) - math.pi
                    max_turn = self.turn_rate / FPS
                    self.heading += max(-max_turn, min(angle_diff, max_turn))
        
        # Apply Brownian rotational noise (rotational diffusion)
        # Reference: Rotational diffusion coefficient scales roughly ~ 1 / r^3
        # D_rot = kT / (8 pi eta r^3). A standard value for 1 um sphere in water at room temp is ~0.16 rad^2/s.
        base_diffusion_coefficient = 0.16
        rotational_diffusion = base_diffusion_coefficient / max(0.1, self.base_size ** 3)
        self.heading += random.gauss(0, rotational_diffusion / FPS)
        
        # Forward translation
        self.x += math.cos(self.heading) * spd
        self.y += math.sin(self.heading) * spd
            
        # Clamp position
        self.x = min(max(self.x, grid_min), grid_max)
        self.y = min(max(self.y, grid_min), grid_max)

class BacteriaSpecies:
    """
    Represents one species, tracking its parameters,
    individual members, and population history.
    """
    def __init__(self, name, init, rep, aging, so, sf, speed, ph, color,
                 food_capacity=10.0, metabolic_rate=0.1, consumption_speed=1.0, greediness=0.5,
                 aggression=0.5, defense=0.5, rep_thresh=0.8, share_thresh=0.6,
                 predation_thresh=0.9, base_size=1.0, sensor_radius=200.0):
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
        self.aggression = max(0.0, min(float(aggression), 1.0))
        self.defense = max(0.0, min(float(defense), 1.0))
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

    def step(self, env, cap, w, h, food_nodes):
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
        curr_ph = env.get('ph', 7.0)
        ph_diff = curr_ph - self.preferred_ph if self.preferred_ph is not None else 0
        ph_tolerance = 1.5  # Standard deviation of pH tolerance
        ph_suit = math.exp(-0.5 * (ph_diff / ph_tolerance)**2)

        # O2/CO2 suitability
        o2_level = env.get('o2', 0.21)
        if self.o2_sens > 0:
            # Aerobes: Need O2 to survive. Scale up to 10% O2
            o2_suit = min(1.0, o2_level / 0.10) * self.o2_sens + (1 - self.o2_sens)
        elif self.o2_sens < 0:
            # Anaerobes: Poisoned by O2. Toxicity kicks in heavily over 2% O2
            toxicity = min(1.0, o2_level / 0.02) * abs(self.o2_sens)
            o2_suit = 1.0 - toxicity
        else:
            o2_suit = 1.0  # Indifferent

        # Overall environmental fitness [0, 1]
        envf = ph_suit * max(0.0, o2_suit)

        grid_min = -32768
        grid_max = 32768
        for b in list(self.members):
            b.step(food_nodes, grid_min, grid_max)
            # Reproduction trial based on the Verhulst (Logistic) growth model modified by environmental fitness.
            # Reference: Verhulst (1838). dp/dt = r*p*(1 - p/K). Here, `capf` represents the limiting (1 - p/K) term.
            rep_prob = (self.base_rep / FPS) * capf * envf
            rep_prob = max(0, min(rep_prob, 1))  # Clamp to [0,1]
            
            # Binary fission occurs only when sufficient biomass is accumulated.
            # Reference: Cell size control and division checkpoints in bacteria (adder/sizer models).
            if rep_prob > 0 and random.random() < rep_prob and b.food_level > b.food_capacity * self.rep_thresh:
                # Thermodynamic loss splits CURRENT food pool to avoid instant death on low thresholds.
                # E.g., parent keeps 50% of its current food, child gets 50%
                parent_cost = b.food_level * 0.5
                child_food = b.food_level * 0.5
                b.food_level -= parent_cost
                births.append(Bacterium(self, b.x, b.y, size=max(0.1, random.gauss(self.base_size, self.base_size * 0.2)), initial_food=child_food))
            
            # Senescence modeled via Gompertz law of mortality: mu(x) = A * e^(alpha * x) + B
            # Reference: Gompertz (1825). Exponential increase in death rate due to intrinsic cellular senescence.
            # Alpha governs the steepness of the decay curve. An alpha of 0.04 represents mortality doubling approx every 17s.
            gompertz_alpha = 0.04 / max(0.001, b.lifespan_variance)
            age_death_penalty = (self.base_aging / FPS) * math.exp(gompertz_alpha * b.age)
            # Cap the penalty to avoid overflow probabilities per simulation tick
            age_death_penalty = min(1.0, age_death_penalty)
            
            # Death trial: extrinsic (environmental toxicity) + intrinsic (Gompertz senescence) + starvation
            toxic_death_prob = ((1.0 - envf) * 0.05) / FPS + age_death_penalty
            if b.food_level <= 0 or random.random() < toxic_death_prob:
                deaths.append(b)
                # Decomposition: drop food node (scavenging)
                food_value = b.food_level + b.food_capacity * 0.2 # remaining food + 20% base chassis mass
                if food_value > 0:
                    food_nodes.append(FoodNode(b.x, b.y, food_value))
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

class Simulation:
    """
    Manages multiple species, environment, and global stats.
    """
    def __init__(self, spp, env, cap, steps=500):
        # give each species a reference to the full list
        for sp in spp:
            setattr(sp, 'species_list', spp)
        self.species = spp
        self.env = env           # environment dict: oxygen/CO2/food
        self.cap = cap           # carrying capacity
        self.steps = steps       # max time steps (NOT USED ANYMORE)
        # stats collected each step
        self.stats = {'time': [], 'total': [], 'births': [], 'deaths': []}
        self.species_step_stats = []  # List of dicts: [{'births': [...], 'deaths': [...]}]
        self.step_idx = 0        # current time step
        self.running = False     # is simulation active?
        self.food_nodes = []
        self.drop_frequency = env.get('drop_frequency', 5)
        self.drop_range = env.get('drop_range', (2, 8))

    def reset(self):
        """
        Reset populations to initial values and clear stats.
        """
        for sp in self.species:
            initial = int(sp.history[0] if sp.history and sp.history[0] else len(sp.members))
            sp.members = [Bacterium(sp, random.randint(0,800), random.randint(0,600),
                                    size=max(0.5, random.gauss(1.0, 0.2)))
                          for _ in range(initial)]
            sp.history = []
            sp.total_births = 0
            sp.total_deaths = 0
        self.stats = {k: [] for k in self.stats}
        self.species_step_stats = []
        self.step_idx = 0
        self.running = False

    def update(self, w, h):
        """
        Advance one time step if running:
        - Let each species step
        - Enforce hard carrying capacity
        - Record stats
        - Update environment (deplete resources)
        """
        if not self.running:
            return
        # Drop new food nodes based on percentage chance (scaled by FPS to maintain constant drops/sec)
        if random.random() * 100 < (self.drop_frequency * 144.0 / FPS):
            x = random.randint(0, w)
            y = random.randint(0, h)
            amount = random.uniform(*self.drop_range)
            self.food_nodes.append(FoodNode(x, y, amount))
        # Remove empty food nodes
        self.food_nodes = [node for node in self.food_nodes if node.amount > 0]
        # --- Optimization: Cap food nodes to most recent 1000 ---
        MAX_FOOD_NODES = 1000
        if len(self.food_nodes) > MAX_FOOD_NODES:
            self.food_nodes = self.food_nodes[-MAX_FOOD_NODES:]
        
        # --- Spatial Hash for Bacterial Interactions ---
        CELL_SIZE = 20
        spatial_hash = {}
        for sp in self.species:
            for b in sp.members:
                cx = int(b.x // CELL_SIZE)
                cy = int(b.y // CELL_SIZE)
                if (cx, cy) not in spatial_hash:
                    spatial_hash[(cx, cy)] = []
                spatial_hash[(cx, cy)].append(b)

        # Process Inter/Intra species interactions
        for sp in self.species:
            for b in sp.members:
                cx = int(b.x // CELL_SIZE)
                cy = int(b.y // CELL_SIZE)
                # Check 9 surrounding cells
                for nx in (cx - 1, cx, cx + 1):
                    for ny in (cy - 1, cy, cy + 1):
                        neighbors = spatial_hash.get((nx, ny), [])
                        for other in neighbors:
                            if other is b:
                                continue
                            dist_sq = (b.x - other.x)**2 + (b.y - other.y)**2
                            interact_dist = b.size + other.size + 2.0
                            if dist_sq < interact_dist**2:
                                if b.species == other.species:
                                    # Intraspecies cooperation: share food slowly to help starving kin
                                    # ONLY share if we have plenty of food ourselves (>share_thresh of capacity)
                                    if b.food_level > b.food_capacity * sp.share_thresh:
                                        diff = b.food_level - other.food_level
                                        if diff > 0:
                                            # Balanced sharing rate: transfer 100% of difference over 1 full second
                                            transfer = diff / FPS
                                            b.food_level -= transfer
                                            other.food_level += transfer
                                else:
                                    # Interspecies competition: damage via aggression/defense
                                    # Only attack if we are below our predation threshold capacity
                                    if b.food_level < b.food_capacity * sp.predation_thresh:
                                        # Phagocytosis (Engulfment): if much larger, instantly consume
                                        # We use 2.0x size advantage (which means 8x volume advantage)
                                        if b.size > other.size * 2.0:
                                            # Instantly kill and eat their food pool (plus bit of biomass)
                                            b.food_level = min(b.food_capacity, b.food_level + other.food_level + other.food_capacity * 0.2)
                                            other.food_level = -9999.0
                                        else:
                                            # Base damage scales with surface area (size^2) and aggression
                                            dmg = (sp.aggression * (b.size ** 2) * 5.0) / FPS
                                            # Fractional defense mitigation so swarms can eventually whittle down a giant
                                            mitigation = min(0.95, (other.species.defense * other.size) / (b.size + 0.01))
                                            actual_dmg = dmg * (1.0 - mitigation)
                                            other.food_level -= actual_dmg

        bsum, dsum = 0, 0
        births_list, deaths_list = [], []
        for sp in self.species:
            b, d = sp.step(self.env, self.cap, w, h, self.food_nodes)
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
        self.stats['time'].append(self.step_idx)
        self.stats['total'].append(sum(len(sp.members) for sp in self.species))
        self.stats['births'].append(bsum)
        sum_forced = sum(forced_deaths_per_species.values())
        self.stats['deaths'].append(dsum + sum_forced)
        self.species_step_stats.append({'births': births_list, 'deaths': deaths_list})
        for sp in self.species:
            if sp.history:
                sp.history[-1] = len(sp.members)
            else:
                sp.history.append(len(sp.members))
                
        # Real-world Environmental Exchange & Metabolism
        # Unit Definitions Context: Let 1 food unit map to 1 picogram (pg) of Glucose.
        # Environment is modeled as 1 mL of fluid (10^12 um^3).
        # A 0.21 partial volume (percentage) of O2 dissolved in 1 mL water represents ~8.4 ug of O2 (8.4e6 pg).
        # `env['o2']` of 1.0 represents ~40 mg/L dissolved O2 (4e7 pg/mL limit).
        # `env['co2']` of 1.0 represents ~1.5 g/L dissolved CO2 (1.5e9 pg/mL limit).
        ENV_100_O2_PG = 4.0e7
        ENV_100_CO2_PG = 1.5e9
        
        o2_consumption_fraction = 0.0
        co2_production_fraction = 0.0
        acid_production = 0.0
        
        for sp in self.species:
            pop = len(sp.members)
            if pop == 0: continue
            
            # total_glucose_pg: Total mass of glucose metabolized by this species this tick
            meta_pg = sp.metabolic_rate / FPS
            total_glucose_pg = meta_pg * pop
            
            if sp.o2_sens > 0:
                # Aerobic Respiration: C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O
                # Reference: Stoichiometry mapping cell metabolism to gas usage exactly. 
                # 1 mole Glucose (180.16 g) consumes 6 moles O2 (192.0 g) and produces 6 moles CO2 (264.0 g).
                o2_consumed_pg = total_glucose_pg * (192.0 / 180.16) * sp.o2_sens
                co2_produced_pg = total_glucose_pg * (264.0 / 180.16) * sp.o2_sens
                
                o2_consumption_fraction += o2_consumed_pg / ENV_100_O2_PG
                co2_production_fraction += co2_produced_pg / ENV_100_CO2_PG
                
                # Carbonic acid (H2CO3) weakly dissociates (pKa 6.35). Moderate effect on environmental pH.
                acid_production += (co2_produced_pg / ENV_100_CO2_PG) * 0.02
                
            elif sp.o2_sens < 0:
                # Anaerobic Fermentation (Lactic Acid Fermentation): C6H12O6 -> 2 Lactic Acid
                # Reference: Stoichiometry mapping homolactic fermentation for anaerobic life.
                # 1 mole Glucose (180.16 g) produces 2 moles Lactic Acid (180.16 g, conservation of mass).
                lactic_acid_pg = total_glucose_pg * abs(sp.o2_sens)
                
                # Small trace amounts of CO2 depending on fermentation pathway (e.g. heterolactic models).
                co2_produced_pg = total_glucose_pg * 0.1 * abs(sp.o2_sens)
                co2_production_fraction += co2_produced_pg / ENV_100_CO2_PG
                
                # Lactic acid strongly drops pH due to highly acidic dissociation equilibrium (pKa 3.86).
                # Reference: Henderson-Hasselbalch representation inside small buffered environments.
                acid_production += (lactic_acid_pg / ENV_100_O2_PG) * 2.5

        self.env['o2'] -= o2_consumption_fraction
        self.env['co2'] += co2_production_fraction
        self.env['ph'] -= acid_production

        # Atmospheric gas dissolution and physiological pH buffering mechanics
        # Reference: Fick's first law of diffusion (J = -D * dC/dx). Concentration gradient drives gas toward equilibrium (0.21 O2, 0.0 CO2).
        # Assuming a shallow 1 mL volume open to sea-level atmospheric pressure, equilibrium typically reaches ~30% per second.
        # Gas dissolution rate modeled biologically as 0.288 per second.
        GAS_DIFFUSION_RATE_SEC = 0.288
        diffusion_rate = GAS_DIFFUSION_RATE_SEC / FPS
        self.env['o2'] += (0.21 - self.env['o2']) * diffusion_rate
        self.env['co2'] += (0.00 - self.env['co2']) * diffusion_rate
        
        # Buffer system (e.g. bicarbonate buffer system naturally driving pH towards neutral 7 in open systems)
        # Reference: Henderson-Hasselbalch equation and Le Chatelier's equilibrium principle in aqueous environments.
        # Approximated restoration rate of 0.144 pH units towards neutral per second.
        BUFFER_RESTORATION_RATE_SEC = 0.144
        buffer_rate = BUFFER_RESTORATION_RATE_SEC / FPS
        self.env['ph'] += (7.0 - self.env['ph']) * buffer_rate

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
            writer.writerow(['time', 'total', 'births', 'deaths'])
            for i in range(len(self.stats['time'])):
                writer.writerow([
                    self.stats['time'][i],
                    self.stats['total'][i],
                    self.stats['births'][i],
                    self.stats['deaths'][i]
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
        # Parse global stats
        stats_start = 1
        for i, row in enumerate(lines):
            if row and row[0] == 'time' and len(row) >= 4 and row[1] == 'total':
                stats_start = i + 1
                break
        
        stats_end = stats_start
        while stats_end < len(lines) and lines[stats_end]:
            stats_end += 1
        self.stats = {'time': [], 'total': [], 'births': [], 'deaths': []}
        for row in lines[stats_start:stats_end]:
            if len(row) >= 4:
                self.stats['time'].append(int(row[0]))
                self.stats['total'].append(int(row[1]))
                self.stats['births'].append(int(row[2]))
                self.stats['deaths'].append(int(row[3]))
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

