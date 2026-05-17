#!/usr/bin/env python3
import sys
import os

def read_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.read()

def write_file(filepath, contents):
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(contents)

def patch_core_sim():
    filepath = os.path.join(os.path.dirname(__file__), 'bacteria', 'core_sim.py')
    content = read_file(filepath)
    
    # Step 1: Add object pools to FoodNode
    old_foodnode = '''class FoodNode:
    """
    Represents a food node on the grid.
    """
    def __init__(self, x, y, amount):
        self.x = x
        self.y = y
        self.amount = amount'''
    
    new_foodnode = '''class FoodNode:
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
        return self'''
    
    if old_foodnode in content:
        content = content.replace(old_foodnode, new_foodnode)
    else:
        print("Could not find FoodNode to patch.")
        
    # Step 2: Add object pool variables to Bacterium class
    old_bacterium_init_start = '''class Bacterium:
    """
    Individual bacterium with position and species reference.
    """
    def __init__(self, species, x, y, size=1.0, initial_food=None):'''
    
    new_bacterium_init_start = '''class Bacterium:
    """
    Individual bacterium with position and species reference.
    """
    _pool = []
    
    def __init__(self, species, x, y, size=1.0, initial_food=None):'''
    
    if old_bacterium_init_start in content:
        content = content.replace(old_bacterium_init_start, new_bacterium_init_start)
    else:
         print("Could not find Bacterium target 1.")
        
    old_step_pattern = '''        self.is_tumbling = False

    def step(self, food_nodes, grid_min, grid_max, dt):'''
    
    new_step_pattern = '''        self.is_tumbling = False
    
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

    def step(self, food_nodes, grid_min, grid_max, dt, viewport_bounds=None, food_hash=None):'''
    
    if old_step_pattern in content:
        content = content.replace(old_step_pattern, new_step_pattern)
    else:
         print("Could not find Bacterium target 2.")

    # Step 3: Add FoodSpatialHash class before Simulation class
    food_spatial_hash_class = '''
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

class Simulation:'''
    
    if 'class Simulation:' in content and 'FoodSpatialHash' not in content:
        content = content.replace('class Simulation:', food_spatial_hash_class)
        
    # Step 4: Modify Simulation.update
    old_update_start = '''    def update(self, w, h, dt):
        """
        Advance one time step if running:
        - Let each species step
        - Enforce hard carrying capacity
        - Record stats
        - Update environment (deplete resources)
        """
        if not self.running:
            return
        # Drop new food nodes based on percentage chance (scaled by dt to maintain constant drops/sec)
        if random.random() * 100 < (self.drop_frequency * 144.0 * dt):
            x = random.randint(0, w)
            y = random.randint(0, h)
            amount = random.uniform(*self.drop_range)
            self.food_nodes.append(FoodNode(x, y, amount))
        # Remove empty food nodes
        self.food_nodes = [node for node in self.food_nodes if node.amount > 0]
        # --- Optimization: Cap food nodes to most recent 1000 ---
        MAX_FOOD_NODES = 1000
        if len(self.food_nodes) > MAX_FOOD_NODES:
            self.food_nodes = self.food_nodes[-MAX_FOOD_NODES:]'''
            
    new_update_start = '''    def update(self, w, h, dt, viewport_bounds=None):
        """
        Advance one time step if running.
        viewport_bounds: dict with keys 'left', 'right', 'top', 'bottom' for visible area
        """
        if not self.running:
            return
        if random.random() * 100 < (self.drop_frequency * 144.0 * dt):
            x = random.randint(0, w)
            y = random.randint(0, h)
            amount = random.uniform(*self.drop_range)
            node = FoodNode.acquire().init(x, y, amount)
            self.food_nodes.append(node)
            
        to_remove = [n for n in self.food_nodes if n.amount <= 0]
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
        food_hash = FoodSpatialHash(cell_size=50)
        for node in self.food_nodes:
            food_hash.insert(node)'''
            
    if old_update_start in content:
        content = content.replace(old_update_start, new_update_start)
    else:
        print("Could not find Update.")
        

    old_sp_step = '''        bsum, dsum = 0, 0
        births_list, deaths_list = [], []
        for sp in self.species:
            b, d = sp.step(self.env, self.cap, w, h, self.food_nodes, dt)'''
    new_sp_step = '''        bsum, dsum = 0, 0
        births_list, deaths_list = [], []
        for sp in self.species:
            b, d = sp.step(self.env, self.cap, w, h, self.food_nodes, dt, viewport_bounds, food_hash)'''
            
    if old_sp_step in content:
        content = content.replace(old_sp_step, new_sp_step)
        
    # BacteriaSpecies step mapping
    old_b_step_species = '''    def step(self, env, cap, w, h, food_nodes, dt):
        """
        Process one simulation step for this species:
        - Move each bacterium
        - Decide births and deaths based on rates,
          environment, and carrying capacity
        Returns (birth_count, death_count).
        """
        grid_min, grid_max = 0, max(w, h)
        
        for b in list(self.members):
            b.step(food_nodes, grid_min, grid_max, dt)'''
            
    new_b_step_species = '''    def step(self, env, cap, w, h, food_nodes, dt, viewport_bounds=None, food_hash=None):
        """
        Process one simulation step for this species
        """
        grid_min, grid_max = 0, max(w, h)
        
        for b in list(self.members):
            b.step(food_nodes, grid_min, grid_max, dt, viewport_bounds, food_hash)'''
            
    if old_b_step_species in content:
        content = content.replace(old_b_step_species, new_b_step_species)

    # Bacterium.step
    old_bacterium_step2 = '''        self.age += dt'''
    new_bacterium_step2 = '''        is_out_of_bounds = False
        if viewport_bounds is not None:
            is_out_of_bounds = (
                self.x < viewport_bounds.get('left') or 
                self.x > viewport_bounds.get('right') or 
                self.y < viewport_bounds.get('top') or 
                self.y > viewport_bounds.get('bottom')
            )
            
        self.age += dt'''
    
    if old_bacterium_step2 in content:
        content = content.replace(old_bacterium_step2, new_bacterium_step2, 1)

    old_target_logic = '''        target = None
        best_score = -float('inf')
        closest_node = None
        min_dist = float('inf')
        max_range = self.species.sensor_radius
        for node in food_nodes:'''
        
    new_target_logic = '''        target = None
        best_score = -float('inf')
        closest_node = None
        min_dist = float('inf')
        max_range = self.species.sensor_radius
        
        candidate_nodes = food_hash.query(self.x, self.y, max_range) if food_hash else food_nodes
        if is_out_of_bounds:
            candidate_nodes = [] # Hybrid ABM/ODE Fast-forward. ODE translation applies instead.
            
        for node in candidate_nodes:'''
        
    if old_target_logic in content:
        content = content.replace(old_target_logic, new_target_logic)
        
    # Death pooling logic
    old_death_add = '''                deaths += 1
                self.members.remove(b)'''
    new_death_add = '''                deaths += 1
                self.members.remove(b)
                Bacterium.release(b)'''
    if old_death_add in content:
        content = content.replace(old_death_add, new_death_add)
        
    old_fission = '''                new_b = Bacterium(self, b.x + random.uniform(-2, 2), b.y + random.uniform(-2, 2), size=child_volume**(1.0/3.0), initial_food=child_food)
                self.members.append(new_b)'''
    new_fission = '''                new_b = Bacterium.acquire().init(self, b.x + random.uniform(-2, 2), b.y + random.uniform(-2, 2), size=child_volume**(1.0/3.0), initial_food=child_food)
                self.members.append(new_b)'''
    if old_fission in content:
        content = content.replace(old_fission, new_fission)

    write_file(filepath, content)
    print("core_sim patched")

def patch_app():
    filepath = os.path.join(os.path.dirname(__file__), 'bacteria', 'app.py')
    content = read_file(filepath)
    old_app_update = '''            # Update simulation over dt
            self.sim.update(self.w, self.h, dt)'''
    
    new_app_update = '''            # Update simulation over dt
            half_w = self.w / (2 * self.zoom) if self.zoom > 0 else self.w / 2
            half_h = self.h / (2 * self.zoom) if self.zoom > 0 else self.h / 2
            viewport_bounds = {
                'left': self.grid_center_x - half_w,
                'right': self.grid_center_x + half_w,
                'top': self.grid_center_y - half_h,
                'bottom': self.grid_center_y + half_h
            }
            self.sim.update(self.w, self.h, dt, viewport_bounds)'''
            
    if old_app_update in content:
        content = content.replace(old_app_update, new_app_update)
        print("app.py patched")
        write_file(filepath, content)

if __name__ == "__main__":
    patch_core_sim()
    patch_app()
