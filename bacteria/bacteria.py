import pygame
import sys
import random
import csv
import os
import re
import json
import math
from pygame.locals import *

# ------------------ Constants ------------------
# Frames per second for the simulation
FPS = 144
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

# ------------------ Data Manager ------------------
class SpeciesDataManager:
    """
    Handles loading, saving, and managing species parameters
    stored in a CSV file.
    """
    def __init__(self, filename=DATA_FILE):
        self.filename = filename
        # Updated CSV columns: name, initial count, rates, O2 sens, food sens, speed, preferred pH, food system
        self.fieldnames = ['name', 'init', 'rep', 'death', 'so', 'sf', 'speed', 'ph', 'food_capacity', 'metabolic_rate', 'consumption_speed', 'greediness']
        self.species = []
        self.load()

    def load(self):
        """
        Read species parameters from CSV if it exists,
        otherwise initialize empty list.
        """
        if os.path.exists(self.filename):
            with open(self.filename, newline='') as f:
                self.species = list(csv.DictReader(f))
        else:
            self.species = []

    def save(self):
        """
        Write current species list to CSV file.
        """
        with open(self.filename, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=self.fieldnames)
            writer.writeheader()
            writer.writerows(self.species)

    def add_blank(self):
        """
        Add a new blank species entry with defaults of '0' and food system defaults.
        """
        blank = {k: '0' for k in self.fieldnames}
        blank['name'] = ''  # name left empty for user to fill
        blank['food_capacity'] = '10'
        blank['metabolic_rate'] = '0.1'
        blank['consumption_speed'] = '1.0'
        blank['greediness'] = '0.5'
        self.species.append(blank)
        self.save()

    def remove(self, index):
        """
        Remove species at given index (if valid) and save.
        """
        if 0 <= index < len(self.species):
            self.species.pop(index)
            self.save()

# ------------------ Simulation Entities ------------------
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
    def __init__(self, species, x, y, size=1.0):
        self.species = species
        self.x = x
        self.y = y
        self.size = size
        self.food_capacity = getattr(species, 'food_capacity', 10.0)
        self.food_level = self.food_capacity
        self.metabolic_rate = getattr(species, 'metabolic_rate', 0.1)
        self.consumption_speed = getattr(species, 'consumption_speed', 1.0)
        # Wander state
        self.wander_dir = random.uniform(0, 2*math.pi) if 'math' in globals() else random.random()
        self.wander_change_timer = random.randint(15, 90)

    def move(self, grid_min, grid_max):
        """
        Move by -speed to +speed pixels per frame, scaled by FPS
        """
        spd = self.species.speed / FPS
        dx = random.uniform(-spd, spd)
        dy = random.uniform(-spd, spd)
        self.x = min(max(self.x + dx, grid_min), grid_max)
        self.y = min(max(self.y + dy, grid_min), grid_max)

    def move_towards(self, target_x, target_y, speed):
        dx = target_x - self.x
        dy = target_y - self.y
        dist = (dx**2 + dy**2) ** 0.5
        if dist > 0:
            # Add freedom: move partially toward target plus some noise
            step = min(speed / FPS, dist)
            # Steering factor reduces perfect tweening
            steer = 0.7
            nx = dx / dist
            ny = dy / dist
            # Small random jitter
            jx = random.uniform(-1, 1) * 0.3
            jy = random.uniform(-1, 1) * 0.3
            vx = (nx * steer + jx * (1 - steer)) * step
            vy = (ny * steer + jy * (1 - steer)) * step
            self.x += vx
            self.y += vy

    def step(self, food_nodes, grid_min, grid_max):
        # Metabolism: lose food each tick
        self.food_level = max(0, self.food_level - self.metabolic_rate)
        # Compute effective greediness that rises toward 1 as stomach empties
        base_gx = getattr(self.species, 'greediness', 0.5)
        cap = self.food_capacity if self.food_capacity > 0 else 1.0
        fullness = max(0.0, min(1.0, self.food_level / cap))
        hunger = 1.0 - fullness
        gx_eff = max(0.0, min(1.0, base_gx + (1.0 - base_gx) * hunger))
        # Wander timer update
        self.wander_change_timer -= 1
        if self.wander_change_timer <= 0:
            self.wander_change_timer = random.randint(30, 120)
            self.wander_dir += random.uniform(-0.5, 0.5)
        # Pathfind with greediness: score food nodes by amount and distance
        target = None
        best_score = -float('inf')
        max_range = 200 + 800 * gx_eff  # willing travel distance scales with effective greediness
        for node in food_nodes:
            dx = node.x - self.x
            dy = node.y - self.y
            dist = (dx*dx + dy*dy) ** 0.5
            if dist > max_range:
                continue
            # Score function: amount weighted minus distance penalty
            amount_weight = 1.0 + 4.0 * gx_eff  # up to 5x emphasis on rich nodes
            distance_penalty = (1.0 - gx_eff)  # greedy species care less about distance
            score = node.amount * amount_weight - dist * distance_penalty
            if score > best_score:
                best_score = score
                target = node
        # Movement: if we have a target, steer toward it with noise; otherwise wander
        if target:
            self.move_towards(target.x, target.y, self.species.speed)
            # If close enough, eat
            d = ((self.x - target.x)**2 + (self.y - target.y)**2) ** 0.5
            if d < 6:
                eat_amt = min(self.consumption_speed, target.amount, self.food_capacity - self.food_level)
                self.food_level += eat_amt
                target.amount -= eat_amt
        else:
            # Wander: forward step plus noise
            spd = self.species.speed / FPS
            self.x += math.cos(self.wander_dir) * spd + random.uniform(-0.5, 0.5)
            self.y += math.sin(self.wander_dir) * spd + random.uniform(-0.5, 0.5)
        # Clamp position
        self.x = min(max(self.x, grid_min), grid_max)
        self.y = min(max(self.y, grid_min), grid_max)

class BacteriaSpecies:
    """
    Represents one species, tracking its parameters,
    individual members, and population history.
    """
    def __init__(self, name, init, rep, death, so, sf, speed, ph, color,
                 food_capacity=10.0, metabolic_rate=0.1, consumption_speed=1.0, greediness=0.5):
        # Store parameters
        self.name = name
        self.base_rep = rep       # base reproduction probability
        self.base_death = death   # base death probability
        self.o2_sens = so  # -1 to 1
        self.food_sens = sf
        self.speed = speed  # pixels/sec
        self.preferred_ph = ph
        self.color = color        # draw color for this species
        self.food_capacity = food_capacity
        self.metabolic_rate = metabolic_rate
        self.consumption_speed = consumption_speed
        # New attribute: greediness in [0,1]
        # 0 = never prioritize far rich nodes; 1 = strongly prefer highest food even if far
        self.greediness = max(0.0, min(float(greediness), 1.0))
        # Each bacterium gets a random size (mean=1.0, std=0.2)
        self.members = [Bacterium(self, random.randint(0,800), random.randint(0,600),
                                  size=max(0.5, random.gauss(1.0, 0.2)))
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

        # pH suitability: 1.0 if preferred, drops off as difference
        ph_diff = abs(env['ph'] - self.preferred_ph) if env and 'ph' in env and self.preferred_ph is not None else 0
        ph_suit = max(0, 1 - ph_diff / 7.0) if ph_diff is not None else 1.0

        # O2/CO2 sensitivity
        o2_sens_abs = abs(self.o2_sens)
        denom = o2_sens_abs + 1
        if denom == 0:
            denom = 1
        if self.o2_sens >= 0:
            o2_suit = env['oxygen'] * self.o2_sens if env and 'oxygen' in env else 0
        else:
            o2_suit = (1 - env['co2']) * o2_sens_abs if env and 'co2' in env else 0
        envf = (o2_suit + ph_suit) / denom

        grid_min = -32768
        grid_max = 32768
        for b in list(self.members):
            b.step(food_nodes, grid_min, grid_max)
            # Reproduction trial
            rep_prob = self.base_rep * capf * envf
            rep_prob = max(0, min(rep_prob, 1))  # Clamp to [0,1]
            # Debug print for reproduction probability
            # print(f"Species {self.name}: rep_prob={rep_prob}, base_rep={self.base_rep}, capf={capf}, envf={envf}")
            if rep_prob > 0 and random.random() < rep_prob and b.food_level > b.food_capacity * 0.5:
                births.append(Bacterium(self, b.x, b.y, size=max(0.5, random.gauss(1.0, 0.2))))
                b.food_level -= b.food_capacity * 0.5
            # Death trial (higher if food is low for big bacteria)
            if b.food_level <= 0:
                deaths.append(b)
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
        self.stats = {k: [] for k in self.stats}
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
        # Drop new food nodes (50x less than before)
        drop_count = max(1, self.drop_frequency // 50)
        for _ in range(drop_count):
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
            pass  # deaths already handled
        self.stats['time'].append(self.step_idx)
        self.stats['total'].append(sum(len(sp.members) for sp in self.species))
        self.stats['births'].append(bsum)
        self.stats['deaths'].append(dsum + (excess if excess > 0 else 0))
        self.species_step_stats.append({'births': births_list, 'deaths': deaths_list})
        for sp in self.species:
            if sp.history:
                sp.history[-1] = len(sp.members)
            else:
                sp.history.append(len(sp.members))
        # Remove old food logic
        total_o2_change = sum(sp.o2_sens * len(sp.members) for sp in self.species) * 0.00005
        self.env['oxygen'] = max(0, min(1, self.env['oxygen'] - total_o2_change))
        self.env['co2'] = max(0, min(1, self.env['co2'] + total_o2_change))
        ph_drift = sum(sp.o2_sens * len(sp.members) for sp in self.species) * 0.0001
        self.env['ph'] = max(0, min(14, self.env['ph'] + ph_drift))
        self.step_idx += 1

    def export(self):
        """
        Write collected stats and per-species data to CSV file for analysis.
        """
        if not self.stats['time']:
            print("No data to export.")
            return
        with open(EXPORT_FILE, 'w', newline='') as f:
            writer = csv.writer(f)
            # Write header for global stats
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
            writer.writerow(['Species', 'Step', 'Count', 'Births', 'Deaths', 'Overall Births', 'Overall Deaths'])
            for idx, sp in enumerate(self.species):
                for step in range(len(sp.history)):
                    births = self.species_step_stats[step]['births'][idx] if step < len(self.species_step_stats) else ''
                    deaths = self.species_step_stats[step]['deaths'][idx] if step < len(self.species_step_stats) else ''
                    writer.writerow([
                        sp.name,
                        step,
                        sp.history[step] if step < len(sp.history) else '',
                        births,
                        deaths,
                        sp.total_births,
                        sp.total_deaths
                    ])

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
        for row in lines[species_start:]:
            if len(row) >= 7:
                name, step, count, births, deaths, overall_births, overall_deaths = row
                # Find species by name
                for sp in self.species:
                    if sp.name == name:
                        step_idx = int(step)
                        while len(sp.history) <= step_idx:
                            sp.history.append(0)
                        sp.history[step_idx] = int(count)
                        sp.total_births = int(overall_births)
                        sp.total_deaths = int(overall_deaths)
                        # Update species_step_stats
                        while len(self.species_step_stats) <= step_idx:
                            self.species_step_stats.append({'births': [0]*len(self.species), 'deaths': [0]*len(self.species)})
                        idx = self.species.index(sp)
                        self.species_step_stats[step_idx]['births'][idx] = int(births)
                        self.species_step_stats[step_idx]['deaths'][idx] = int(deaths)
        print(f"Imported stats from {filename}.")

def load_env_config():
    if os.path.exists(ENV_CONFIG_FILE):
        try:
            with open(ENV_CONFIG_FILE, 'r') as f:
                content = f.read().strip()
                if not content:
                    return None
                return json.loads(content)
        except Exception:
            # If file is corrupt or invalid JSON, ignore and use defaults
            return None
    return None

def save_env_config(env):
    with open(ENV_CONFIG_FILE, 'w') as f:
        json.dump({
            'oxygen': env['oxygen'],
            'co2': env['co2'],
            'ph': env['ph'],
            'capacity': env['capacity'],
            'drop_frequency': env['drop_frequency'],
            'drop_range_min': env['drop_range_min'],
            'drop_range_max': env['drop_range_max']
        }, f)

# ------------------ UI Components ------------------
class Button:
    """
    Clickable button with rectangle bounds, text, and action callback.
    """
    def __init__(self, rect, text, action, enabled=True):
        self.rect = pygame.Rect(rect)
        self.text = text
        self.action = action
        self.enabled = enabled

    def draw(self, s):
        mx, my = pygame.mouse.get_pos()
        # gray out when disabled, otherwise highlight on hover
        color = (50, 50, 50) if not self.enabled else (
            BUTTON_HOVER if self.rect.collidepoint((mx, my)) else BUTTON_COLOR
        )
        pygame.draw.rect(s, color, self.rect, border_radius=5)
        txt = FONT.render(self.text, True, WHITE)
        s.blit(txt, txt.get_rect(center=self.rect.center))

    def click(self, pos):
        if self.enabled and self.rect.collidepoint(pos):
            self.action()

class Slider:
    """
    Horizontal slider controlling an environment variable; draggable knob.
    Supports float (0-1) or integer ranges.
    """
    def __init__(self, rect, label, env_key, minval=0.0, maxval=1.0, is_int=False):
        self.rect = pygame.Rect(rect)
        self.label = label
        self.key = env_key
        self.drag = False
        self.minval = minval
        self.maxval = maxval
        self.is_int = is_int

    def draw(self, s):
        pygame.draw.rect(s, SLIDER_COLOR, self.rect, border_radius=3)
        val = env[self.key]
        norm = (val - self.minval) / (self.maxval - self.minval)
        x = self.rect.x + norm * self.rect.w
        pygame.draw.circle(s, WHITE, (int(x), self.rect.centery), 6)
        disp_val = int(val) if self.is_int else val
        label = FONT.render(f"{self.label}: {disp_val}", True, WHITE)
        s.blit(label, (self.rect.x, self.rect.y - 25))

    def handle_mouse_down(self, pos):
        if self.rect.collidepoint(pos):
            self.drag = True
            return True
        return False

    def handle_mouse_up(self):
        self.drag = False

    def handle_mouse_motion(self, pos):
        if self.drag:
            relx = min(max(pos[0], self.rect.x), self.rect.x + self.rect.w)
            norm = (relx - self.rect.x) / self.rect.w
            val = self.minval + norm * (self.maxval - self.minval)
            if self.is_int:
                val = int(round(val))
            env[self.key] = val

class CellEditor:
    """
    Handles user text entry for editing species table cells.
    Shows only the current input with a blinking caret.
    """
    def __init__(self, col, rect, value, row_idx, col_idx):
        self.col = col
        self.row_idx = row_idx
        self.col_idx = col_idx
        self.rect = pygame.Rect(rect)
        self.text = str(value)
        self.active = True
        self.caret_pos = len(self.text)  # caret at end by default
        self.last_blink = pygame.time.get_ticks()
        self.caret_visible = True

    def handle(self, e):
        if e.type == KEYDOWN:
            if e.key == K_RETURN:
                self.active = False
            elif e.key == K_BACKSPACE:
                if self.caret_pos > 0 and len(self.text) > 0:
                    self.text = self.text[:self.caret_pos-1] + self.text[self.caret_pos:]
                    self.caret_pos -= 1
                # Ensure caret doesn't go negative
                if self.caret_pos < 0:
                    self.caret_pos = 0
            elif e.key == K_LEFT:
                if self.caret_pos > 0:
                    self.caret_pos -= 1
            elif e.key == K_RIGHT:
                if self.caret_pos < len(self.text):
                    self.caret_pos += 1
            elif e.unicode and e.unicode.isprintable():
                if len(self.text) < 10:
                    self.text = self.text[:self.caret_pos] + e.unicode + self.text[self.caret_pos:]
                    self.caret_pos += 1

    def draw(self, s):
        pygame.draw.rect(s, WHITE, self.rect, 2)
        display_text = self.text
        # Blinking caret logic
        now = pygame.time.get_ticks()
        if now - self.last_blink > 500:
            self.caret_visible = not self.caret_visible
            self.last_blink = now
        before = display_text[:self.caret_pos]
        after = display_text[self.caret_pos:]
        txt_before = FONT.render(before, True, WHITE)
        x = self.rect.x + 5 + txt_before.get_width()
        y = self.rect.y + 5
        s.blit(FONT.render(before + after, True, WHITE), (self.rect.x + 5, y))
        if self.active and self.caret_visible:
            caret_height = FONT.get_height()
            pygame.draw.line(s, WHITE, (x, y), (x, y + caret_height), 2)

    def handle_mouse_down(self, pos):
        """Set caret position based on mouse click within the cell."""
        if not self.rect.collidepoint(pos):
            return
        x = pos[0] - self.rect.x - 5
        display_text = self.text
        for i in range(len(display_text) + 1):
            width = FONT.size(display_text[:i])[0]
            if width > x:
                self.caret_pos = max(0, i - 1)
                break
        else:
            self.caret_pos = len(display_text)

class App:
    """
    Main application class: sets up pygame, UI, simulation,
    and dispatches events to draw/update accordingly.
    """
    def __init__(self):
        pygame.init()
        global FONT
        FONT = pygame.font.SysFont('Consolas', 20)
        # create a resizable fullscreen window
        info = pygame.display.Info()
        flags = RESIZABLE | FULLSCREEN
        self.w, self.h = info.current_w, info.current_h
        self.screen = pygame.display.set_mode((self.w, self.h), flags)
        self.clock = pygame.time.Clock()

        # Data manager for species parameters
        self.dm = SpeciesDataManager()
        # Field descriptions for species tab
        self.descriptions = {
            'name': 'Species name (string)',
            'init': 'Initial count (int ≥ 0)',
            'rep': 'Reproduction rate (0–1)',
            'death': 'Death rate (0–1)',
            'so': 'O₂ sensitivity (-1 to 1; <0 needs CO₂, >0 needs O₂)',
            'sf': 'Food sensitivity (0–1)',
            'speed': 'Speed (pixels/sec)',
            'ph': 'Preferred pH (0–14)'
        }
        # Add new fields to species table
        self.descriptions.update({
            'food_capacity': 'Max food storage (units)',
            'metabolic_rate': 'Food used per tick',
            'consumption_speed': 'Food eaten per tick',
            'greediness': '0.0–1.0 preference for rich food even if far'
        })
        if 'food_capacity' not in self.dm.fieldnames:
            self.dm.fieldnames += ['food_capacity', 'metabolic_rate', 'consumption_speed', 'greediness']

        # Undo and redo stacks for editing species table
        self.undo_stack = []
        self.redo_stack = []

        # Abbreviated fieldnames for species tab
        self.ABBREV_FIELDNAMES = {
            'name': 'Name',
            'init': 'Init',
            'rep': 'Rep',
            'death': 'Death',
            'so': 'O2',
            'sf': 'Food',
            'speed': 'Spd',
            'ph': 'pH',
            'food_capacity': 'Cap',
            'metabolic_rate': 'Meta',
            'consumption_speed': 'Eat',
            'greediness': 'Greed'
        }

        # Add scroll state for species tab
        self.species_scroll_x = 0
        self.species_scroll_y = 0
        self.species_max_visible_cols = max(1, (self.w - 100) // 120)
        self.species_max_visible_rows = max(1, (self.h - 180) // 30)
        self.species_scrollbar_thickness = 16

        # Add scroll state for stats tab
        self.stats_scroll_y = 0
        self.stats_max_visible_rows = max(1, (self.h - 180) // 28)
        self.stats_scrollbar_thickness = 16

        # Show/hide species attribute descriptions dialog
        self.show_species_desc = False
        # Button in bottom right corner
        self.species_desc_btn = Button((self.w - 200, self.h - 60, 180, 30), 'Show Desc', self.toggle_species_desc)
        # Dialog position will be set dynamically in draw_species_tab
        self.species_desc_dialog_rect = None

        self.init_sim()
        # Create buttons, tabs, sliders
        # Use consistent horizontal spacing for buttons
        button_w, button_h, button_gap = 80, 30, 10
        button_y = 20
        button_xs = [20 + i * (button_w + button_gap) for i in range(8)]
        self.play_btn = Button((button_xs[0], button_y, button_w, button_h), 'Play', self.play)
        self.import_btn = Button((button_xs[1], button_y, button_w, button_h), 'Import', self.import_data)
        self.pause_btn = Button((button_xs[2], button_y, button_w, button_h), 'Pause', self.toggle_pause)
        self.export_btn = Button((button_xs[3], button_y, button_w, button_h), 'Export', self.sim.export, enabled=False)
        self.add_btn = Button((button_xs[4], button_y, button_w, button_h), 'Add', self.add_species)
        self.remove_btn = Button((button_xs[5], button_y, button_w, button_h), 'Remove', self.remove_species)
        self.undo_btn = Button((button_xs[6], button_y, button_w, button_h), 'Undo', self.undo, enabled=False)
        self.redo_btn = Button((button_xs[7], button_y, button_w, button_h), 'Redo', self.redo, enabled=False)

        self.buttons = [
            self.play_btn,
            self.import_btn,
            self.pause_btn,
            self.export_btn,
            self.add_btn,
            self.remove_btn,
            self.undo_btn,
            self.redo_btn
        ]

        self.tabs = {
            'grid': Button((650, 20, 80, 30), 'Grid', lambda: self.set_tab('grid')),
            'graph': Button((740, 20, 80, 30), 'Graph', lambda: self.set_tab('graph')),
            'stats': Button((830, 20, 80, 30), 'Stats', lambda: self.set_tab('stats')),
            'species': Button((920, 20, 80, 30), 'Species', lambda: self.set_tab('species'))
        }
        self.tab = 'grid'
        # Align sliders to top right, flush with button text
        slider_x = self.w - 220  # 200 width + 20 margin
        slider_y_start = 60      # Increased from 20 to 60 for visibility
        slider_gap = 40
        # Add new sliders for drop frequency and drop range
        self.sliders = [
            Slider((slider_x, slider_y_start + 0*slider_gap, 200, 10), 'O2', 'oxygen'),
            Slider((slider_x, slider_y_start + 1*slider_gap, 200, 10), 'CO2', 'co2'),
            Slider((slider_x, slider_y_start + 2*slider_gap, 200, 10), 'Drop Freq', 'drop_frequency', minval=1, maxval=20, is_int=True),
            Slider((slider_x, slider_y_start + 3*slider_gap, 200, 10), 'Drop Range Min', 'drop_range_min', minval=1, maxval=20, is_int=True),
            Slider((slider_x, slider_y_start + 4*slider_gap, 200, 10), 'Drop Range Max', 'drop_range_max', minval=1, maxval=50, is_int=True),
            Slider((slider_x, slider_y_start + 5*slider_gap, 200, 10), 'pH', 'ph', minval=0, maxval=14)
        ]
        self.editors = []  # active cell editors
        self.selected_row = 0
        self.selected_col = 0
        self.dragged_slider = None

        # --- Grid/scrolling state ---
        self.GRID_SIZE = 65536
        self.GRID_HALF = self.GRID_SIZE // 2
        self.grid_center_x = 0
        self.grid_center_y = 0
        self.scroll_speed = 40  # pixels per keypress/mousewheel
        self.scrollbar_thickness = 16

        # For mouse drag scrolling
        self.dragging = False
        self.drag_start = None
        self.grid_start = None

        # Minimap state
        self.minimap_size = 160  # pixels (square)
        self.minimap_margin = 20
        self.minimap_dragging = False
        self.grid_dragging = False
        self.grid_drag_start = None
        self.grid_center_start = None

    def clamp_grid_center(self):
        # Clamp so the view never goes outside the grid
        half_w = self.w // 2
        half_h = self.h // 2
        min_x = -self.GRID_HALF + half_w
        max_x = self.GRID_HALF - half_w
        min_y = -self.GRID_HALF + half_h
        max_y = self.GRID_HALF - half_h
        self.grid_center_x = max(min(self.grid_center_x, max_x), min_x)
        self.grid_center_y = max(min(self.grid_center_y, max_y), min_y)

    def handle_grid_scroll(self, dx, dy):
        self.grid_center_x += dx
        self.grid_center_y += dy
        self.clamp_grid_center()

    def draw_grid(self):
        # Calculate top-left of view in grid coordinates
        half_w = self.w // 2
        half_h = self.h // 2
        view_left = self.grid_center_x - half_w
        view_top = self.grid_center_y - half_h
        view_right = self.grid_center_x + half_w
        view_bottom = self.grid_center_y + half_h

        # Draw white border if visible
        border_rect = pygame.Rect(
            self.w//2 - self.grid_center_x - self.GRID_HALF,
            self.h//2 - self.grid_center_y - self.GRID_HALF,
            self.GRID_SIZE, self.GRID_SIZE
        )
        pygame.draw.rect(self.screen, WHITE, border_rect, 2)

        # Optionally draw grid lines or bacteria here

        # Draw scrollbars
        self.draw_scrollbars()

        # Draw current center coordinates at bottom right
        coord_text = f"Center: ({self.grid_center_x}, {self.grid_center_y})"
        txt = FONT.render(coord_text, True, WHITE)
        self.screen.blit(txt, (self.w - txt.get_width() - 20, self.h - txt.get_height() - 20))

    def draw_scrollbars(self):
        # Horizontal scrollbar
        bar_w = self.w - self.scrollbar_thickness
        bar_x = self.scrollbar_thickness // 2
        bar_y = self.h - self.scrollbar_thickness
        pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, bar_w, self.scrollbar_thickness), border_radius=6)
        # Thumb position
        thumb_w = max(bar_w * self.w / self.GRID_SIZE, 40)
        thumb_x = bar_x + (self.grid_center_x + self.GRID_HALF - self.w//2) / (self.GRID_SIZE - self.w) * (bar_w - thumb_w)
        pygame.draw.rect(self.screen, WHITE, (thumb_x, bar_y+2, thumb_w, self.scrollbar_thickness-4), border_radius=6)

        # Vertical scrollbar
        bar_h = self.h - self.scrollbar_thickness
        bar_x = self.w - self.scrollbar_thickness
        bar_y = self.scrollbar_thickness // 2
        pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, self.scrollbar_thickness, bar_h), border_radius=6)
        thumb_h = max(bar_h * self.h / self.GRID_SIZE, 40)
        thumb_y = bar_y + (self.grid_center_y + self.GRID_HALF - self.h//2) / (self.GRID_SIZE - self.h) * (bar_h - thumb_h)
        pygame.draw.rect(self.screen, WHITE, (bar_x+2, thumb_y, self.stats_scrollbar_thickness-4, thumb_h), border_radius=6)

    def draw_minimap(self):
        # Position minimap above horizontal scrollbar, bottom left
        x = self.minimap_margin
        y = self.h - self.scrollbar_thickness - self.minimap_size - self.minimap_margin
        size = self.minimap_size
        grid_min = -self.GRID_HALF
        grid_max = self.GRID_HALF
        grid_size = self.GRID_SIZE

        # Draw minimap background
        pygame.draw.rect(self.screen, (40, 40, 40), (x, y, size, size), border_radius=8)
        pygame.draw.rect(self.screen, WHITE, (x, y, size, size), 2)

        # Draw bacteria as dots
        for sp in self.sim.species:
            for b in sp.members:
                mx = x + int((b.x - grid_min) / grid_size * size)
                my = y + int((b.y - grid_min) / grid_size * size)
                if 0 <= mx < x + size and 0 <= my < y + size:
                    pygame.draw.circle(self.screen, sp.color, (mx, my), 2)

        # Draw current view rectangle
        half_w = self.w // 2
        half_h = self.h // 2
        left = self.grid_center_x - half_w
        top = self.grid_center_y - half_h
        right = self.grid_center_x + half_w
        bottom = self.grid_center_y + half_h
        rect_x = x + int((left - grid_min) / grid_size * size)
        rect_y = y + int((top - grid_min) / grid_size * size)
        rect_w = max(2, int((right - left) / grid_size * size))
        rect_h = max(2, int((bottom - top) / grid_size * size))
        pygame.draw.rect(self.screen, WHITE, (rect_x, rect_y, rect_w, rect_h), 2)

        # Store minimap rect for event handling
        self._minimap_rect = pygame.Rect(x, y, size, size)

    def minimap_to_grid(self, mx, my):
        # Convert minimap pixel to grid center coordinates
        x = self._minimap_rect.x
        y = self._minimap_rect.y
        size = self._minimap_rect.w
        grid_min = -self.GRID_HALF
        grid_max = self.GRID_HALF
        gx = grid_min + (mx - x) / size * self.GRID_SIZE
        gy = grid_min + (my - y) / size * self.GRID_SIZE
        return int(gx), int(gy)

    def handle_grid_drag(self, mx, my):
        if self.grid_drag_start is not None and self.grid_center_start is not None:
            dx = mx - self.grid_drag_start[0]
            dy = my - self.grid_drag_start[1]
            self.grid_center_x = self.grid_center_start[0] - dx
            self.grid_center_y = self.grid_center_start[1] - dy
            self.clamp_grid_center()

    def init_sim(self):
        cols = [(200, 50, 50), (50, 200, 50), (50, 50, 200), (200, 200, 50)]
        spp = []
        for i, d in enumerate(self.dm.species):
            try:
                spp.append(BacteriaSpecies(
                    d['name'],
                    int(d['init']),
                    float(d['rep']),
                    float(d['death']),
                    float(d['so']),
                    float(d['sf']),
                    float(d.get('speed', 40)),
                    float(d.get('ph', 7)),
                    cols[i % len(cols)],
                    float(d.get('food_capacity', 10)),
                    float(d.get('metabolic_rate', 0.1)),
                    float(d.get('consumption_speed', 1.0)),
                    float(d.get('greediness', 0.5))
                ))
            except Exception:
                continue
        if not spp:
            spp.append(BacteriaSpecies('Ecoli', 50, 0.1, 0.05, 1, 1, 40, 7, cols[0], 10, 0.1, 1.0, 0.5))
        global env
        loaded_env = load_env_config()
        if loaded_env:
            env = loaded_env
        else:
            env = {
                'oxygen': 0.8,
                'co2': 0.2,
                'ph': 7.0,
                'capacity': 500,
                'drop_frequency': 5,
                'drop_range_min': 2,
                'drop_range_max': 8
            }
        env['drop_range'] = (env['drop_range_min'], env['drop_range_max'])
        self.sim = Simulation(spp, env, cap=env['capacity'], steps=500)

    def add_species(self):
        """Append new blank species and re-init simulation, with undo/redo support."""
        blank = {k: '0' for k in self.dm.fieldnames}
        blank['name'] = ''
        idx = len(self.dm.species)
        self.dm.species.append(blank)
        self.dm.save()
        self.init_sim()
        # Undo support
        self.undo_stack.append({'action': 'add_row', 'index': idx, 'row': blank})
        self.redo_stack.clear()
        self.undo_btn.enabled = True
        self.redo_btn.enabled = bool(self.redo_stack)

    def remove_species(self):
        """Remove currently selected species row and re-init, with undo/redo support."""
        idx = self.selected_row if self.selected_row is not None else -1
        if 0 <= idx < len(self.dm.species):
            removed_row = self.dm.species[idx].copy()
            self.dm.remove(idx)
            self.init_sim()
            if self.dm.species:
                self.selected_row = min(self.selected_row, len(self.dm.species) - 1)
            else:
                self.selected_row = 0
            self.undo_stack.append({'action': 'remove_row', 'index': idx, 'row': removed_row})
            self.redo_stack.clear()
            self.undo_btn.enabled = True
            self.redo_btn.enabled = bool(self.redo_stack)

    def undo(self):
        """Undo last species table edit or row add/remove."""
        if not self.undo_stack:
            return
        action = self.undo_stack.pop()
        if action.get('action') == 'add_row':
            idx = action['index']
            if 0 <= idx < len(self.dm.species):
                removed_row = self.dm.species[idx].copy()
                self.dm.remove(idx)
                self.init_sim()
                if self.dm.species:
                    self.selected_row = min(self.selected_row, len(self.dm.species) - 1)
                else:
                    self.selected_row = 0
                self.redo_stack.append({'action': 'add_row', 'index': idx, 'row': action['row']})
        elif action.get('action') == 'remove_row':
            idx = action['index']
            if 0 <= idx <= len(self.dm.species):
                self.dm.species.insert(idx, action['row'])
                self.dm.save()
                self.init_sim()
                self.selected_row = idx
                self.redo_stack.append({'action': 'remove_row', 'index': idx, 'row': action['row']})
        else:
            row, col, old = action['row'], action['col'], action['old']
            if 0 <= row < len(self.dm.species) and col in self.dm.fieldnames:
                new = self.dm.species[row][col]
                self.dm.species[row][col] = old
                self.dm.save()
                self.init_sim()
                self.redo_stack.append({'action': 'cell_edit', 'row': row, 'col': col, 'old': old, 'new': new})
        self.undo_btn.enabled = bool(self.undo_stack)
        self.redo_btn.enabled = bool(self.redo_stack)

    def redo(self):
        """Redo last undone add/remove/edit."""
        if not self.redo_stack:
            return
        action = self.redo_stack.pop()
        if action.get('action') == 'add_row':
            idx = action['index']
            self.dm.species.insert(idx, action['row'])
            self.dm.save()
            self.init_sim()
            self.selected_row = idx
            self.undo_stack.append({'action': 'add_row', 'index': idx, 'row': action['row']})
        elif action.get('action') == 'remove_row':
            idx = action['index']
            if 0 <= idx < len(self.dm.species):
                removed_row = self.dm.species[idx].copy()
                self.dm.remove(idx)
                self.init_sim()
                if self.dm.species:
                    self.selected_row = min(self.selected_row, len(self.dm.species) - 1)
                else:
                    self.selected_row = 0
                self.undo_stack.append({'action': 'remove_row', 'index': idx, 'row': removed_row})
        elif action.get('action') == 'cell_edit':
            row, col, new = action['row'], action['col'], action['new']
            if 0 <= row < len(self.dm.species) and col in self.dm.fieldnames:
                old = self.dm.species[row][col]
                self.dm.species[row][col] = new
                self.dm.save()
                self.init_sim()
                self.undo_stack.append({'action': 'cell_edit', 'row': row, 'col': col, 'old': old, 'new': new})
        self.undo_btn.enabled = bool(self.undo_stack)
        self.redo_btn.enabled = bool(self.redo_stack)

    def play(self):
        """Start simulation run if species exist and force unpause."""
        if not self.dm.species:
            print("Add at least one species before playing.")
            return
        self.init_sim()
        self.sim.running = True
        self.pause_btn.text = 'Pause'  # Ensure Pause button shows correct state

    def toggle_pause(self):
        """Toggle simulation pause/unpause and update button label."""
        self.sim.running = not self.sim.running
        if self.sim.running:
            self.pause_btn.text = 'Pause'
        else:
            self.pause_btn.text = 'Unpause'

    def set_tab(self, t):
        """Switch UI tab (grid, graph, stats, species)."""
        self.tab = t
        self.editors = []
        # Auto-hide species desc dialog when leaving species tab
        if t != 'species':
            self.show_species_desc = False

    def draw_species_tab(self):
        start_x, start_y = 50, 100
        cell_w, cell_h = 120, 30
        cols = self.dm.fieldnames
        abbrev = self.ABBREV_FIELDNAMES
        # Calculate visible columns/rows
        total_cols = len(cols)
        total_rows = len(self.dm.species)
        max_cols = self.species_max_visible_cols
        max_rows = self.species_max_visible_rows
        col_offset = self.species_scroll_x
        row_offset = self.species_scroll_y
        # Draw headers
        for j in range(col_offset, min(col_offset + max_cols, total_cols)):
            col = cols[j]
            head = FONT.render(abbrev.get(col, col), True, WHITE)
            self.screen.blit(head, (start_x + (j-col_offset)*cell_w, start_y - cell_h))
        # Draw cells
        for i in range(row_offset, min(row_offset + max_rows, total_rows)):
            row = self.dm.species[i]
            for j in range(col_offset, min(col_offset + max_cols, total_cols)):
                col = cols[j]
                rect = (start_x + (j-col_offset)*cell_w, start_y + (i-row_offset)*cell_h, cell_w, cell_h)
                highlight = (i == self.selected_row and j == self.selected_col)
                pygame.draw.rect(
                    self.screen,
                    SLIDER_COLOR if highlight else WHITE,
                    rect,
                    2 if not highlight else 0
                )
                if highlight and not self.editors:
                    overlay = pygame.Surface((cell_w,cell_h), pygame.SRCALPHA)
                    overlay.fill((50,255,50,50))
                    self.screen.blit(overlay, (rect[0], rect[1]))
                if not (self.editors and self.editors[0].row_idx == i and self.editors[0].col_idx == j):
                    txt = FONT.render(str(row[col]), True, WHITE)
                    self.screen.blit(txt, (rect[0] + 5, rect[1] + 5))
        # Draw horizontal scrollbar if needed
        if total_cols > max_cols:
            bar_w = max_cols * cell_w
            bar_x = start_x
            bar_y = start_y + max_rows * cell_h
            pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, bar_w, self.species_scrollbar_thickness), border_radius=6)
            thumb_w = max(bar_w * max_cols / total_cols, 40)
            thumb_x = bar_x + (col_offset / (total_cols - max_cols)) * (bar_w - thumb_w) if total_cols > max_cols else bar_x
            pygame.draw.rect(self.screen, WHITE, (thumb_x, bar_y+2, thumb_w, self.species_scrollbar_thickness-4), border_radius=6)
        # Draw vertical scrollbar if needed
        if total_rows > max_rows:
            bar_h = max_rows * cell_h
            bar_x = start_x + max_cols * cell_w
            bar_y = start_y
            pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, self.species_scrollbar_thickness, bar_h), border_radius=6)
            thumb_h = max(bar_h * max_rows / total_rows, 40)
            thumb_y = bar_y + (row_offset / (total_rows - max_rows)) * (bar_h - thumb_h) if total_rows > max_rows else bar_y
            pygame.draw.rect(self.screen, WHITE, (bar_x+2, thumb_y, self.stats_scrollbar_thickness-4, thumb_h), border_radius=6)
        # Draw show/hide descriptions button in bottom right
        self.species_desc_btn.rect.x = self.w - 200
        self.species_desc_btn.rect.y = self.h - 60
        # Ensure button text always matches dialog state
        self.species_desc_btn.text = 'Hide Desc' if self.show_species_desc else 'Show Desc'
        self.species_desc_btn.draw(self.screen)
        # Draw descriptions dialog if shown, above and to the right of button
        if self.show_species_desc:
            dialog_w, dialog_h = 400, 400
            dialog_x = self.species_desc_btn.rect.x + self.species_desc_btn.rect.w + 10
            dialog_y = self.species_desc_btn.rect.y - dialog_h - 10
            # Clamp dialog to window bounds
            if dialog_x + dialog_w > self.w:
                dialog_x = self.w - dialog_w - 10
            if dialog_y < 0:
                dialog_y = 10
            self.species_desc_dialog_rect = pygame.Rect(dialog_x, dialog_y, dialog_w, dialog_h)
            pygame.draw.rect(self.screen, (40,40,40), self.species_desc_dialog_rect, border_radius=8)
            pygame.draw.rect(self.screen, WHITE, self.species_desc_dialog_rect, 2)
            x, y = self.species_desc_dialog_rect.x + 20, self.species_desc_dialog_rect.y + 20
            for k, field in enumerate(cols):
                desc = self.descriptions.get(field, '')
                label = f"{abbrev.get(field, field)}: {desc}"
                text = FONT.render(label, True, WHITE)
                self.screen.blit(text, (x, y + k*32))

    def handle_species_keys(self, e):
        """Keyboard navigation within species table, with scroll offsets."""
        if self.editors:
            return
        rows = len(self.dm.species)
        cols = len(self.dm.fieldnames)
        if rows == 0:
            return
        max_rows = self.species_max_visible_rows
        max_cols = self.species_max_visible_cols
        # Navigation
        if e.key == K_DOWN and self.selected_row < rows-1:
            self.selected_row += 1
            if self.selected_row >= self.species_scroll_y + max_rows:
                self.species_scroll_y = min(self.species_scroll_y + 1, rows - max_rows)
        if e.key == K_UP and self.selected_row > 0:
            self.selected_row -= 1
            if self.selected_row < self.species_scroll_y:
                self.species_scroll_y = max(self.species_scroll_y - 1, 0)
        if e.key == K_RIGHT and self.selected_col < cols-1:
            self.selected_col += 1
            if self.selected_col >= self.species_scroll_x + max_cols:
                self.species_scroll_x = min(self.species_scroll_x + 1, cols - max_cols)
        if e.key == K_LEFT and self.selected_col > 0:
            self.selected_col -= 1
            if self.selected_col < self.species_scroll_x:
                self.species_scroll_x = max(self.species_scroll_x - 1, 0)
        if e.key == K_RETURN:
            i, j = self.selected_row, self.selected_col
            if 0 <= i < rows and 0 <= j < cols:
                val = self.dm.species[i][self.dm.fieldnames[j]]
                rect = (50 + (j-self.species_scroll_x)*120, 100 + (i-self.species_scroll_y)*30, 120, 30)
                self.editors = [CellEditor(self.dm.fieldnames[j], rect, val, i, j)]

    def is_mouse_over_ui(self, mx, my):
        # Check buttons
        for b in self.buttons:
            if b.rect.collidepoint((mx, my)):
                return True
        # Check tabs
        for b in self.tabs.values():
            if b.rect.collidepoint((mx, my)):
                return True
        # Check sliders
        for s in self.sliders:
            if s.rect.collidepoint((mx, my)):
                return True
        # Check minimap
        if hasattr(self, '_minimap_rect') and self._minimap_rect.collidepoint(mx, my):
            return True
        # Add other dialogs here as needed
        return False

    @staticmethod
    def normalize_numeric_string(s):
        s = s.strip()
        # Add leading zero if string starts with '.' or '-.'
        if s.startswith('.'):
            s = '0' + s
        elif s.startswith('-.'):
            s = '-0.' + s[2:]
        # Remove leading zeros for numbers (but keep single zero if needed)
        # Use regex to match valid float/int
        match = re.match(r'^(-?)(0*)(\d*)(\.?\d*)$', s)
        if match:
            sign, zeros, intpart, decpart = match.groups()
            # If only zeros before decimal, keep one zero
            if intpart == '' and zeros:
                intpart = '0'
            num_str = f"{sign}{intpart}{decpart}"
            try:
                if '.' in num_str:
                    num_str = str(float(num_str))
                else:
                    num_str = str(int(num_str))
                return num_str
            except Exception:
                pass
        # Fallback: try float/int conversion
        try:
            if '.' in s:
                return str(float(s))
            else:
                return str(int(s))
        except Exception:
            return s

    def run(self):
        """Main loop: handle events, update sim, and draw UI."""
        while True:
            for e in pygame.event.get():
                if e.type == QUIT:
                    pygame.quit(); sys.exit()
                if e.type == VIDEORESIZE:
                    self.w, self.h = e.w, e.h
                    self.screen = pygame.display.set_mode((self.w, self.h), RESIZABLE)
                    self.clamp_grid_center()
                if e.type == MOUSEBUTTONDOWN:
                    mx, my = e.pos
                    # --- Always process button and tab clicks ---
                    for b in self.buttons + list(self.tabs.values()):
                        b.click(e.pos)
                    # --- Always process slider drag ---
                    for s in self.sliders:
                        if s.handle_mouse_down(e.pos):
                            self.dragged_slider = s
                    # --- Minimap click/drag ---
                    if self.tab == 'grid' and hasattr(self, '_minimap_rect') and self._minimap_rect.collidepoint(mx, my):
                        self.minimap_dragging = True
                        gx, gy = self.minimap_to_grid(mx, my)
                        self.grid_center_x, self.grid_center_y = gx, gy
                        self.clamp_grid_center()
                    # --- Main grid clickdrag (only if not over UI) ---
                    elif self.tab == 'grid' and not self.is_mouse_over_ui(mx, my):
                        self.grid_dragging = True
                        self.grid_drag_start = (mx, my)
                        self.grid_center_start = (self.grid_center_x, self.grid_center_y)
                    # --- Scrollbar drag ---
                    if self.tab == 'grid':
                        if my >= self.h - self.scrollbar_thickness:
                            self.dragging = 'h'
                            self.drag_start = mx
                            self.grid_start = self.grid_center_x
                        elif mx >= self.w - self.scrollbar_thickness:
                            self.dragging = 'v'
                            self.drag_start = my
                            self.grid_start = self.grid_center_y
                    # --- Species table cell click ---
                    if self.tab == 'species':
                        if self.species_desc_btn.rect.collidepoint(e.pos):
                            self.species_desc_btn.click(e.pos)
                            continue
                        if self.editors:
                            self.editors[0].handle_mouse_down(e.pos)
                        else:
                            self.handle_mouse_species(e.pos)
                if e.type == MOUSEBUTTONUP:
                    if self.dragged_slider:
                        self.dragged_slider.handle_mouse_up()
                        self.dragged_slider = None
                    self.dragging = False
                    self.minimap_dragging = False
                    self.grid_dragging = False
                    self.grid_drag_start = None
                    self.grid_center_start = None
                if e.type == MOUSEMOTION:
                    mx, my = e.pos
                    if self.dragged_slider:
                        self.dragged_slider.handle_mouse_motion(e.pos)
                    # Minimap drag
                    if self.minimap_dragging and self.tab == 'grid':
                        gx, gy = self.minimap_to_grid(mx, my)
                        self.grid_center_x, self.grid_center_y = gx, gy
                        self.clamp_grid_center()
                    # Main grid drag
                    if self.grid_dragging and self.tab == 'grid':
                        self.handle_grid_drag(mx, my)
                    # Scrollbar drag
                    if self.dragging and self.tab == 'grid':
                        if self.dragging == 'h':
                            bar_w = self.w - self.scrollbar_thickness
                            delta = mx - self.drag_start
                            grid_range = self.GRID_SIZE - self.w
                            self.grid_center_x = self.grid_start + int(delta * grid_range / (bar_w - max(bar_w * self.w / self.GRID_SIZE, 40)))
                            self.clamp_grid_center()
                        elif self.dragging == 'v':
                            bar_h = self.h - self.scrollbar_thickness
                            delta = my - self.drag_start
                            grid_range = self.GRID_SIZE - self.h
                            self.grid_center_y = self.grid_start + int(delta * grid_range / (bar_h - max(bar_h * self.h / self.GRID_SIZE, 40)))
                            self.clamp_grid_center()
                if e.type == KEYDOWN:
                    # Ctrl+Z undo
                    if e.key == K_z and pygame.key.get_mods() & KMOD_CTRL and not (pygame.key.get_mods() & KMOD_SHIFT):
                        if self.tab == 'species' and self.undo_btn.enabled and not self.editors:
                            self.undo()
                        continue
                    # Ctrl+Shift+Z redo
                    if e.key == K_z and pygame.key.get_mods() & KMOD_CTRL and pygame.key.get_mods() & KMOD_SHIFT:
                        if self.tab == 'species' and self.redo_btn.enabled and not self.editors:
                            self.redo()
                        continue
                    # Remove species with Delete
                    if e.key == K_DELETE:
                        if self.tab == 'species' and self.remove_btn.enabled and not self.editors:
                            self.remove_species()
                        continue
                    # Add species with "+" key (main and numpad)
                    if e.key in (K_PLUS, K_KP_PLUS):
                        if self.tab == 'species' and self.add_btn.enabled and not self.editors:
                            self.add_species()
                        continue

                    # --- NEW LOGIC: Press Enter to start editing ---
                    if self.tab == 'species' and not self.editors and e.key == K_RETURN:
                        i, j = self.selected_row, self.selected_col
                        rect = (50 + j*120, 100 + i*30, 120, 30)
                        val = self.dm.species[i][self.dm.fieldnames[j]]
                        self.editors = [CellEditor(self.dm.fieldnames[j], rect, val, i, j)]
                        continue

                    # Navigation keys (arrow keys, etc.)
                    if self.tab == 'grid':
                        if e.key == K_LEFT or e.key == K_a:
                            self.handle_grid_scroll(-self.scroll_speed, 0)
                        if e.key == K_RIGHT or e.key == K_d:
                            self.handle_grid_scroll(self.scroll_speed, 0)
                        if e.key == K_UP or e.key == K_w:
                            self.handle_grid_scroll(0, -self.scroll_speed)
                        if e.key == K_DOWN or e.key == K_s:
                            self.handle_grid_scroll(0, self.scroll_speed)
                    if self.tab == 'species' and not self.editors:
                        self.handle_species_keys(e)

                    # Editing logic
                    for ed in self.editors:
                        ed.handle(e)
                    if e.key == K_ESCAPE and self.editors:
                        self.editors = []
                    if e.key == K_RETURN and self.editors and not self.editors[0].active:
                        ed = self.editors[0]
                        val = ed.text
                        i, j = ed.row_idx, ed.col_idx
                        col = ed.col
                        old_val = self.dm.species[i][col]
                        # --- Normalize numeric input for non-name columns ---
                        new_val = val
                        if col != 'name':
                            new_val = App.normalize_numeric_string(val)
                        if col == 'name':
                            self.undo_stack.append({'action': 'cell_edit', 'row': i, 'col': col, 'old': old_val, 'new': new_val})
                            self.dm.species[i][col] = new_val
                        else:
                            try:
                                # Validation for new columns
                                if col == 'init':
                                    num = int(new_val)
                                    if num < 0:
                                        raise ValueError
                                elif col in ['rep', 'death', 'sf']:
                                    num = float(new_val)
                                    if not 0 <= num <= 1:
                                        raise ValueError
                                elif col == 'so':
                                    num = float(new_val)
                                    if not -1 <= num <= 1:
                                        raise ValueError
                                elif col == 'speed':
                                    num = float(new_val)
                                    if num < 0:
                                        raise ValueError
                                elif col == 'ph':
                                    num = float(new_val)
                                    if not 0 <= num <= 14:
                                        raise ValueError
                                elif col == 'greediness':
                                    num = float(new_val)
                                    if not 0 <= num <= 1:
                                        raise ValueError
                                # apply change if valid
                                self.undo_stack.append({'action': 'cell_edit', 'row': i, 'col': col, 'old': old_val, 'new': new_val})
                                self.dm.species[i][col] = new_val
                            except Exception:
                                pass
                        self.dm.save()
                        self.editors = []
                        self.redo_stack.clear()
                        self.undo_btn.enabled = True

                if e.type == MOUSEWHEEL:
                    if self.tab == 'grid':
                        mods = pygame.key.get_mods()
                        # Shift for horizontal scroll
                        if mods & KMOD_SHIFT:
                            self.handle_grid_scroll(e.y * self.scroll_speed, 0)
                        else:
                            self.handle_grid_scroll(0, -e.y * self.scroll_speed)
                        # Ctrl for left/right scroll (optional)
                        if mods & KMOD_CTRL:
                            self.handle_grid_scroll(e.x * self.scroll_speed, 0)
                # --- Species tab scroll handling ---
                if self.tab == 'species':
                    # Mouse wheel vertical scroll
                    if e.type == MOUSEWHEEL:
                        if e.y != 0:
                            total_rows = len(self.dm.species)
                            max_rows = self.species_max_visible_rows
                            self.species_scroll_y = max(0, min(self.species_scroll_y - e.y, max(0, total_rows - max_rows)))
                        if e.x != 0:
                            total_cols = len(self.dm.fieldnames)
                            max_cols = self.species_max_visible_cols
                            self.species_scroll_x = max(0, min(self.species_scroll_x - e.x, max(0, total_cols - max_cols)))
                    # Mouse drag scrollbars
                    # TODO: Add mouse drag logic for scrollbars
                # --- Stats tab scroll handling ---
                if self.tab == 'stats':
                    if e.type == MOUSEWHEEL:
                        total_rows = len(self.sim.species) * (5 + 1) + 6
                        max_rows = self.stats_max_visible_rows
                        self.stats_scroll_y = max(0, min(self.stats_scroll_y - e.y, max(0, total_rows - max_rows)))

            # --- Enable/disable export button based on data ---
            self.export_btn.enabled = bool(self.sim.stats['time'])
            # --- Update drop range in env from sliders ---
            # Clamp drop_range_min and drop_range_max to maintain max >= min
            if env['drop_range_min'] > env['drop_range_max']:
                env['drop_range_max'] = env['drop_range_min']
            if env['drop_range_max'] < env['drop_range_min']:
                env['drop_range_min'] = env['drop_range_max']
            env['drop_range'] = (env['drop_range_min'], env['drop_range_max'])
            self.sim.drop_frequency = int(env['drop_frequency'])
            self.sim.drop_range = env['drop_range']
            # --- Always update simulation every frame ---
            save_env_config(env)  # Save env after any change
            self.sim.update(self.w, self.h)

            # Draw background
            self.screen.fill(BACKGROUND)
            # Draw content based on current tab
            if self.tab == 'grid':
                self.draw_grid()
                self.draw_minimap()
                # Draw bacteria in grid coordinates
                for sp in self.sim.species:
                    for b in sp.members:
                        # Convert grid coordinates to screen
                        sx = self.w//2 + (b.x - self.grid_center_x)
                        sy = self.h//2 + (b.y - self.grid_center_y)
                        # Only draw if visible
                        if 0 <= sx < self.w and 0 <= sy < self.h:
                            pygame.draw.circle(self.screen, sp.color, (int(sx), int(sy)), 4)
                # Draw food nodes in grid
                for node in self.sim.food_nodes:
                    sx = self.w//2 + (node.x - self.grid_center_x)
                    sy = self.h//2 + (node.y - self.grid_center_y)
                    if 0 <= sx < self.w and 0 <= sy < self.h:
                        pygame.draw.circle(self.screen, (200,200,70), (int(sx), int(sy)), 6)
            elif self.tab == 'graph' and self.sim.stats['time']:
                mx, my = int(self.w * 0.05), int(self.h * 0.2)
                w, h = int(self.w * 0.6), int(self.h * 0.6)
                pygame.draw.rect(self.screen, WHITE, (mx, my, w, h), 2)

                times = self.sim.stats['time']
                total_vals = self.sim.stats['total']
                # Patch: avoid division by zero for max_t and max_v
                max_t = max(times) if times else 1
                if max_t == 0:
                    max_t = 1
                max_v_candidates = [max(total_vals)] if total_vals else [1]
                for sp in self.sim.species:
                    if sp.history:
                        max_v_candidates.append(max(sp.history))
                max_v = max(max_v_candidates)
                if max_v == 0:
                    max_v = 1

                font = FONT
                num_xticks = min(10, max_t + 1)
                for i in range(num_xticks + 1):
                    step = int(i * max_t / num_xticks)
                    # Patch: avoid division by zero for max_t
                    x = mx + (step / max_t) * w if max_t else mx
                    pygame.draw.line(self.screen, WHITE, (x, my + h), (x, my + h + 8), 2)
                    label = font.render(str(step), True, WHITE)
                    self.screen.blit(label, (x - label.get_width() // 2, my + h + 12))
                num_yticks = 8
                for i in range(num_yticks + 1):
                    pop = int(i * max_v / num_yticks)
                    # Patch: avoid division by zero for max_v
                    y = my + h - (pop / max_v) * h if max_v else my + h
                    pygame.draw.line(self.screen, WHITE, (mx - 8, y), (mx, y), 2)
                    label = font.render(str(pop), True, WHITE)
                    self.screen.blit(label, (mx - label.get_width() - 12, y - label.get_height() // 2))

                if len(times) > 1:
                    for i in range(1, len(times)):
                        x1 = mx + times[i-1] / max_t * w if max_t else mx
                        y1 = my + h - (total_vals[i-1] / max_v) * h if max_v else my + h
                        x2 = mx + times[i] / max_t * w if max_t else mx
                        y2 = my + h - (total_vals[i] / max_v) * h if max_v else my + h
                        pygame.draw.line(self.screen, WHITE, (x1, y1), (x2, y2), 2)
                    label = font.render("Total", True, WHITE)
                    self.screen.blit(label, (mx + w - label.get_width(), my + 8))

                for sp in self.sim.species:
                    color = sp.color
                    history = sp.history
                    if len(history) > 1:
                        for i in range(1, len(history)):
                            x1 = mx + times[i-1] / max_t * w if max_t else mx
                            y1 = my + h - (history[i-1] / max_v) * h if max_v else my + h
                            x2 = mx + times[i] / max_t * w if max_t else mx
                            y2 = my + h - (history[i] / max_v) * h if max_v else my + h
                            pygame.draw.line(self.screen, color, (x1, y1), (x2, y2), 2)
                        label = font.render(sp.name, True, color)
                        self.screen.blit(label, (mx + w - label.get_width(), my + 28 + 22 * self.sim.species.index(sp)))
            elif self.tab == 'stats':
                self.draw_stats_tab()
            elif self.tab == 'species':
                self.draw_species_tab()
                for ed in self.editors:
                    ed.draw(self.screen)

            for b in self.buttons:
                b.draw(self.screen)
            for t, b in self.tabs.items():
                b.draw(self.screen)
            if self.tab != 'species':
                for s in self.sliders:
                    s.draw(self.screen)

            # --- Draw species attribute descriptions dialog if toggled ---
            if self.show_species_desc:
                pygame.draw.rect(self.screen, (50, 50, 50, 200), self.species_desc_dialog_rect)
                pygame.draw.rect(self.screen, WHITE, self.species_desc_dialog_rect, 2)
                title_surf = FONT.render("Species Attribute Descriptions", True, WHITE)
                self.screen.blit(title_surf, (self.species_desc_dialog_rect.x + 10, self.species_desc_dialog_rect.y + 10))
                # Draw each description line
                for idx, (field, desc) in enumerate(self.descriptions.items()):
                    field_abbrev = self.ABBREV_FIELDNAMES.get(field, field)
                    text = FONT.render(f"{field_abbrev}: {desc}", True, WHITE)
                    self.screen.blit(text, (self.species_desc_dialog_rect.x + 10, self.species_desc_dialog_rect.y + 50 + idx * 30))

            pygame.display.flip()
            self.clock.tick(FPS)

    def draw_stats_tab(self):
        """Render stats tab with dynamic vertical scrolling."""
        total_count = self.sim.stats['total'][-1] if self.sim.stats['total'] else 0
        total_births = self.sim.stats['births'][-1] if self.sim.stats['births'] else 0
        total_deaths = self.sim.stats['deaths'][-1] if self.sim.stats['deaths'] else 0
        overall_births = sum(self.sim.stats['births']) if self.sim.stats['births'] else 0
        overall_deaths = sum(self.sim.stats['deaths']) if self.sim.stats['deaths'] else 0
        total_texts = [
            f"Total Count: {total_count}",
            f"Total Births (last step): {total_births}",
            f"Total Deaths (last step): {total_deaths}",
            f"Overall Births: {overall_births}",
            f"Overall Deaths: {overall_deaths}",
            f"Step: {self.sim.step_idx}"
        ]
        stat_labels = [
            "Count",
            "Births (last step)",
            "Deaths (last step)",
            "Overall Births",
            "Overall Deaths"
        ]
        margin = 24
        start_x = 50
        start_y = 100
        row_h = 28
        # Calculate visible rows
        total_rows = len(self.sim.species) * (len(stat_labels) + 1) + len(total_texts)
        max_rows = self.stats_max_visible_rows
        row_offset = self.stats_scroll_y
        # Draw total stats
        for i, txt in enumerate(total_texts):
            if row_offset <= i < row_offset + max_rows:
                self.screen.blit(FONT.render(txt, True, WHITE), (start_x, start_y + (i-row_offset)*row_h))
        # Draw species stats
        row_idx = len(total_texts)
        for idx, sp in enumerate(self.sim.species):
            if row_offset <= row_idx < row_offset + max_rows:
                name_label = FONT.render(sp.name, True, sp.color)
                self.screen.blit(name_label, (start_x, start_y + (row_idx-row_offset)*row_h))
            row_idx += 1
            count = len(sp.members)
            # Use tracked births/deaths per step
            births = 0
            deaths = 0
            if self.sim.species_step_stats and self.sim.step_idx > 0:
                births = self.sim.species_step_stats[self.sim.step_idx-1]['births'][idx]
                deaths = self.sim.species_step_stats[self.sim.step_idx-1]['deaths'][idx]
            overall_births = sp.total_births
            overall_deaths = sp.total_deaths
            stat_values = [count, births, deaths, overall_births, overall_deaths]
            for stat_idx, (label, value) in enumerate(zip(stat_labels, stat_values)):
                if row_offset <= row_idx < row_offset + max_rows:
                    stat_label = FONT.render(f"{label}: {value}", True, sp.color)
                    self.screen.blit(stat_label, (start_x, start_y + (row_idx-row_offset)*row_h))
                row_idx += 1
        # Draw vertical scrollbar if needed
        if total_rows > max_rows:
            bar_h = max_rows * row_h
            bar_x = start_x + 400
            bar_y = start_y
            pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, self.stats_scrollbar_thickness, bar_h), border_radius=6)
            thumb_h = max(bar_h * max_rows / total_rows, 40)
            thumb_y = bar_y + (row_offset / (total_rows - max_rows)) * (bar_h - thumb_h) if total_rows > max_rows else bar_y
            pygame.draw.rect(self.screen, WHITE, (bar_x+2, thumb_y, self.stats_scrollbar_thickness-4, thumb_h), border_radius=6)

    def toggle_species_desc(self):
        self.show_species_desc = not self.show_species_desc
        self.species_desc_btn.text = 'Hide Desc' if self.show_species_desc else 'Show Desc'
        # Reposition button in case window size changed
        self.species_desc_btn.rect.x = self.w - 200
        self.species_desc_btn.rect.y = self.h - 60

    def handle_mouse_species(self, pos):
        """Select species table cell on mouse click, with scroll offsets."""
        if self.editors:
            return
        start_x, start_y = 50, 100
        cell_w, cell_h = 120, 30
        cols = self.dm.fieldnames
        col_offset = self.species_scroll_x
        row_offset = self.species_scroll_y
        max_cols = self.species_max_visible_cols
        max_rows = self.species_max_visible_rows
        for i in range(row_offset, min(row_offset + max_rows, len(self.dm.species))):
            for j in range(col_offset, min(col_offset + max_cols, len(cols))):
                rect = pygame.Rect(start_x + (j-col_offset)*cell_w, start_y + (i-row_offset)*cell_h, cell_w, cell_h)
                if rect.collidepoint(pos):
                    self.selected_row, self.selected_col = i, j
                    return

    def import_data(self):
        """Import simulation stats from CSV and update simulation state."""
        self.sim.import_stats()

# Move this block outside the App class
if __name__ == '__main__':
    App().run()