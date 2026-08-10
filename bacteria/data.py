import os, csv, json
import sys
import os

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'

from .constants import *

class SpeciesDataManager:
    """
    Handles loading, saving, and managing species parameters
    stored in a CSV file.
    """
    def __init__(self, filename=DATA_FILE):
        self.filename = filename
        # Updated CSV columns: name, initial count, rates, O2 sens, food sens, speed, preferred pH, food system, interactions
        self.fieldnames = ['name', 'init', 'rep', 'aging', 'so', 'sf', 'speed', 'ph', 'food_capacity', 'metabolic_rate', 'consumption_speed', 'greediness', 't6ss_activity', 'bacteriocin_production', 'immunity_factor', 'rep_thresh', 'share_thresh', 'predation_thresh', 'base_size', 'sensor_radius', 'temp_min', 'temp_max', 'acetate_secretion', 'acetate_uptake', 'baranyi_q0', 'persister_switch', 'persister_enter_thresh', 'persister_exit_thresh', 'phage_resistance', 'hgt_competence', 'aff_calorie', 'aff_glycan', 'aff_mdcf']

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
            
            # Fill missing keys with defaults for backward compatibility
            defaults = {
                'food_capacity': '10',
                'metabolic_rate': '0.1',
                'consumption_speed': '1.0',
                'greediness': '0.5',
                't6ss_activity': '0.5',
                'bacteriocin_production': '0.5',
                'immunity_factor': '0.5',
                'rep_thresh': '0.8',
                'share_thresh': '0.6',
                'predation_thresh': '0.9',
                'base_size': '1.0',
                'sensor_radius': '200.0',
                'temp_min': '4.0',
                'temp_max': '47.0',
                'acetate_secretion': '0.0',
                'acetate_uptake': '0.0',
                'baranyi_q0': '0.05',
                'persister_switch': '0.0',
                'persister_enter_thresh': '0.15',
                'persister_exit_thresh': '0.5',
                'phage_resistance': '0.0',
                'hgt_competence': '0.0',
                # Missing legacy affinities default to dimensionless calorie-only uptake-capacity
                # multipliers for backward compatibility; see DOCUMENTATION.md section 10.4.
                'aff_calorie': '1.0',
                'aff_glycan': '0.0',
                'aff_mdcf': '0.0'
            }
            normalized_species = []
            for sp in self.species:
                if not sp.get('t6ss_activity') and sp.get('aggression'):
                    sp['t6ss_activity'] = sp.get('aggression')
                if not sp.get('immunity_factor') and sp.get('defense'):
                    sp['immunity_factor'] = sp.get('defense')
                normalized_species.append({col: sp.get(col) or defaults.get(col, '0') for col in self.fieldnames})
            self.species = normalized_species
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
        blank.update({
            'name': '',
            'food_capacity': '10',
            'metabolic_rate': '0.1',
            'consumption_speed': '1.0',
            't6ss_activity': '0.5',
            'bacteriocin_production': '0.5',
            'immunity_factor': '0.5',
            'rep_thresh': '0.8',
            'share_thresh': '0.6',
            'predation_thresh': '0.9',
            'base_size': '1.0',
            'sensor_radius': '200.0',
            'temp_min': '4.0',
            'temp_max': '47.0',
            'acetate_secretion': '0.0',
            'acetate_uptake': '0.0',
            'baranyi_q0': '0.05',
            'persister_switch': '0.0',
            'persister_enter_thresh': '0.15',
            'persister_exit_thresh': '0.5',
            'phage_resistance': '0.0',
            'hgt_competence': '0.0'
        })
        self.species.append(blank)
        self.save()

    def remove(self, index):
        """
        Remove species at given index (if valid) and save.
        """
        if 0 <= index < len(self.species):
            self.species.pop(index)
            self.save()


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
            'o2': env['o2'],
            'co2': env['co2'],
            'ph': env['ph'],
            'temperature': env.get('temperature', 37.0),
            'drop_frequency': env['drop_frequency'],
            'drop_range_min': env['drop_range_min'],
            'drop_range_max': env['drop_range_max'],
            'capacity': env['capacity'],
            # M5: persist free-phage seeding rate so a configured outbreak survives reload
            'phage_inflow': env.get('phage_inflow', 0),
            # Phase I: persist the active diet profile (substrate-class drop mix) across reloads.
            'diet_profile': env.get('diet_profile', DEFAULT_DIET_PROFILE),
            # Phase I: persist the RNG seed (null = nondeterministic) for reproducible baselines.
            'seed': env.get('seed', None)
        }, f)

