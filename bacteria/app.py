import math
import re
import pygame, sys, os
from pygame.locals import *
import sys
import os
import csv

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'

from .constants import *
from . import constants
from .data import SpeciesDataManager, load_env_config, save_env_config
from .core_sim import Bacterium, BacteriaSpecies, Simulation, FoodNode
from .ui_components import Button, Slider, CellEditor

class App:
    """
    Main application class: sets up pygame, UI, simulation,
    and dispatches events to draw/update accordingly.
    """
    def __init__(self):
        pygame.init()
        constants.FONT = pygame.font.SysFont('Consolas', 20)
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
            'rep': 'Reproduction rate per second (0–1)',
            'aging': 'Aging rate multiplier (base age penalty scalar)',
            'so': 'O₂ sensitivity (-1 to 1; <0 poisoned by O₂, >0 needs O₂)',
            'sf': 'Food sensitivity (0–1)',
            'speed': 'Base speed (divides by size)',
            'ph': 'Preferred pH (0–14; strict survival bell curve)'
        }
        # Add new fields to species table
        self.descriptions.update({
            'food_capacity': 'Base max food (pg Glucose, scales w/ \u03bcm³)',
            'metabolic_rate': 'Base food used (pg/tick, scales w/ Kleiber)',
            'consumption_speed': 'Base eat speed (pg/tick, scales w/ \u03bcm²)',
            'greediness': '0.0–1.0 preference for rich food even if far',
            'aggression': 'Toxin power against other species (0–1)',
            'defense': 'Toxin resistance (0–1)',
            'rep_thresh': 'Min food capacity % to replicate (0–1)',
            'share_thresh': 'Min food capacity % to share (0–1)',
            'predation_thresh': 'Max food capacity % to predate (0–1)',
            'base_size': 'Base size (\u03bcm, default 1.0)',
            'sensor_radius': 'Food detection radius (\u03bcm, default 200.0)'
        })
        if 'food_capacity' not in self.dm.fieldnames:
            self.dm.fieldnames += ['food_capacity', 'metabolic_rate', 'consumption_speed', 'greediness', 'aggression', 'defense']
        if 'rep_thresh' not in self.dm.fieldnames:
            self.dm.fieldnames += ['rep_thresh', 'share_thresh']
        if 'predation_thresh' not in self.dm.fieldnames:
            self.dm.fieldnames += ['predation_thresh']
        if 'base_size' not in self.dm.fieldnames:
            self.dm.fieldnames += ['base_size', 'sensor_radius']

        # Undo and redo stacks for editing species table
        self.undo_stack = []
        self.redo_stack = []
        
        self.global_undo_stack = []
        self.global_redo_stack = []

        # Abbreviated fieldnames for species tab
        self.ABBREV_FIELDNAMES = {
            'name': 'Name',
            'init': 'Init',
            'rep': 'Rep',
            'aging': 'Aging',
            'so': 'O2',
            'sf': 'Food',
            'speed': 'Spd',
            'ph': 'pH',
            'food_capacity': 'Cap',
            'metabolic_rate': 'Meta',
            'consumption_speed': 'Eat',
            'greediness': 'Greed',
            'aggression': 'Atk',
            'defense': 'Def',
            'rep_thresh': 'RepT',
            'share_thresh': 'ShrT',
            'predation_thresh': 'PredT',
            'base_size': 'Size',
            'sensor_radius': 'SensR'
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

        # Globals descriptions
        self.globals_descriptions = {
            'o2': 'O₂ level (0.0 to 1.0)',
            'co2': 'CO₂ level (0.0 to 1.0)',
            'ph': 'Environment pH (0.0 to 14.0)',
            'capacity': 'Max bacteria pop. (≥ 0)',
            'drop_frequency': 'Food drop frequency (0 to 100)',
            'drop_range_min': 'Min food drop amount (≥ 1)',
            'drop_range_max': 'Max food drop amount (≥ drop_range_min)'
        }
        self.show_globals_desc = False
        self.globals_desc_btn = Button((self.w - 200, self.h - 60, 180, 30), 'Show Desc', self.toggle_globals_desc)
        self.globals_desc_dialog_rect = None

        # Load env config as default env
        loaded_env = load_env_config()
        if loaded_env:
            self.default_env = loaded_env
        else:
            self.default_env = {
                'o2': 0.8,
                'co2': 0.2,
                'ph': 7.0,
                'drop_frequency': 5,
                'drop_range_min': 2,
                'drop_range_max': 8,
                'capacity': 500
            }
        self.default_env['drop_range'] = (self.default_env['drop_range_min'], self.default_env['drop_range_max'])

        self.init_sim()
        # Create buttons, tabs, sliders
        # Use consistent horizontal spacing for buttons
        button_w, button_h, button_gap = 80, 30, 10
        button_y = 20
        button_xs = [20 + i * (button_w + button_gap) for i in range(6)]
        self.play_btn = Button((button_xs[0], button_y, button_w, button_h), 'Play', self.play)
        self.import_btn = Button((button_xs[1], button_y, button_w, button_h), 'Import', self.import_data)
        self.pause_btn = Button((button_xs[2], button_y, button_w, button_h), 'Pause', self.toggle_pause)
        self.export_btn = Button((button_xs[3], button_y, button_w, button_h), 'Export', lambda: self.sim.export(self.dm.species, self.default_env), enabled=False)
        self.add_btn = Button((button_xs[4], button_y, button_w, button_h), 'Add', self.add_species, enabled=False)
        self.remove_btn = Button((button_xs[5], button_y, button_w, button_h), 'Remove', self.remove_species, enabled=False)
        
        self.zoom_in_btn = Button((0, 0, 40, 40), '+', lambda: self.set_zoom(1.0))
        self.zoom_out_btn = Button((0, 0, 40, 40), '-', lambda: self.set_zoom(-1.0))

        self.buttons = [
            self.play_btn,
            self.import_btn,
            self.pause_btn,
            self.export_btn,
            self.add_btn,
            self.remove_btn
        ]

        self.tabs = {
            'grid': Button((560, 20, 80, 30), 'Grid', lambda: self.set_tab('grid')),
            'graph': Button((650, 20, 80, 30), 'Graph', lambda: self.set_tab('graph')),
            'stats': Button((740, 20, 80, 30), 'Stats', lambda: self.set_tab('stats')),
            'species': Button((830, 20, 80, 30), 'Species', lambda: self.set_tab('species')),
            'globals': Button((920, 20, 80, 30), 'Globals', lambda: self.set_tab('globals'))
        }
        self.tab = 'grid'
        # Align sliders to top right, flush with button text
        slider_x = self.w - 220  # 200 width + 20 margin
        slider_y_start = 60      # Increased from 20 to 60 for visibility
        slider_gap = 40
        
        # Make a deep copy for live environment
        self.active_env = self.default_env.copy()

        # Add new sliders for drop frequency and drop range
        self.sliders = [
            Slider((slider_x, slider_y_start + 0*slider_gap, 200, 10), 'O2', 'o2', lambda: self.active_env),
            Slider((slider_x, slider_y_start + 1*slider_gap, 200, 10), 'CO2', 'co2', lambda: self.active_env),
            Slider((slider_x, slider_y_start + 2*slider_gap, 200, 10), 'pH', 'ph', lambda: self.active_env, minval=0, maxval=14),
            Slider((slider_x, slider_y_start + 3*slider_gap, 200, 10), 'Drop Freq', 'drop_frequency', lambda: self.active_env, minval=0, maxval=100, is_int=True),
            Slider((slider_x, slider_y_start + 4*slider_gap, 200, 10), 'Drop Range Min', 'drop_range_min', lambda: self.active_env, minval=1, maxval=20, is_int=True),
            Slider((slider_x, slider_y_start + 5*slider_gap, 200, 10), 'Drop Range Max', 'drop_range_max', lambda: self.active_env, minval=1, maxval=50, is_int=True)
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
        self.zoom = 1.0
        self.zoom_target = 1.0
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
        half_w = self.w / (2 * self.zoom)
        half_h = self.h / (2 * self.zoom)
        min_x = -self.GRID_HALF + half_w
        max_x = self.GRID_HALF - half_w
        min_y = -self.GRID_HALF + half_h
        max_y = self.GRID_HALF - half_h
        
        # If the view is larger than the grid, center it
        if min_x > max_x:
            self.grid_center_x = 0
        else:
            self.grid_center_x = max(min(self.grid_center_x, max_x), min_x)
            
        if min_y > max_y:
            self.grid_center_y = 0
        else:
            self.grid_center_y = max(min(self.grid_center_y, max_y), min_y)

    def handle_grid_scroll(self, dx, dy):
        self.grid_center_x += dx / self.zoom
        self.grid_center_y += dy / self.zoom
        self.clamp_grid_center()

    def set_zoom(self, dz):
        # smooth zoom approach (exponential step)
        self.zoom_target = max(0.01, min(self.zoom_target * (1.2 ** dz), 20.0))

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
            self.w//2 - (self.grid_center_x + self.GRID_HALF) * self.zoom,
            self.h//2 - (self.grid_center_y + self.GRID_HALF) * self.zoom,
            self.GRID_SIZE * self.zoom, self.GRID_SIZE * self.zoom
        )
        pygame.draw.rect(self.screen, WHITE, border_rect, 2)

        # Optionally draw grid lines or bacteria here

        # Draw scrollbars
        self.draw_scrollbars()

        # Draw current center coordinates at bottom right
        coord_text = f"Center: ({int(self.grid_center_x)}, {int(self.grid_center_y)}) | Zoom: {self.zoom:.2f}x"
        txt = constants.FONT.render(coord_text, True, WHITE)
        self.screen.blit(txt, (self.w - txt.get_width() - 20, self.h - txt.get_height() - 20))
        
        # Zoom buttons above the coordinate text
        btn_y = self.h - txt.get_height() - 60
        btn_w, btn_h = 40, 30
        self.zoom_in_btn.rect = pygame.Rect(self.w - btn_w - 20, btn_y, btn_w, btn_h)
        self.zoom_out_btn.rect = pygame.Rect(self.w - 2*btn_w - 30, btn_y, btn_w, btn_h)
        self.zoom_in_btn.draw(self.screen)
        self.zoom_out_btn.draw(self.screen)

    def draw_scrollbars(self):
        # Horizontal scrollbar
        bar_w = self.w - self.scrollbar_thickness
        bar_x = self.scrollbar_thickness // 2
        bar_y = self.h - self.scrollbar_thickness
        pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, bar_w, self.scrollbar_thickness), border_radius=6)
        # Thumb position
        view_w = self.w / self.zoom
        thumb_w = max(bar_w * view_w / self.GRID_SIZE, 40)
        
        # Avoid division by zero when grid is smaller than screen
        if self.GRID_SIZE > view_w:
            thumb_x = bar_x + (self.grid_center_x + self.GRID_HALF - view_w/2) / (self.GRID_SIZE - view_w) * (bar_w - thumb_w)
        else:
            thumb_x = bar_x + (bar_w - thumb_w) / 2
        pygame.draw.rect(self.screen, WHITE, (thumb_x, bar_y+2, thumb_w, self.scrollbar_thickness-4), border_radius=6)

        # Vertical scrollbar
        bar_h = self.h - self.scrollbar_thickness
        bar_x = self.w - self.scrollbar_thickness
        bar_y = self.scrollbar_thickness // 2
        pygame.draw.rect(self.screen, (80,80,80), (bar_x, bar_y, self.scrollbar_thickness, bar_h), border_radius=6)
        view_h = self.h / self.zoom
        thumb_h = max(bar_h * view_h / self.GRID_SIZE, 40)
        
        if self.GRID_SIZE > view_h:
            thumb_y = bar_y + (self.grid_center_y + self.GRID_HALF - view_h/2) / (self.GRID_SIZE - view_h) * (bar_h - thumb_h)
        else:
            thumb_y = bar_y + (bar_h - thumb_h) / 2
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
        half_w = self.w / (2 * self.zoom)
        half_h = self.h / (2 * self.zoom)
        left = self.grid_center_x - half_w
        top = self.grid_center_y - half_h
        right = self.grid_center_x + half_w
        bottom = self.grid_center_y + half_h
        rect_x = x + int((left - grid_min) / grid_size * size)
        rect_y = y + int((top - grid_min) / grid_size * size)
        rect_w = max(2, int((right - left) / grid_size * size))
        rect_h = max(2, int((bottom - top) / grid_size * size))
        
        # Clamp rect to minimap bounds visually (optional, but prevents drawing outside if zoomed out too much)
        # Actually it's fine just drawing it 
        
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
                    int(d.get('init') or 0),
                    float(d.get('rep') or 0.0),
                    float(d.get('aging') or 0.0),
                    float(d.get('so') or 0.0),
                    float(d.get('sf') or 0.0),
                    float(d.get('speed') or 40),
                    float(d.get('ph') or 7),
                    cols[i % len(cols)],
                    float(d.get('food_capacity') or 10),
                    float(d.get('metabolic_rate') or 0.1),
                    float(d.get('consumption_speed') or 1.0),
                    float(d.get('greediness') or 0.5),
                    float(d.get('aggression') or 0.5),
                    float(d.get('defense') or 0.5),
                    float(d.get('rep_thresh') or 0.8),
                    float(d.get('share_thresh') or 0.6),
                    float(d.get('predation_thresh') or 0.9),
                    float(d.get('base_size') or 1.0),
                    float(d.get('sensor_radius') or 200.0)
                ))
            except Exception:
                continue
        if not spp:
            spp.append(BacteriaSpecies('Ecoli', 50, 0.1, 0.05, 1, 1, 40, 7, cols[0], 10, 0.1, 1.0, 0.5, 0.5, 0.5, 0.8, 0.6, 0.9, 1.0, 200.0))
            
        self.active_env = self.default_env.copy()
        self.sim = Simulation(spp, self.active_env, cap=self.active_env['capacity'], steps=500)

    def add_species(self):
        """Append new blank species and re-init simulation, with undo/redo support."""
        blank = {k: '0' for k in self.dm.fieldnames}
        blank.update({
            'name': '',
            'food_capacity': '10',
            'metabolic_rate': '0.1',
            'consumption_speed': '1.0',
            'greediness': '0.5',
            'aggression': '0.5',
            'defense': '0.5',
            'rep_thresh': '0.8',
            'share_thresh': '0.6',
            'predation_thresh': '0.9',
            'base_size': '1.0',
            'sensor_radius': '200.0'
        })
        idx = len(self.dm.species)
        self.dm.species.append(blank)
        self.dm.save()
        self.init_sim()
        # Undo support
        self.undo_stack.append({'action': 'add_row', 'index': idx, 'row': blank})
        self.redo_stack.clear()

    def duplicate_species(self):
        """Duplicate currently selected species row and insert below it, with undo/redo support."""
        idx = self.selected_row if self.selected_row is not None else -1
        if 0 <= idx < len(self.dm.species):
            dup_row = self.dm.species[idx].copy()
            # If name is present, optionally change it? Keep it the same for now
            new_idx = idx + 1
            self.dm.species.insert(new_idx, dup_row)
            self.dm.save()
            self.init_sim()
            self.selected_row = new_idx
            self.undo_stack.append({'action': 'add_row', 'index': new_idx, 'row': dup_row})
            self.redo_stack.clear()

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

    def global_undo(self):
        """Undo last globals edit."""
        if not self.global_undo_stack:
            return
        action = self.global_undo_stack.pop()
        col, old, new = action['col'], action['old'], action['new']
        self.default_env[col] = old
        save_env_config(self.default_env)
        self.global_redo_stack.append({'action': 'cell_edit', 'col': col, 'old': old, 'new': new})

    def global_redo(self):
        """Redo last undone globals edit."""
        if not self.global_redo_stack:
            return
        action = self.global_redo_stack.pop()
        col, old, new = action['col'], action['old'], action['new']
        self.default_env[col] = new
        save_env_config(self.default_env)
        self.global_undo_stack.append({'action': 'cell_edit', 'col': col, 'old': old, 'new': new})

    def play(self):
        """Start simulation run if species exist and force unpause."""
        if not self.dm.species:
            print("Add at least one species before playing.")
            return
        self.init_sim()
        self.sim.running = True
        self.pause_btn.enabled = True
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
        
        if self.tab != 'species':
            self.show_species_desc = False
            self.add_btn.enabled = False
            self.remove_btn.enabled = False
        else:
            self.add_btn.enabled = True
            self.remove_btn.enabled = True
        
        if self.tab != 'globals':
            self.show_globals_desc = False

    def draw_globals_tab(self):
        start_x, start_y = 50, 100
        cell_w, cell_h = 240, 30
        env_keys = ['o2', 'co2', 'ph', 'capacity', 'drop_frequency', 'drop_range_min', 'drop_range_max']
        
        # Determine current target
        # Show default parameters, not live environment,
        # so it acts as "starting conditions"
        env_target = self.default_env
        
        for i, key in enumerate(env_keys):
            # Draw label
            label = constants.FONT.render(key, True, WHITE)
            self.screen.blit(label, (start_x, start_y + i*cell_h + 5))
            
            # Draw value cell
            rect = (start_x + 200, start_y + i*cell_h, cell_w, cell_h)
            highlight = (i == self.selected_row and self.selected_col == 0 and self.tab == 'globals')
            pygame.draw.rect(self.screen, SLIDER_COLOR if highlight else WHITE, rect, 2 if not highlight else 0)
            
            if highlight and not self.editors:
                overlay = pygame.Surface((cell_w,cell_h), pygame.SRCALPHA)
                overlay.fill((50,255,50,50))
                self.screen.blit(overlay, (rect[0], rect[1]))
            
            # Show edit context text if not being edited currently
            if not (self.editors and self.editors[0].row_idx == i and self.editors[0].col == key and self.tab == 'globals'):
                txt = constants.FONT.render(str(env_target.get(key, 0)), True, WHITE)
                self.screen.blit(txt, (rect[0] + 5, rect[1] + 5))

        # Draw show/hide descriptions button in bottom right
        self.globals_desc_btn.rect.x = self.w - 200
        self.globals_desc_btn.rect.y = self.h - 60
        # Ensure button text always matches dialog state
        self.globals_desc_btn.text = 'Hide Desc' if self.show_globals_desc else 'Show Desc'
        self.globals_desc_btn.draw(self.screen)
        # Draw descriptions dialog if shown, above and to the right of button
        if self.show_globals_desc:
            # Dynamically calculate dialog dimensions
            dialog_w = 0
            dialog_h = 40 + len(env_keys) * 32
            render_cache = []
            for field in env_keys:
                desc = self.globals_descriptions.get(field, '')
                label = f"{field}: {desc}"
                text = constants.FONT.render(label, True, WHITE)
                render_cache.append(text)
                dialog_w = max(dialog_w, text.get_width() + 40)
                
            dialog_x = self.globals_desc_btn.rect.right - dialog_w
            dialog_y = self.globals_desc_btn.rect.top - dialog_h - 10
            
            # Clamp dialog to window bounds
            if dialog_x < 10:
                dialog_x = 10
            if dialog_x + dialog_w > self.w - 10:
                dialog_x = self.w - dialog_w - 10
            if dialog_y < 10:
                dialog_y = 10
                
            self.globals_desc_dialog_rect = pygame.Rect(dialog_x, dialog_y, dialog_w, dialog_h)
            pygame.draw.rect(self.screen, (40,40,40), self.globals_desc_dialog_rect, border_radius=8)
            pygame.draw.rect(self.screen, WHITE, self.globals_desc_dialog_rect, 2)
            
            x, y = self.globals_desc_dialog_rect.x + 20, self.globals_desc_dialog_rect.y + 20
            for k, text in enumerate(render_cache):
                self.screen.blit(text, (x, y + k*32))

    def handle_globals_keys(self, e):
        if self.editors:
            return
        keys = ['o2', 'co2', 'ph', 'capacity', 'drop_frequency', 'drop_range_min', 'drop_range_max']
        if e.key == K_DOWN and self.selected_row < len(keys)-1:
            self.selected_row += 1
        if e.key == K_UP and self.selected_row > 0:
            self.selected_row -= 1
        if e.key == K_RETURN:
            key = keys[self.selected_row]
            val = self.default_env.get(key, 0)
            rect = (50 + 200, 100 + self.selected_row*30, 240, 30)
            self.editors = [CellEditor(key, rect, val, self.selected_row, 0)]

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
            head = constants.FONT.render(abbrev.get(col, col), True, WHITE)
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
                    txt = constants.FONT.render(str(row[col]), True, WHITE)
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
            # Dynamically calculate dialog dimensions
            dialog_w = 0
            dialog_h = 40 + len(cols) * 32
            render_cache = []
            for field in cols:
                desc = self.descriptions.get(field, '')
                label = f"{abbrev.get(field, field)}: {desc}"
                text = constants.FONT.render(label, True, WHITE)
                render_cache.append(text)
                dialog_w = max(dialog_w, text.get_width() + 40)
                
            dialog_x = self.species_desc_btn.rect.right - dialog_w
            dialog_y = self.species_desc_btn.rect.top - dialog_h - 10
            
            # Clamp dialog to window bounds
            if dialog_x < 10:
                dialog_x = 10
            if dialog_x + dialog_w > self.w - 10:
                dialog_x = self.w - dialog_w - 10
            if dialog_y < 10:
                dialog_y = 10
                
            self.species_desc_dialog_rect = pygame.Rect(dialog_x, dialog_y, dialog_w, dialog_h)
            pygame.draw.rect(self.screen, (40,40,40), self.species_desc_dialog_rect, border_radius=8)
            pygame.draw.rect(self.screen, WHITE, self.species_desc_dialog_rect, 2)
            
            x, y = self.species_desc_dialog_rect.x + 20, self.species_desc_dialog_rect.y + 20
            for k, text in enumerate(render_cache):
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
        if self.tab == 'grid':
            if self.zoom_in_btn.rect.collidepoint((mx, my)) or self.zoom_out_btn.rect.collidepoint((mx, my)):
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
            # Dynamically place sliders so they don't get truncated by right border
            slider_max_w = 200
            for s in self.sliders:
                # Use max possible logical value as a static worst-case width
                worst_case_val = f"{int(s.maxval)}" if s.is_int else f"{s.maxval:.2f}"
                txt_w = constants.FONT.render(f"{s.label}: {worst_case_val}", True, WHITE).get_width()
                slider_max_w = max(slider_max_w, txt_w)
                
            for s in self.sliders:
                s.rect.x = self.w - max(200, slider_max_w) - 30

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
                    if self.tab == 'grid':
                        self.zoom_in_btn.click(e.pos)
                        self.zoom_out_btn.click(e.pos)
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
                    if self.tab == 'globals':
                        if self.globals_desc_btn.rect.collidepoint(e.pos):
                            self.globals_desc_btn.click(e.pos)
                            continue
                        if self.editors:
                            self.editors[0].handle_mouse_down(e.pos)
                        else:
                            self.handle_mouse_globals(e.pos)
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
                    # Ctrl + + zoom in
                    if e.key in (K_PLUS, K_KP_PLUS, K_EQUALS) and pygame.key.get_mods() & KMOD_CTRL:
                        self.set_zoom(1.0)
                    # Ctrl + - zoom out 
                    if e.key in (K_MINUS, K_KP_MINUS) and pygame.key.get_mods() & KMOD_CTRL:
                        self.set_zoom(-1.0)

                    # Ctrl+Z undo
                    if e.key == K_z and pygame.key.get_mods() & KMOD_CTRL and not (pygame.key.get_mods() & KMOD_SHIFT):
                        if self.tab == 'species' and self.undo_stack and not self.editors:
                            self.undo()
                        elif self.tab == 'globals' and self.global_undo_stack and not self.editors:
                            self.global_undo()
                        continue
                    # Ctrl+Shift+Z redo
                    if e.key == K_z and pygame.key.get_mods() & KMOD_CTRL and pygame.key.get_mods() & KMOD_SHIFT:
                        if self.tab == 'species' and self.redo_stack and not self.editors:
                            self.redo()
                        elif self.tab == 'globals' and self.global_redo_stack and not self.editors:
                            self.global_redo()
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
                    # Duplicate species with Ctrl+D
                    if e.key == K_d and pygame.key.get_mods() & KMOD_CTRL and not (pygame.key.get_mods() & KMOD_SHIFT):
                        if self.tab == 'species' and not self.editors:
                            self.duplicate_species()
                        continue

                    # --- NEW LOGIC: Press Enter to start editing ---
                    if self.tab == 'species' and not self.editors and e.key == K_RETURN:
                        i, j = self.selected_row, self.selected_col
                        rect = (50 + j*120, 100 + i*30, 120, 30)
                        val = self.dm.species[i][self.dm.fieldnames[j]]
                        self.editors = [CellEditor(self.dm.fieldnames[j], rect, val, i, j)]
                        continue

                    if self.tab == 'globals' and not self.editors and e.key == K_RETURN:
                        self.handle_globals_keys(e)
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
                        if self.tab == 'globals':
                            col = ed.col
                            old_val = self.default_env.get(col, 0)
                            new_val = App.normalize_numeric_string(val)
                            try:
                                num = float(new_val)
                                if col in ['o2', 'co2']:
                                    if not 0.0 <= num <= 1.0: raise ValueError
                                elif col == 'ph':
                                    if not 0.0 <= num <= 14.0: raise ValueError
                                elif col == 'capacity':
                                    num = int(num)
                                    if num < 0: raise ValueError
                                elif col == 'drop_frequency':
                                    num = int(num)
                                    if not 0 <= num <= 100: raise ValueError
                                elif col == 'drop_range_min':
                                    num = int(num)
                                    if num < 1: raise ValueError
                                elif col == 'drop_range_max':
                                    num = int(num)
                                    if num < 1: raise ValueError
                                
                                if str(old_val) != str(num):
                                    self.global_undo_stack.append({'action': 'cell_edit', 'col': col, 'old': old_val, 'new': num})
                                    self.global_redo_stack.clear()
                                    self.default_env[col] = num
                                    save_env_config(self.default_env)
                            except Exception:
                                pass
                            self.editors = []
                        else:
                            i, j = ed.row_idx, ed.col_idx
                            col = ed.col
                            old_val = self.dm.species[i][col]
                            # --- Normalize numeric input for non-name columns ---
                            new_val = val
                            if col != 'name':
                                new_val = App.normalize_numeric_string(val)
                            if col == 'name':
                                if str(old_val) != str(new_val):
                                    self.undo_stack.append({'action': 'cell_edit', 'row': i, 'col': col, 'old': old_val, 'new': new_val})
                                    self.dm.species[i][col] = new_val
                                    self.redo_stack.clear()
                            else:
                                try:
                                    # Validation for new columns
                                    if col == 'init':
                                        num = int(new_val)
                                        if num < 0:
                                            raise ValueError
                                    elif col in ['rep', 'aging', 'sf']:
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
                                    elif col == 'food_capacity':
                                        num = float(new_val)
                                        if num < 0:
                                            raise ValueError
                                    elif col == 'metabolic_rate':
                                        num = float(new_val)
                                        if num < 0:
                                            raise ValueError
                                    elif col == 'consumption_speed':
                                        num = float(new_val)
                                        if num < 0:
                                            raise ValueError
                                    # apply change if valid
                                    if str(old_val) != str(new_val):
                                        self.undo_stack.append({'action': 'cell_edit', 'row': i, 'col': col, 'old': old_val, 'new': new_val})
                                        self.dm.species[i][col] = new_val
                                        self.redo_stack.clear()
                                except Exception:
                                    pass
                            self.dm.save()
                            self.editors = []

                if e.type == MOUSEWHEEL:
                    if self.tab == 'grid':
                        mods = pygame.key.get_mods()
                        if mods & KMOD_CTRL:
                            self.set_zoom(e.y)
                        else:
                            # Shift for horizontal scroll
                            if mods & KMOD_SHIFT:
                                self.handle_grid_scroll(-e.y * self.scroll_speed, 0)
                            else:
                                self.handle_grid_scroll(0, -e.y * self.scroll_speed)
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
            self.import_btn.enabled = not self.sim.running
            # --- Update drop range in env from sliders ---
            # Clamp drop_range_min and drop_range_max to maintain max >= min
            if self.active_env['drop_range_min'] > self.active_env['drop_range_max']:
                self.active_env['drop_range_max'] = self.active_env['drop_range_min']
            if self.active_env['drop_range_max'] < self.active_env['drop_range_min']:
                self.active_env['drop_range_min'] = self.active_env['drop_range_max']
            self.active_env['drop_range'] = (self.active_env['drop_range_min'], self.active_env['drop_range_max'])
            self.sim.drop_frequency = int(self.active_env['drop_frequency'])
            self.sim.drop_range = self.active_env['drop_range']
            # --- Always update simulation every frame ---
            self.sim.update(self.w, self.h)
            
            # Smooth zoom update
            if abs(self.zoom - self.zoom_target) > 0.001:
                self.zoom += (self.zoom_target - self.zoom) * 0.2
                self.clamp_grid_center()

            # Draw background
            self.screen.fill(BACKGROUND)
            # Draw content based on current tab
            if self.tab == 'grid':
                self.draw_grid()
                self.draw_minimap()
                # Draw food nodes in grid (Draw first so they are under bacteria)
                for node in self.sim.food_nodes:
                    sx = self.w//2 + (node.x - self.grid_center_x) * self.zoom
                    sy = self.h//2 + (node.y - self.grid_center_y) * self.zoom
                    if 0 <= sx < self.w and 0 <= sy < self.h:
                        # Scale food visual radius strictly by geometric cubic root (volume -> radius)
                        # amount = V (pg = um^3), r = cbrt(3V / 4pi)
                        food_r_um = max(0.5, ((3.0 * max(0.1, node.amount)) / (4.0 * math.pi)) ** (1.0 / 3.0))
                        radius = max(1, int(4 * self.zoom * food_r_um))
                        pygame.draw.circle(self.screen, (200,200,70), (int(sx), int(sy)), radius)

                # Draw bacteria in grid coordinates
                for sp in self.sim.species:
                    for b in sp.members:
                        # Convert grid coordinates to screen
                        sx = self.w//2 + (b.x - self.grid_center_x) * self.zoom
                        sy = self.h//2 + (b.y - self.grid_center_y) * self.zoom
                        # Only draw if visible
                        if 0 <= sx < self.w and 0 <= sy < self.h:
                            radius = max(1, int(4 * self.zoom * b.size))
                            pygame.draw.circle(self.screen, sp.color, (int(sx), int(sy)), radius)
                
                # --- HOVER TOOLTIP ---
                mx, my = pygame.mouse.get_pos()
                hovered_b = None
                hovered_food = None
                
                if not self.sim.running:
                    for sp in self.sim.species:
                        for b in sp.members:
                            sx = self.w//2 + (b.x - self.grid_center_x) * self.zoom
                            sy = self.h//2 + (b.y - self.grid_center_y) * self.zoom
                            radius = max(1, int(4 * self.zoom * b.size))
                            # Check distance to mouse
                            if (mx - sx)**2 + (my - sy)**2 <= max(radius, 5)**2:  # Use a min radius of 5 for easier hovering
                                hovered_b = b
                                break
                        if hovered_b:
                            break
                        
                if not hovered_b:
                    # Check food nodes if no bacterium is hovered
                    for node in self.sim.food_nodes:
                        sx = self.w//2 + (node.x - self.grid_center_x) * self.zoom
                        sy = self.h//2 + (node.y - self.grid_center_y) * self.zoom
                        food_r_um = max(0.5, ((3.0 * max(0.1, node.amount)) / (4.0 * math.pi)) ** (1.0 / 3.0))
                        radius = max(1, int(4 * self.zoom * food_r_um))
                        if (mx - sx)**2 + (my - sy)**2 <= max(radius, 5)**2:
                            hovered_food = node
                            break
                
                if hovered_b or hovered_food:
                    # Prepare text lines
                    if hovered_b:
                        lines = [
                            f"Species: {hovered_b.species.name}",
                            f"Age: {hovered_b.age:.1f}s",
                            f"Food: {hovered_b.food_level:.1f}/{hovered_b.food_capacity:.1f}",
                            f"Size: {hovered_b.size:.2f} (Base: {hovered_b.base_size:.2f})",
                            f"Metabolism: {hovered_b.metabolic_rate:.2f}/s",
                            f"Eating Spd: {hovered_b.consumption_speed:.2f}/s",
                            f"Speed: {hovered_b.speed:.2f} px/s",
                            f"Greediness: {hovered_b.species.greediness:.2f}",
                            f"Base Rep/s: {hovered_b.species.base_rep:.3f}",
                            f"Age Multiplier: {hovered_b.species.base_aging:.3f}"
                        ]
                        outline_color = hovered_b.species.color
                    else:
                        lines = [
                            f"Food Node",
                            f"Amount: {hovered_food.amount:.2f}",
                            f"Pos: ({int(hovered_food.x)}, {int(hovered_food.y)})"
                        ]
                        outline_color = (200, 200, 70)
                    
                    # Render texts to find dimensions
                    rendered_lines = [constants.FONT.render(text, True, WHITE) for text in lines]
                    box_w = max(r.get_width() for r in rendered_lines) + 20
                    box_h = len(rendered_lines) * 24 + 10
                    
                    # Dynamic position (clamp to screen edges)
                    box_x, box_y = mx + 15, my + 15
                    if box_x + box_w > self.w:
                        box_x = mx - box_w - 15
                    if box_y + box_h > self.h:
                        box_y = my - box_h - 15
                        
                    # Draw box backing + border
                    tooltip_rect = pygame.Rect(box_x, box_y, box_w, box_h)
                    overlay = pygame.Surface((box_w, box_h), pygame.SRCALPHA)
                    overlay.fill((30, 30, 30, 230))
                    self.screen.blit(overlay, (box_x, box_y))
                    pygame.draw.rect(self.screen, outline_color, tooltip_rect, 2, border_radius=4)
                    
                    # Blit text
                    curr_y = box_y + 10
                    for rline in rendered_lines:
                        self.screen.blit(rline, (box_x + 10, curr_y))
                        curr_y += 24

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

                font = constants.FONT
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
                    self.screen.blit(label, (mx + w + 10, my + 8))

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
                        self.screen.blit(label, (mx + w + 10, my + 30 + 22 * self.sim.species.index(sp)))
            elif self.tab == 'stats':
                self.draw_stats_tab()
            elif self.tab == 'species':
                self.draw_species_tab()
                for ed in self.editors:
                    ed.draw(self.screen)

            elif self.tab == 'globals':
                self.draw_globals_tab()
                for ed in self.editors:
                    ed.draw(self.screen)

            for b in self.buttons:
                b.draw(self.screen)
            for t, b in self.tabs.items():
                b.draw(self.screen)
            if self.tab not in ('species', 'globals'):
                for s in self.sliders:
                    s.draw(self.screen)

            # --- Draw species attribute descriptions dialog if toggled ---
            if self.show_species_desc:
                pygame.draw.rect(self.screen, (50, 50, 50, 200), self.species_desc_dialog_rect)
                pygame.draw.rect(self.screen, WHITE, self.species_desc_dialog_rect, 2)
                title_surf = constants.FONT.render("Species Attribute Descriptions", True, WHITE)
                self.screen.blit(title_surf, (self.species_desc_dialog_rect.x + 10, self.species_desc_dialog_rect.y + 10))
                # Draw each description line
                for idx, (field, desc) in enumerate(self.descriptions.items()):
                    field_abbrev = self.ABBREV_FIELDNAMES.get(field, field)
                    text = constants.FONT.render(f"{field_abbrev}: {desc}", True, WHITE)
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
                self.screen.blit(constants.FONT.render(txt, True, WHITE), (start_x, start_y + (i-row_offset)*row_h))
        # Draw species stats
        row_idx = len(total_texts)
        for idx, sp in enumerate(self.sim.species):
            if row_offset <= row_idx < row_offset + max_rows:
                name_label = constants.FONT.render(sp.name, True, sp.color)
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
                    stat_label = constants.FONT.render(f"{label}: {value}", True, sp.color)
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

    def toggle_globals_desc(self):
        self.show_globals_desc = not self.show_globals_desc
        self.globals_desc_btn.text = 'Hide Desc' if self.show_globals_desc else 'Show Desc'
        # Reposition button in case window size changed
        self.globals_desc_btn.rect.x = self.w - 200
        self.globals_desc_btn.rect.y = self.h - 60

    def toggle_species_desc(self):
        self.show_species_desc = not self.show_species_desc
        self.species_desc_btn.text = 'Hide Desc' if self.show_species_desc else 'Show Desc'
        # Reposition button in case window size changed
        self.species_desc_btn.rect.x = self.w - 200
        self.species_desc_btn.rect.y = self.h - 60

    def handle_mouse_globals(self, pos):
        """Select globals table cell on mouse click."""
        if self.editors:
            return
        start_x, start_y = 50, 100
        cell_w, cell_h = 240, 30
        env_keys = ['o2', 'co2', 'ph', 'capacity', 'drop_frequency', 'drop_range_min', 'drop_range_max']
        
        for i in range(len(env_keys)):
            rect = pygame.Rect(start_x + 200, start_y + i*cell_h, cell_w, cell_h)
            if rect.collidepoint(pos):
                self.selected_row, self.selected_col = i, 0
                return

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
        if self.sim.running:
            print("Cannot import while simulation is running. Pause first.")
            return

        if not os.path.exists(EXPORT_FILE):
            print(f"File {EXPORT_FILE} not found.")
            return

        with open(EXPORT_FILE, 'r', newline='') as f:
            reader = csv.reader(f)
            lines = list(reader)

        env_config = None
        species_config = None

        for i, row in enumerate(lines):
            if not row: continue
            if row[0] == '--- Environment Configuration ---':
                env_config = {}
                j = i + 1
                while j < len(lines) and lines[j] and lines[j][0] not in ('--- Species Configuration ---', '--- Global Stats ---'):
                    if len(lines[j]) >= 2:
                        k, v = lines[j][0], lines[j][1]
                        try:
                            val = float(v)
                            if val.is_integer():
                                val = int(val)
                            env_config[k] = val
                        except ValueError:
                            if k == 'drop_range':
                                try:
                                    val = eval(v)
                                    env_config[k] = list(val) if isinstance(val, tuple) else val
                                except Exception:
                                    env_config[k] = v
                            else:
                                env_config[k] = v
                    j += 1
            elif row[0] == '--- Species Configuration ---':
                species_config = []
                j = i + 1
                if j < len(lines) and lines[j]:
                    headers = lines[j]
                    j += 1
                    while j < len(lines) and lines[j] and lines[j][0] not in ('--- Environment Configuration ---', '--- Global Stats ---'):
                        sp_dict = dict(zip(headers, lines[j]))
                        species_config.append(sp_dict)
                        j += 1

        if env_config is not None:
            self.default_env.update(env_config)
            
            env = self.default_env
            save_env_config(self.active_env)

        if species_config is not None:
            self.dm.species = species_config
            self.dm.save()

        if env_config is not None or species_config is not None:
            self.init_sim()

        # Pre-scan for old format back-compat if no species config was found
        if species_config is None:
            stats_start = 1
            for idx, row in enumerate(lines):
                if row and row[0] == 'time' and len(row) >= 4 and row[1] == 'total':
                    stats_start = idx + 1
                    break

            stats_end = stats_start
            while stats_end < len(lines) and lines[stats_end]:
                stats_end += 1
            species_start = stats_end + 2
            
            if species_start < len(lines) and lines[species_start] and lines[species_start][0] != '--- Species Configuration ---':
                imported_species = {}
                for row in lines[species_start:]:
                    if not row or row[0].startswith('---'): continue
                    if len(row) >= 8 and row[0].isdigit():
                        idx, name = int(row[0]), row[1]
                        if idx not in imported_species:
                            imported_species[idx] = name
                    elif len(row) >= 7 and not row[0].isdigit():
                        name = row[0]
                        if name not in imported_species.values():
                            imported_species[len(imported_species)] = name

                needs_init = False
                while len(self.dm.species) < len(imported_species):
                    self.dm.add_blank()
                    needs_init = True

                for idx, name in imported_species.items():
                    if idx < len(self.dm.species):
                        self.dm.species[idx]['name'] = name

                if needs_init:
                    self.dm.save()
                    self.init_sim()

        self.sim.import_stats()
        self.pause_btn.enabled = False
        print("Graph loaded! Click 'Play' to start a new simulation.")

