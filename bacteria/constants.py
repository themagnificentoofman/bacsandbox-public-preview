from pygame.locals import *
import pygame

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

