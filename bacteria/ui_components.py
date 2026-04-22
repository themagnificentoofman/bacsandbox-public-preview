import pygame
from pygame.locals import *
import sys
import os

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'

from .constants import *
from . import constants

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
        txt = constants.FONT.render(self.text, True, WHITE)
        s.blit(txt, txt.get_rect(center=self.rect.center))

    def click(self, pos):
        if self.enabled and self.rect.collidepoint(pos):
            self.action()

class Slider:
    """
    Horizontal slider controlling an environment variable; draggable knob.
    Supports float (0-1) or integer ranges.
    """
    def __init__(self, rect, label, env_key, get_env, minval=0.0, maxval=1.0, is_int=False):
        self.get_env = get_env
        self.rect = pygame.Rect(rect)
        self.label = label
        self.key = env_key
        self.drag = False
        self.minval = minval
        self.maxval = maxval
        self.is_int = is_int

    def draw(self, s):
        pygame.draw.rect(s, SLIDER_COLOR, self.rect, border_radius=3)
        val = self.get_env()[self.key]
        norm = (val - self.minval) / (self.maxval - self.minval)
        x = self.rect.x + norm * self.rect.w
        pygame.draw.circle(s, WHITE, (int(x), self.rect.centery), 6)
        disp_val = int(val) if self.is_int else f"{val:.2f}"
        label = constants.FONT.render(f"{self.label}: {disp_val}", True, WHITE)
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
            self.get_env()[self.key] = val

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
        txt_before = constants.FONT.render(before, True, WHITE)
        x = self.rect.x + 5 + txt_before.get_width()
        y = self.rect.y + 5
        s.blit(constants.FONT.render(before + after, True, WHITE), (self.rect.x + 5, y))
        if self.active and self.caret_visible:
            caret_height = constants.FONT.get_height()
            pygame.draw.line(s, WHITE, (x, y), (x, y + caret_height), 2)

    def handle_mouse_down(self, pos):
        """Set caret position based on mouse click within the cell."""
        if not self.rect.collidepoint(pos):
            return
        x = pos[0] - self.rect.x - 5
        display_text = self.text
        for i in range(len(display_text) + 1):
            width = constants.FONT.size(display_text[:i])[0]
            if width > x:
                self.caret_pos = max(0, i - 1)
                break
        else:
            self.caret_pos = len(display_text)

