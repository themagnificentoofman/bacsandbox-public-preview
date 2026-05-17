import sys
import os

if __package__ is None and not hasattr(sys, 'frozen'):
    # Allows execution via VS Code's "Run" button directly
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    from bacteria.app import App
else:
    from .app import App

if __name__ == '__main__':
    App().run()