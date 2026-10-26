"""Bounded synthetic viewer for BacSandbox's existing legacy simulation engine.

This presentation wrapper does not modify scientific code or run a registered
experiment. All values below are illustrative, not fitted biological parameters.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--headless", action="store_true", help="Run without a display or real-time pacing")
    parser.add_argument("--seconds", type=float, default=20, help="Simulated duration, >0 and <=120 (default: 20)")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output", type=Path, help="New CSV path; existing files are never overwritten")
    args = parser.parse_args()
    if not math.isfinite(args.seconds) or not 0 < args.seconds <= 120:
        parser.error("--seconds must be finite, greater than 0, and at most 120")
    if args.output and args.output.exists():
        parser.error("--output already exists; choose a new filename")
    if args.headless:
        os.environ.setdefault("SDL_VIDEODRIVER", "dummy")
        os.environ.setdefault("SDL_AUDIODRIVER", "dummy")

    import pygame
    from bacteria.core_sim import BacteriaSpecies, Simulation
    from bacteria.constants import SUBSTRATE_COLORS

    # Two invented populations share lifecycle values and differ only in motion
    # speed. This is a visual example, not an experiment testing a hypothesis.
    common = dict(init=20, rep=0.04, aging=0, so=0, sf=1, ph=7,
                  food_capacity=100, metabolic_rate=0.1, consumption_speed=2,
                  t6ss_activity=0, bacteriocin_production=0, immunity_factor=0,
                  acetate_secretion=0, acetate_uptake=0, baranyi_q0=10,
                  persister_switch=0, hgt_competence=0,
                  aff_calorie=1, aff_glycan=1, aff_mdcf=1)
    populations = [
        BacteriaSpecies("Synthetic A", speed=30, color=(82, 205, 201), **common),
        BacteriaSpecies("Synthetic B", speed=60, color=(242, 153, 93), **common),
    ]
    environment = dict(o2=0.2, co2=0.003, ph=7, temperature=37,
                       drop_frequency=10, drop_range=(20, 50), phage_inflow=0,
                       diet_profile="standard", spawn_mode="clustered",
                       species_arrangement="mixed", spawn_cluster_count=3,
                       spawn_cluster_spread=0.06, spawn_lumen_bias=0.9)
    sim = Simulation(populations, environment, cap=120, seed=args.seed,
                     physics_backend="python", ecm_backend="python")
    dt = 1 / 60
    ticks = max(1, math.ceil(args.seconds / dt))
    sim.update(600, 600, dt)  # Establish geometry and seed founders while paused.
    sim.running = True
    screen = clock = font = title = None
    if not args.headless:
        pygame.display.init()
        pygame.font.init()
        screen = pygame.display.set_mode((960, 640))
        pygame.display.set_caption("BacSandbox | synthetic software demonstration")
        clock = pygame.time.Clock()
        font = pygame.font.Font(None, 24)
        title = pygame.font.Font(None, 32)
    rows = [dict(time_s=0, synthetic_a=len(populations[0].members),
                 synthetic_b=len(populations[1].members), births=0, deaths=0)]
    completed_ticks = 0
    paused = False
    cancelled = False
    try:
        while completed_ticks < ticks:
            if screen:
                for event in pygame.event.get():
                    if event.type == pygame.QUIT or (event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE):
                        cancelled = True
                    elif event.type == pygame.KEYDOWN and event.key == pygame.K_SPACE:
                        paused = not paused
                if cancelled:
                    break
            if not paused:
                step = min(dt, args.seconds - completed_ticks * dt)
                if step <= 0:
                    break
                sim.update(600, 600, step)
                completed_ticks += 1
                rows.append(dict(time_s=min(completed_ticks * dt, args.seconds),
                                 synthetic_a=len(populations[0].members),
                                 synthetic_b=len(populations[1].members),
                                 births=sum(sp.total_births for sp in populations),
                                 deaths=sum(sp.total_deaths for sp in populations)))
            if screen:
                screen.fill((15, 24, 35))
                center = (310, 310)
                scale = 276 / sim.petri_radius
                pygame.draw.circle(screen, (26, 42, 52), center, 276)
                pygame.draw.circle(screen, (80, 105, 118), center, 276, 1)
                def point(x, y):
                    return (round(center[0] + x * scale), round(center[1] + y * scale))
                for node in sim.food_nodes:
                    pygame.draw.circle(screen, SUBSTRATE_COLORS[node.substrate_id], point(node.x, node.y), 2)
                for population in populations:
                    for cell in population.members:
                        pygame.draw.circle(screen, population.color, point(cell.x, cell.y), 3)
                lines = [
                    ("BacSandbox", title, (238, 243, 247)),
                    ("Synthetic software demo", font, (190, 207, 217)),
                    (f"Time: {rows[-1]['time_s']:.1f} / {args.seconds:g} s", font, (238, 243, 247)),
                    (f"Synthetic A: {len(populations[0].members)}", font, populations[0].color),
                    (f"Synthetic B: {len(populations[1].members)}", font, populations[1].color),
                    (f"Births: {rows[-1]['births']}    Deaths: {rows[-1]['deaths']}", font, (238, 243, 247)),
                    (f"Seed: {args.seed} | fixed 1/60 s steps", font, (190, 207, 217)),
                    ("Python physics + ECM backends", font, (190, 207, 217)),
                    ("Space: pause/resume | Esc: quit", font, (190, 207, 217)),
                    ("PAUSED" if paused else "Running existing simulation code", font, (238, 243, 247)),
                ]
                for i, (text, face, color) in enumerate(lines):
                    screen.blit(face.render(text, True, color), (615, 56 + i * 35))
                screen.blit(font.render("Illustrative inputs; no biological validation claim. Cells enlarged for visibility.", True,
                                        (190, 207, 217)), (25, 605))
                pygame.display.flip()
                clock.tick(60)
    finally:
        if screen:
            pygame.quit()
    if args.output:
        # Exclusive creation also guards a file appearing after the initial check.
        with args.output.open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps({"kind": "synthetic_software_demo", "seed": args.seed,
                      "requested_seconds": args.seconds,
                      "completed_seconds": rows[-1]["time_s"],
                      "status": "cancelled" if cancelled else "completed",
                      "steps": completed_ticks,
                      "physics_backend": "python", "ecm_backend": "python",
                      "final_populations": {sp.name: len(sp.members) for sp in populations},
                      "births": rows[-1]["births"], "deaths": rows[-1]["deaths"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
