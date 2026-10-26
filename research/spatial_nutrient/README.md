# Spatial nutrient experiment infrastructure

This directory contains source for a headless, auditable simulation workflow. It supports registered configurations and seeds, isolated process execution, retained attempts and RNG state, event/trajectory accounting, deterministic replay, feasibility selection, numerical qualification and gated freeze/registration.

| Engineering component | Source |
| --- | --- |
| Registration, isolated execution and recovery | `runner.py` |
| Accounting and replay checks | `analysis.py`, `runner.py` |
| Lifecycle controls | `controls_v03.py` |
| Candidate selection and stopping rules | `study_v03.py` |
| Qualification and freeze contracts | `qualification_run.py`, `qualification_v03.py` |
| Explicit protocols and candidate configurations | `protocol_v02.json`, `study_plan_v03.json`, `feasibility_v03/` |
| Disposable verification fixtures | `test_*.py` |

The protocols are retained as inspectable scientific configuration, not as permission to launch a registered experiment. Qualification/freeze gates remain enforced by the implementation. The public portfolio demo in `examples/portfolio_demo.py` is a separate bounded synthetic viewer.

Private run evidence, checkpoints and standalone operational result reports are not distributed in this preview. Historical technical documentation describes its source revision; the preview does not track intermediate research findings or assert biological predictive accuracy. See the root README and PORTFOLIO.md for runnable entry points and architecture.
