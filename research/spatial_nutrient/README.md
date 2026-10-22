# Spatial nutrient integrity experiment

This source-only research workflow implements the proposed nutrient-only legacy-pg comparison. Read [the governing scientific contract](../../DOCUMENTATION.md#108-spatial-nutrient-integrity-experiment) and [protocol_v02.json](protocol_v02.json). The code does not imply that qualification passed or that main outcomes exist. See `EXECUTION_REPORT.md` for retained status.

## Setup and stages

From the repository root, using the standard Python 3.14.7 environment:

```bash
python -m pip install -r requirements.txt -r research/spatial_nutrient/requirements-spatial.txt
python -m research.spatial_nutrient.qualification_run tests --output OUTPUT/engineering.xml
python -m research.spatial_nutrient.runner register OUTPUT/pilot --cohort pilot
python -m research.spatial_nutrient.runner run OUTPUT/pilot --workers 4
python -m research.spatial_nutrient.runner analyze OUTPUT/pilot
```

`OUTPUT` is a user-chosen parent directory. Registration creates a new directory and rejects overwrites. All 20 pilot seeds are registered before execution; do not pass a shortened seed list as a pilot. An arm runs in its own process and shared-global threads never evolve simulations. Main seeds are rejected for diagnostic or qualification cohorts.

Each attempted arm retains `started.json`, progress, event/trajectory logs, RNG state and result. A completion marker references the validated result; the result references hashes and sizes of all final output files. A failed arm is not extinction. Use `run --retry` to restart the identical saved inputs in a new attempt after investigating a transient failure. If a worker was killed, `runner recover OUTPUT/pilot --seed SEED --arm U` releases only a verified dead local worker's lock and retains a recovery record. Never delete failure records or replace seeds. A changed executable or runtime requires a newly registered, consistently rerun cohort.

## Numerical qualification and freeze

Only after pilot mechanism and accounting checks are satisfactory, create a fine protocol copy changing `clock.dt_model_seconds` to 0.05 and `clock.steps` to 12000; retain `numerical_gates.refined_dt=0.05`. All scientific parameters and physical feed times remain fixed.

```bash
python -m research.spatial_nutrient.runner register OUTPUT/coarse --cohort numerical_coarse
python -m research.spatial_nutrient.runner register OUTPUT/fine --cohort numerical_fine --protocol FINE_PROTOCOL.json
python -m research.spatial_nutrient.runner run OUTPUT/coarse --workers 4
python -m research.spatial_nutrient.runner run OUTPUT/fine --workers 4
python -m research.spatial_nutrient.qualification_run qualify --pilot OUTPUT/pilot --tests OUTPUT/engineering.xml --coarse OUTPUT/coarse --fine OUTPUT/fine --output OUTPUT/qualification.json
python -m research.spatial_nutrient.qualification_run freeze --evidence OUTPUT/qualification.json --output OUTPUT/freeze.json
python -m research.spatial_nutrient.runner register OUTPUT/main --cohort main --freeze OUTPUT/freeze.json
```

A pilot-only report can omit both `--coarse` and `--fine`; it will keep the numerical gate pending and refuse freeze. Numerical looks use registered seeds 2000 onward at cumulative sizes 50, 200 and 1000. Register a new complete cumulative numerical cohort with the exact corresponding `--seeds` list when required; never select cases by outcome. The confidence decision spends the error budget across the prespecified looks. An inconclusive interval remains inconclusive. Diagnose individual-arm survival, births/deaths, food balance and exposure as well as H–U differences.

No main registration succeeds without the exact qualified protocol, code, runtime, retained evidence and test identity. The 400 main pairs quantify Monte Carlo uncertainty conditional on the model; they do not establish biological validation.

## Sensitivity and interpretation

Make explicit versioned protocol copies and use distinct `sensitivity_NAME` cohorts/seeds for founder counts 1/10, feed multipliers 0.5/2, speed 0/80, equal focal uptake, capacity 10000, source contrast, maintenance aging and daughter-size policy. Keep all changes visible in the saved protocol; do not pool these outcomes with the primary estimate. Container reordering and neutral label splitting are software invariants tested separately, not an expected zero biological effect.

The lifecycle intentionally retains legacy reserve clipping: unpaid maintenance is reported and forgiven before uptake. It is not energy-deficit death. Low-lifespan focal cells may never accumulate enough to divide. Resource-field contrast after extinction cannot establish experienced exposure; inspect occupancy, uptake and residence alongside the field metrics. A pilot without births in either arm cohort blocks qualification of the intended birth-and-death mechanism.

The separate opt-in live-FBA `BiomassMaintenanceFbaLifecyclePolicy` uses explicit biomass-maintenance loss and a death threshold. It is covered by independent integration tests and is not substituted into this experiment.
