# v0.3 uniform feasibility diagnostics

These eight inputs use the existing **v0.2 execution schema**. Read [the protocol](../PROTOCOL_V03.md) and [study plan](../study_plan_v03.json) for the prospective design. They are not qualified pilot/main inputs. Complete the declared isolated controls before running this grid.

Register all eight candidates before running any. From the repository root (Bash), choose a new output parent:

```bash
for candidate in C01 C02 C03 C04 C05 C06 C07 C08; do
  python -m research.spatial_nutrient.runner register "OUTPUT/v03_${candidate}" \
    --protocol "research/spatial_nutrient/feasibility_v03/${candidate}.json" \
    --cohort "diagnostic_v03_${candidate}" --seeds {4000..4019}
done
for candidate in C01 C02 C03 C04 C05 C06 C07 C08; do
  python -m research.spatial_nutrient.runner run "OUTPUT/v03_${candidate}" --workers 4
  python -m research.spatial_nutrient.runner analyze "OUTPUT/v03_${candidate}"
done
```

Each candidate has 20 independent seed realizations and two identical uniform arm allocations (`contrast_nodes=0`). Retain the second arm as replay evidence; do not count it as another independent realization. The grid contains 160 candidate/seed combinations and 320 arm executions. Common seeds across candidates are not independent additional evidence.

The declared eligibility decision counts U runs with births in **both** populations, U endpoint survival, accounting and cap culls; the existing `analyze` command supplies only part of that decision. Derive the remaining criteria from retained event/result records and save the complete eligibility table. H–U contrasts here are uniform replay checks, not treatment effects. Run all candidates and select the first eligible by the fixed rank in the study plan.

No v0.3 feasibility executions, selection or freeze existed at publication. Fresh-pilot and freeze support require the implementation amendments listed in the plan; the current v0.2 qualifier cannot certify v0.3 by renaming a cohort.
