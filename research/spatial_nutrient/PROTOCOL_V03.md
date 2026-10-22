# BacSandbox spatial nutrient experiment protocol

**SPATIAL-NUTRIENT-001 · Version 0.3 · October 22, 2026**

**Status: prespecified feasibility plan; not executed or frozen for main inference.** This revision supersedes v0.2 for future study planning. The v0.2 protocol, pilot inputs, results and qualification decision remain unchanged historical evidence. All new settings and acceptance thresholds below are Own-model design choices, not measured biological parameters.

## 1. Decision and supporting evidence

PR #8 was reviewed and squash-merged into `master` at `c98e8d283506e7ba33a6f22af9eb5438cc4a1058`. Its three CI workflows passed at the reviewed head. A focused review rerun passed 73 tests and 46 subtests. Merging the implementation does not qualify the experiment.

The completed v0.2 pilot used 20 matched pairs (seeds 0–19), 40 first-attempt completions, 600 model seconds and 15,300 pg external food per arm. Every arm had zero births, and every focal population was extinct at the endpoint. Focal extinction occurred at 26.5–38.4 seconds in U and 29.7–44.4 seconds in H. Five resident cells remained across four H arms; none remained in U.

The pilot risk difference was zero. Its degenerate bootstrap interval [0, 0] is not evidence of equivalence; the reported conservative 95% interval was approximately [−0.6074, 0.6074]. Matched-input and accounting checks passed. Numerical qualification was not run, freeze was refused, and main seeds remain unused.

The mean paired full-time food-field contrast was 0.4379. Much of that field existed after focal extinction. Neither this statistic nor the alive-community-weighted field statistic establishes focal-local exposure. The main design problem is therefore feasibility of the intended birth-and-death ecology, before estimating an H–U effect.

| Decision | Revision |
| --- | --- |
| Establish a viable reference ecology | Select from a fixed uniform-source feasibility grid, without H outcomes. |
| Separate dose from delivery schedule | Cross four feed rates with two pulse intervals; retain exact total mass within each rate. |
| Reduce lifecycle confounding | Disable maintenance aging and fix daughter size in the candidate baseline. |
| Prevent outcome-driven gatekeeping | Require reproduction in U; do not reject H solely for extinction or absent births. |
| Confirm before inference | Use fresh pilot seeds, then numerical qualification, then a frozen main cohort. |

The scientific question remains: **How does spatial redistribution of equal external nutrient input change finite-time persistence of a rare, weaker-uptake competitor under a declared synthetic model?** Benefit, harm and little detectable change are all admissible findings. No mouse-day mapping or natural-species validation is claimed.

<!-- pagebreak -->

## 2. Fixed model and comparison

The exact inherited settings are defined by `protocol_v02.json`; `study_plan_v03.json` identifies its canonical hash and lists every planned override. The primary engine remains `legacy_pg_experiment` with the Python backend. The live-FBA lifecycle is a separate research path.

| Component | Prespecified setting |
| --- | --- |
| Founders | 500 total: 5 focal and 495 resident; labels independent of area-uniform positions. |
| Founder state | Size 1; age 0; initial reserve 2.5 pg, independent of capacity. Save every realized attribute. |
| Resource competition | Focal consumption speed 0.5; resident 1.0; other species settings inherited. |
| Storage and fission | Capacity 100 pg per unit size cubed; threshold 0.8; existing fission cost and reserve partitioning. |
| Lifecycle overrides | `maintenance_aging=false`; `daughter_size_policy=fixed`; recycling remains enabled. |
| Space and motion | Dish center (0, 0), radius 600 micrometers; speed 40 micrometers/s. |
| Clock and capacity | 600 seconds; dt 0.1 seconds; 6,000 ticks; capacity 5,000 cells. |
| Environment | Temperature 37, pH 7, oxygen 0.21, carbon dioxide 0; no phage inflow. |
| Neutralized pathways | Existing disabled interaction pathways and neutral radial motility, washout and culling factors retained. |

Maintenance aging is distinct from age mortality, which was already disabled. Fixed daughter size and maintenance aging are changed together to define a simpler candidate ecology; this revision does not attribute an improvement to either change separately. Their effects will be examined in declared sensitivities.

Legacy starvation semantics are retained: maintenance is clipped to available reserve before uptake, unpaid demand is forgiven, and starvation is evaluated after uptake. A tiny positive uptake can rescue a cell despite unmet demand. Report demanded, paid and unmet maintenance; do not describe this law as energy-deficit death. Changing that law requires a separate version.

For the eventual comparison, U has two equal-mass nutrient nodes in each half-domain per event. H reflects one poor-half node into the rich half, producing a 3:1 source allocation. Orientation is random and shared within the pair. This is uniform source allocation over time, not a constant or perfectly mixed food field.

Both arms receive identical founder snapshots, feed times, substrate identities and total mass. Only declared nutrient coordinates differ. Use global stable-agent random scheduling, tick-start crowding and next-tick daughter updates. Feed before washout and uptake. Preserve separate named RNG streams, full matching records and event ledgers. Uptake, occupancy, reserve, births and deaths are outcomes, never matching targets.

<!-- pagebreak -->

## 3. Feasibility controls and candidate grid

Before the grid, add controlled fixtures for each species at lifespan multipliers 0.4, 1.0 and 1.6. Use seeds 6000–6019, size 1, reserve 2.5 pg, fixed daughters, maintenance aging off, no movement, and disabled environmental and age mortality. These isolated controls require a dedicated fixture harness; they are not two-species runner cohorts.

Each no-feed fixture must show reserve exhaustion followed by starvation before 600 seconds, with zero births and a closed ledger. Each abundant-local-feed fixture supplies an accessible substrate-0 node at 100 pg for each live cell before each tick, records replenishment as external input, and must show at least one division before 600 seconds with no starvation. Any failure blocks the grid pending a documented diagnosis. Do not weaken a failed control after seeing grid outcomes.

Register all eight candidates before executing any. Each uses 20 seeds, 4000–4019. The same seeds across candidates support comparison but are not additional independent replicates. Set `contrast_nodes=0`: both U and H labels receive the same uniform allocation. Retain both as replay checks, and count only the 20 U runs in candidate criteria.

| Rank / ID | Feed multiplier | Interval (s) | Mass/node (pg) | Rate (pg/s) | Total/arm (pg) |
| --- | --- | --- | --- | --- | --- |
| 1 / C01 | 1 | 1.0 | 6.375 | 25.5 | 15,300 |
| 2 / C02 | 1 | 0.1 | 0.6375 | 25.5 | 15,300 |
| 3 / C03 | 2 | 1.0 | 12.75 | 51 | 30,600 |
| 4 / C04 | 2 | 0.1 | 1.275 | 51 | 30,600 |
| 5 / C05 | 4 | 1.0 | 25.5 | 102 | 61,200 |
| 6 / C06 | 4 | 0.1 | 2.55 | 102 | 61,200 |
| 7 / C07 | 8 | 1.0 | 51 | 204 | 122,400 |
| 8 / C08 | 8 | 0.1 | 5.1 | 204 | 122,400 |

All events have four nodes. Deliver at time zero through the final event strictly before 600 seconds. Mass per node equals 6.375 times the multiplier times the interval in seconds. Thus dose and pulse interval vary independently in the design. Smaller pulses also change node size and capture geometry in this engine; the interval comparison is not a pure test of encounter frequency.

A candidate is eligible only when all 20 pairs complete, matching/replay and accounting checks pass, no capacity culling occurs, at least 10 of 20 U runs each contain both a focal birth and a resident birth, and 4–16 of 20 U focal populations survive to 600 seconds, inclusive. These thresholds choose a nondegenerate operating regime; they are not biological targets or inferential tests.

Run and retain the complete grid. Select the first eligible candidate in the fixed rank order, preferring lower dose and then the original 1-second interval. Never rank by an H–U effect. If none qualifies, stop and record that outcome; any new grid or lifecycle law requires a new amendment and new seeds.

<!-- pagebreak -->

## 4. Fresh pilot, exposure and qualification

After selection, record the selected ID, full resolved configuration and hashes before generating heterogeneous inputs. Use fresh pilot seeds 5000–5019, restoring `contrast_nodes=1` for the intended U 2:2 versus H 3:1 comparison. Retain all 40 arms at 600 seconds, including zeros and failures. Do not recycle the v0.2 seeds or calibration runs into this pilot or the main estimate.

The fresh pilot requires all arms complete, exact matching, closed food/reserve/population ledgers, no capacity culling, and at least 10 of 20 U runs each containing both a focal and resident birth. Do not impose a U survival-count interval on the holdout pilot. Failure stops progression; it does not authorize selecting a different candidate based on pilot H results.

H may legitimately have zero births, complete extinction, or universal survival. Those outcomes do not by themselves fail the mechanism gate. Report the full U/H birth, death and survival distributions regardless of gate status. The v0.2 gate requiring births across both arm cohorts must be amended and tested before v0.3 qualification can be claimed.

Retain the existing full-time field manipulation gate: mean paired H-minus-U normalized food contrast at least 0.05. Also report contrast over the fixed 0–25-second window, chosen from the historical pilot before its earliest focal extinction. This additional diagnostic has no acceptance threshold and must not be selected retrospectively for a stronger contrast.

Add focal-specific rich/poor residence and uptake summaries from event and trajectory records. Normalize uptake and residence comparisons by the corresponding focal cell-time; report zero-denominator summaries as undefined, not zero. Record the time support, number of living focal cells and censoring. Distinguish source allocation, realized food field, and experienced focal exposure. Do not require equal uptake across arms or a favorable focal-exposure direction to pass.

| Gate | Requirement and response |
| --- | --- |
| Engineering | Full required integrity suite on the executing code/runtime; failures block progression. |
| Accounting | Existing 1e-9 scale-relative tick tolerance; exact count reconciliation; pruned sources below 0.001. |
| Fresh pilot | Complete matched cohort; U reproduction; no cap culling; field manipulation as declared above. |
| Numerical | Paired dt 0.1 versus 0.05, identical physical feed schedule and endpoint; seeds 2000 onward. |
| Numerical decision | Cumulative looks at 50, 200, 1,000 pairs; declared conservative staged intervals for H–U and each arm, all within ±0.05. |
| Freeze | All gates pass; bind exact code, runtime, selected config, tests, calibration decision and fresh-pilot/numerical evidence. |

An inconclusive numerical interval remains inconclusive, including at the final look. Retain births/deaths, food residual and exposure diagnostics alongside endpoint equivalence. Never substitute a point estimate or degenerate bootstrap for the numerical decision.

<!-- pagebreak -->

## 5. Main analysis, sensitivities and records

Only a verified v0.3 freeze admits the 400 main pairs, seeds 1000–1399. The estimand is the H survival probability minus the U survival probability for the selected synthetic ecology, where survival means focal count greater than zero at 600 seconds. The selection makes the inference conditional on a viable uniform baseline; it does not estimate an average over all candidate ecologies.

Extinction time is the first completed tick with zero focal agents; surviving lineages are right-censored at the endpoint. Report all four paired survival cells, both arm counts and probabilities, the paired risk difference, 10,000 whole-pair percentile bootstrap resamples, the conservative interval and Wilson arm intervals. Retain the declared tail-window establishment outcomes at thresholds 1, 5 and 10 from v0.2. Always display the conservative interval when the bootstrap is degenerate. Report Monte Carlo uncertainty as conditional on the model, not as uncertainty about natural populations.

| Seed range | Role and separation |
| --- | --- |
| 0–19 | Completed v0.2 pilot; permanently historical. |
| 1000–1399 | Main; reserved and unused. |
| 2000–2999 | Numerical qualification; reserved and unused. |
| 4000–4019 | Eight uniform feasibility candidates; shared across candidates only. |
| 5000–5019 | One fresh pilot after selection. |
| 6000–6019 | Controlled fixture replicates; separate from inference. |
| 7000–8099 | Eleven one-factor sensitivity blocks, 100 pairs each. |

In fixed block order, use focal founders 1 and 10 (total founders stays 500); feed multipliers 0.5 and 2 relative to the selected rate; motility 0 and 80; focal consumption speed 1; capacity 10,000; maintenance aging on; variable daughter size; and recycling off. Each block changes only the named factor from the selected primary setting and uses the next contiguous 100 seeds. No pooling with main results and no post hoc replacement of the primary setting. Any added block, source contrast, alternative starvation law or joint aging/size change needs a separately registered amendment.

Register each cohort and its complete denominator before execution. Preserve code/runtime/configuration identities, full initial snapshots, feed inputs, RNG states, per-attempt records, event logs, common-time trajectories and final file hashes. Numerical failure is not extinction. Retry the same saved realization only after diagnosing an execution failure; retain all attempts. A scientific implementation change requires a new registration and consistent rerun of the affected cohort.

The selected configuration, candidate eligibility table and reasons for exclusions must be published even if no main experiment runs. Keep the original v0.2 report and evidence archive intact. New Drive status entries and source packets must identify v0.3 as planned, distinguish merge status from scientific qualification, and leave historical Zhou evidence unchanged. Notebook ingestion must be confirmed separately before it is reported as current.

<!-- pagebreak -->

## 6. Implementation boundary and live FBA

`study_plan_v03.json` is a study-plan document, not an executable protocol accepted by the current validator. The eight files in `feasibility_v03/` use the existing v0.2 execution schema and already-supported parameter switches. They are runnable diagnostic inputs, with uniform replay in both arm labels. Their metadata identifies this v0.3 amendment; their schema version stays 0.2 deliberately.

Before scientific execution, complete the isolated control harness and validate it. Before fresh-pilot qualification or a main freeze, implement and test the v0.3 study-plan/selection binding, fresh pilot registry (5000–5019), U-only reproduction gate, cap-culling gate and focal-specific exposure summaries. The current CLI's v0.2 qualification path must not be represented as satisfying those new requirements. No selected candidate, v0.3 run, freeze or main result exists at publication.

Register all diagnostic configurations using the existing runner, explicit diagnostic cohort names and the complete 4000–4019 seed list. Then run all candidates and collect both arm labels. Do not use the existing `pilot` shortcut for calibration. `feasibility_v03/README.md` provides the command and explains the replay denominator. Schema validation and feed arithmetic checks are not feasibility execution.

The merged opt-in `BiomassMaintenanceFbaLifecyclePolicy` supplies reduced live-FBA biomass growth, exponential maintenance loss, starvation at a declared biomass fraction, division, ledger fields, event export and rollback. Callers must specify initial biomass, maintenance rate and starvation fraction. These are Own-model parameters; the policy is neither ATP-maintenance stoichiometry nor a genome-scale calibrated model.

A spatial live-FBA protocol needs its own registration before execution. Required decisions include finite-domain boundary conditions, extracellular mesh and diffusion, unit-consistent bolus/mass mapping, uptake kinetics, founder biomass, maintenance/death thresholds and resource recycling/removal semantics. Required evidence includes no-feed starvation, fed division, mass and biomass balance, diffusion/solver and timestep convergence, matched physical inputs, experienced exposure, and the same independence of selection from H outcomes. Test the growth-then-maintenance splitting explicitly. Do not inherit qualification from the legacy pilot, continuous Zhou reference or isolated FBA fixtures, and do not pool engines.

## 7. Authoritative records

The governing contract is `DOCUMENTATION.md`, section 10.8 and Own-design reference [54]. The complete historical execution record is `research/spatial_nutrient/EXECUTION_REPORT.md`; retained pilot files are under `retained/pilot_2026-10-21/`. The archive index is `evidence_archive.json`.

Historical protocol canonical SHA256: `6bb53d83f368f4d9cb7539220a897bd7bf2642342e4100a0c987901b7fcffdff`.

Historical executable identity SHA256: `25fd9afab07e81b69f32946bb6426d75ca84c022ed14f2bf0d717d3d300fc4fb`.

The current machine-readable plan is `study_plan_v03.json`. Its companion diagnostics, this readable protocol and the publication commit together define the prospective amendment. The original pilot's qualification remains inconclusive; the new plan does not retroactively change that decision.
