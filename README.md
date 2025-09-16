# Social Flanker Error (SFE)

## Project Goal
Quantify how **social observation** modulates performance-monitoring ERPs—error-related negativity (ERN) and correct-related negativity (CRN)—at **short (within-block)** and **long (between-block)** timescales using **single-trial mixed-effects modeling** in a two-session social flanker paradigm, and release a fully reproducible analysis workflow.

## Background & Design
ERN and CRN index performance monitoring following error and correct responses, respectively. Prior work shows the ERN is sensitive to the **motivational significance** of errors, which increases under **social observation**; however, most studies rely on **trial-averaged ERPs**, potentially obscuring meaningful temporal dynamics. In this project, participants completed a flanker task **twice** (once **under social observation**, once **alone**). We extracted **single-trial** ERN/CRN amplitudes and used **mixed-effects models** to test whether social observation alters their **trajectories** over short (within blocks) and long (across blocks) timescales. Social observation **selectively shaped short-timescale dynamics**: in observed blocks, **ERN magnitudes increased** across trials while **CRN magnitudes decreased**; these trends were not evident when participants performed the task alone. Over longer timescales, **both ERN and CRN declined** across blocks **regardless of social context**, consistent with a vigilance decrement. To our knowledge, this is the **first demonstration** that social observation influences performance-monitoring trajectories over short timescales. These results underscore the value of **single-trial, time-resolved** analyses (beyond trial averages) and lay the groundwork for testing whether social observation interacts with **individual differences in motivation/affect** to shape performance-monitoring dynamics.

## Roadmap
A visual roadmap of planned releases will live at `docs/roadmap.drawio` (edited via [diagrams.net](https://app.diagrams.net/) with GitHub integration so changes are tracked as commits).

## Contents
Below is a brief guide to what is currently published on `main`:

#### `code/MADE-EEG-preprocessing-pipeline/`
Preprocessing code for EEG using the **MADE** pipeline.  
**Citation:** Debnath, R., Buzzell, G. A., Morales, S., Bowers, M. E., Leach, S. C., & Fox, N. A. (2020). *The Maryland analysis of developmental EEG (MADE) pipeline*. **Psychophysiology, 57(6), e13580**.

#### `code/postprocessing/`
Scripts to tabulate MADE outputs into analysis-ready CSV spreadsheets (e.g., `2_trial_level_erp_compute_control_n200Andn100.m`).

#### `code/statistics/`
R scripts for primary analyses:
- `1_basic_behavior.R` — **Subject-level behavioral data** summaries.
- `2_LMM_trial_erp_400200..._4trialCountNewData_final.R` — **Main ERP models**: trial-level ERN/CRN as a function of observation condition, trial, and block.
- `4_LMM_trial_erp_400200..._forTestACC.R` and `4_LMM_trial_erp_400200..._forTestRT.R` — **Trial-level behavioral outcomes** as a function of accuracy, observation, block number, and trial number.

#### `materials/task/`
PsychoPy task code for the social flanker (see `flanker-basic-v5/`, task spreadsheets, and image assets).

> Watch for our first public release (tagged on `main`) with frozen code, data dictionaries, and figures.

## Work in Development
This `main` branch contains completed releases for this project. For all work-in-progress, please switch over to the `dev` branches.

## Contributors
| Role | Name |
| --- | --- |
| Co-first authors / Analysis | **Yanbin Niu**\*¹, **Kianoosh Hosseini**\*² ³ |
| Task development | **Andy Peña**², **Carlos Rodriguez**² |
| PI / Advisor | **George A. Buzzell**² ³ |

\*Equal contribution. Superscripts denote affiliations (¹²³) as listed in the associated manuscript or lab pages.

Learn more about us [here](https://www.ndclab.com/people).

## Contributing
If you are interested in contributing, please read our [CONTRIBUTING.md](CONTRIBUTING.md) file.
