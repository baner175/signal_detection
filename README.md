# signal_detection

R code and numerical results for **compensator-based signal detection under an unknown background**. [Find the paper here](https://arxiv.org/abs/2605.20508)

The project studies statistical tests and estimators for detecting a signal component in physics data. The methods use a proposal background density and estimates the compensator statistic to estimate and test the signal fraction. Included scripts cover numerical simulations, power studies, generated figures, and an analysis of Fermi-LAT data mimics $\gamma$ -ray emission from a mixture of dark matter annihilation (the signal of interest) and the astrophysical background.

## Method overview

The observed physics sample is modeled as a mixture

$$
f(x; \eta) = (1 - \eta) f_b(x) + \eta f_s(x)
$$

where:

- $\eta$ is the signal fraction.
- $f_s$ is the density of the known or specified signal distribution.
- $f_b$ is the density of the unknown background distribution.
- $g$ or $g_\beta$ is a proposal background density which is either fully specified or has a parametric form, respectively.

The scripts compare proposal backgrounds and estimate the signal fraction and test statistic in simulation and data-analysis settings.

## Repository contents

| Path | Description |
|---|---|
| `power_simulation_functions.R` | Shared functions for running and summarizing power simulations. |
| `Numerical_example_Sec2.1.R` | Likelihood-ratio test simulation for Section 2.1. |
| `Numerical_study_Sec3.R` | Power study for Section 3 with an unknown background. |
| `Numerical_study_Sec4.R` | Power study for Section 4 without a separate background sample. |
| `Numerical_example_Sec2.1/` | Stored likelihood-ratio test results for Section 2.1. |
| `Numerical_example_Sec3/` | Stored power-study results for Section 3. |
| `Numerical_example_Sec4/` | Stored power-study results for Section 4. |
| `Fermi_LAT_*.R` | Fermi-LAT analyses with and without a background-only sample. |
| `Fermi_LAT_physics.txt` | Physics-sample data used in the Fermi-LAT analysis. |
| `Fermi_LAT_bkg_only.txt` | Background-only sample used in the Fermi-LAT analysis. |
| `Figure_4_left.R` through `Figure_9_right.R` | Scripts for reproducing the corresponding manuscript figures. |
| `power_curve_legend_wobkg.pdf` | Legend PDF for a power-curve figure. |

## Requirements

The project requires R and the following packages:

- `truncdist`
- `VGAM`
- `doSNOW`
- `parallel`
- `foreach`
- `knitr` (used to format the Fermi-LAT results)
- `latex2exp` (used for plotting purposes)

Install the commonly needed packages from R:

```r
install.packages(c("truncdist", "VGAM", "doSNOW", "foreach", "knitr", "latex2exp"))
```

Some scripts use 8 worker processes by default. Adjust `makeCluster(8)` to match the available cores.

## Running the simulations

Run scripts from the repository root because several scripts use relative paths and write output directories beneath the repository.

```bash
Rscript Numerical_example_Sec2.1.R
Rscript Numerical_study_Sec3.R
Rscript Numerical_study_Sec4.R
```

The scripts set smaller demonstration values by default:

- `Numerical_example_Sec2.1.R`: `B <- 1e2` and `n_samp <- 5e3`
- `Numerical_study_Sec3.R`: `B <- 1e2`
- `Numerical_study_Sec4.R`: `B <- 1e2`

For the full paper-scale experiments, increase `B` to `1e5` and change `n_phys` to the desired values (`50`, `100`, `250`, `500`, `1000`, or `2000`) where noted in the scripts. The full simulations may be time-consuming and require substantial CPU and memory resources.

## Fermi-LAT example

The Fermi-LAT scripts read the bundled text datasets:

```text
Fermi_LAT_physics.txt
Fermi_LAT_bkg_only.txt
```

Run the two analyses with:

```bash
Rscript Fermi_LAT_Analysis_wbkg.R
Rscript Fermi_LAT_Analysis_wobkg.R
```

The `wbkg` script include a background-only sample, while the `wobkg` script performs the proposed sensitivity analysis and runs a conservative test.

## Reproducing figures

Each `Figure_*.R` script corresponds to a figure in the associated manuscript. Run it from the repository root:

```bash
Rscript Figure_4_left.R
Rscript Figure_4_right.R
Rscript Figure_5_left.R
Rscript "Figure_5(right)_and_Figure_6.R"
```

Repeat as needed for the remaining `Figure_*.R` scripts. Figure scripts may depend on pre-generated simulation output in the `Numerical_example_*` directories.

## Output

Generated simulation files use names that record their parameters, for example:

```text
Numerical_example_Sec2.1/LRT_B(1e+05)_beta0(4)_n_samp(5000)_eta(0)_lambda(0).csv
Numerical_example_Sec3/WBKG__n_phys(250)_beta0(2)_eta(0.01)_B(1e+05)_bkg_to_phys(2).txt
```

Random seeds are set explicitly in the simulation scripts so that the analyses can be rerun consistently.

## Notes

- All core scientific code is written in R.
- The repository is organized around manuscript sections and figures rather than as an installable R package.
- A `set.seed()` value is used for reproducibility, but parallel execution and package/platform differences can still affect floating-point details.
- Stored simulation files document parameter combinations and can be used instead of rerunning every experiment.

## Reproducibility checklist

1. Clone or download the repository.
2. Open R or a terminal in the repository root.
3. Install the required R packages.
4. Inspect the `B`, sample-size, and parameter settings near the top of the selected script.
5. Run the script and check the generated files beneath the corresponding `Numerical_example_*` directory.
6. Run the relevant figure scripts to regenerate manuscript plots.
