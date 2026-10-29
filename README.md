# fungal-highways-reanalysis

Statistical reanalysis of plasmid transfer between bacteria along fungal networks ("fungal highways").

## Question

Do fungal hyphae act as transport vectors that increase horizontal gene transfer between bacteria, and which hyphal traits go with higher transfer? The experiments counted transconjugants (bacteria that received the plasmid) along three fungi, *Botrytis cinerea*, *Cunninghamella elegans* and *Pythium ultimum*, along glass fibres as an abiotic vector, and in a negative control, each at three inoculum levels. A spatial simulation gives the matching model prediction.

## My role

I reanalysed the transconjugant counts, made the figures, and explored how hyphal traits relate to transfer frequency. I did not run the experiments or the simulations.

## Status

Manuscript in preparation with collaborators at UFZ Leipzig. The citation will be added here when it is published.

## What the script does

`hgt_reanalysis.R`, from top to bottom:

1. Fits a linear model of log10(transconjugants + 1) with transport vector and inoculum as main effects (compared by AIC with the interaction model), checks assumptions (Shapiro-Wilk, Levene), and runs a Type II ANOVA with Tukey post-hoc comparisons (`car`, `emmeans`, `multcompView`).
2. Plots transconjugants per vector and inoculum level (`figures/Figure_HGT_frequency`).
3. Plots the simulated transconjugants per total cells for each vector (`figures/Fig5_simulated_transconjugants`).
4. Exploratory: correlates mean transconjugants with hyphal traits (water contact angle, hyphal diameter, water film width, vector and water-film coverage) on log and linear scales (`figures/Figure_hyphal_traits_*`).

## How to run

Last run in October 2026 with R 4.6.0, ggplot2 4.0.3, plyr 1.8.9, dplyr 1.2.1, emmeans 2.0.3, multcomp 1.4.32, multcompView 0.1.12, car 3.1.5 and cowplot 1.2.0.

```r
install.packages(c("ggplot2", "plyr", "dplyr", "emmeans", "multcomp", "multcompView", "car", "cowplot"))
```

1. Put the three data files in `data/` (see Data availability). `data/README.md` describes their columns.
2. Set the working directory to the repository folder and run `hgt_reanalysis.R`. The figures are written to `figures/`.

## Data availability

The data are not included in this repository: the transconjugant counts come from experiments by Tom Berthold, the simulated transconjugants from the modelling, and the hyphal traits from Table 1 of the manuscript. They are available from the authors on request.

## License

The code is released under the MIT License (see `LICENSE`). The license does not cover the data.
