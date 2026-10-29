# Data files

The data are not included in this repository (see the main README). `hgt_reanalysis.R` expects these three files here:

**`Testdatensatz_dreiExperimente_noZero.csv`**: transconjugant counts. Semicolon-separated, decimal comma.
Columns: `Vector` (negativecontrol, botrytis, cunning, pythium, glassfibre), `Inoculum` (low, medium, high), `Replicate`, `Cellcount` (transconjugants), `Sampling`, `Well`, `Hyphaecoverage`, `Waterfilmcoverage`.

**`fig5_data.csv`**: simulated transconjugants. Comma-separated.
Columns: `Fungi` (B, C, P, G for *B. cinerea*, *C. elegans*, *P. ultimum*, glass fibre), `KON`, `minimalDiff`, `Replicate`, `TotalPop`, `TotalHyph`, `TotalTrans`, `TransPerHyph`, `TransPerPop`.

**`hyphal_traits.csv`**: hyphal traits from Table 1 of the manuscript, one row per vector (botrytis, cunning, pythium, glassfibre). Comma-separated.
Columns: `Vector`, `contact_angle`, `hyphal_diameter`, `water_film_width`, `vector_coverage`, `waterfilm_coverage`, `fractal_dimension`, `lacunarity` (`NA` where a trait does not apply).
