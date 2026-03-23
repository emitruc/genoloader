# GenoLoader — Counting and plotting functions

This script ```GENOLOADER_counting_plotting_functions.ipynb``` is a Python Jupyter Notebook which contains three functions for summarizing and visualizing genetic load from `.gt` files produced by the GenoLoader `write_gt` function. These functions compute per-individual derived allele counts stratified by SnpEff functional effect category, generate strip-plot figures, and estimate the Rxy statistic for comparing mutational load between population pairs.

---

## Table of Contents

- [Requirements](#requirements)
- [Input files](#input-files)
- [Functions](#functions)
  - [counts_plots](#counts_plots)
  - [counts_plots_allele_lowcov](#counts_plots_allele_lowcov)
  - [freqchange_stat](#freqchange_stat)

---

## Requirements

- Python ≥ 3.7
- [pandas](https://pandas.pydata.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)

---

## Input files

### `.gt` file
The tab-delimited genotype file produced by `write_gt` (see main GenoLoader README). Must contain columns: `flag`, `effect`, and one column per individual with derived-allele dosage values (`0`, `1`, `2`, or `nan`).

### `.gt.covRandom1` file
The low-coverage companion file produced by `write_gt` when `low_cov='YES'`. Same structure as the `.gt` file but with pseudo-haploid allele values (`0`, `1`, or `.`).

### Population map file (`pop_map`)
A tab-delimited plain-text file with no header, assigning each sample to a population or species label. One assignment per line:

```
ind1    species_A
ind2    species_A
ind3    species_B
ind4    species_B
```

---

## Functions

### `counts_plots`

Computes per-individual counts of homozygous ancestral (`0`), heterozygous (`1`), and homozygous derived (`2`) genotypes across all four SnpEff effect categories, for loci passing a user-defined polarization flag filter. Results are saved to a counts table and visualized as strip plots grouped by population.

```python
counts_plots(f_gt, flag, pop_map)
```

| Parameter | Type  | Description |
|-----------|-------|-------------|
| `f_gt`    | `str` | Path to the `.gt` file produced by `write_gt` |
| `flag`    | `str` | Substring to filter polarization flags (e.g. `'unfolded'` matches all unfolded flags; see [Polarization flag filtering](#polarization-flag-filtering)) |
| `pop_map` | `str` | Path to the population map file |

**Effect categories summarized:** `MODIFIER`, `LOW`, `MODERATE`, `HIGH`

**Output columns in counts table:** `MF_0`, `MF_1`, `MF_2` (MODIFIER), `LO_0`, `LO_1`, `LO_2` (LOW), `MD_0`, `MD_1`, `MD_2` (MODERATE), `HI_0`, `HI_1`, `HI_2` (HIGH)

> **Note:** Loci with missing genotype data (`nan`) in any of the individuals listed in `pop_map` are excluded before counting.

---

### `counts_plots_allele_lowCov`

Equivalent to `counts_plots` but designed for the pseudo-haploid `.gt.covRandom1` output. Counts ancestral (`0`) and derived (`1`) alleles per individual (no heterozygous class) across the four effect categories.

```python
counts_plots_allele_lowCov(f_gt, flag, pop1, pop2, pop_map)
```

| Parameter | Type   | Description |
|-----------|--------|-------------|
| `f_gt`    | `str`  | Path to the `.gt.covRandom1` file produced by `write_gt` |
| `flag`    | `str`  | Substring to filter polarization flags |
| `pop_map` | `str`  | Path to the population map file |

**Output columns in counts table:** `MF_0`, `MF_1` (MODIFIER), `LO_0`, `LO_1` (LOW), `MD_0`, `MD_1` (MODERATE), `HI_0`, `HI_1` (HIGH)

---

### `freqChange_stat`

Computes the **Rxy statistic** (Do et al. 2015) for each SnpEff effect category, quantifying the relative excess of derived alleles in population X compared to population Y. Rxy > 1 indicates a higher load of derived alleles in population X. Also reports the average derived allele frequency in each population. Results are displayed as an inline table and saved as a figure.

```python
freqChange_stat(f_gt, flag, popx, popy, nameX, nameY)
```

| Parameter | Type   | Description |
|-----------|--------|-------------|
| `f_gt`    | `str`  | Path to the `.gt` file produced by `write_gt` |
| `flag`    | `str`  | Substring to filter polarization flags |
| `popx`    | `list` | List of sample names for population X |
| `popy`    | `list` | List of sample names for population Y |
| `nameX`   | `str`  | Display label for population X |
| `nameY`   | `str`  | Display label for population Y |

**Rxy formula** (per effect category):

$$R_{xy} = \frac{\sum_i f_{x,i} \cdot (1 - f_{y,i})}{\sum_i f_{y,i} \cdot (1 - f_{x,i})}$$

where $f_{x,i}$ and $f_{y,i}$ are the derived allele frequencies at locus $i$ in populations X and Y respectively.

> **Note:** Loci with missing genotype data in any individual from `popx` or `popy` are excluded before computation.

---
