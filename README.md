Statistical Analysis of Molecular Docking-Derived Binding Energies

This repository contains the R script used for the statistical analysis of molecular docking-derived binding-energy values across four Janus kinase (JAK) isoforms: JAK1, JAK2, JAK3, and TYK2.

Study objective

The analysis compares the distributions of docking-derived binding-energy values among three predefined molecular groups:

Reference JAK inhibitors vs flavonoids

Other phytochemicals vs reference JAK inhibitors

Flavonoids vs other phytochemicals

The three comparisons are evaluated separately for JAK1, JAK2, JAK3, and TYK2.

Input data

The analysis contains the docking-derived values directly within the R script. The dataset contains 13 ligands:

Reference JAK inhibitors: Baricitinib, Peficitinib, Tofacitinib

Flavonoids: Catechin, Kaempferol, Naringenin, Quercetin

Other phytochemicals: Caffeic acid, Chlorogenic acid, Coumaric acid, Gallic acid, Gingerol, Isoliquiritigenin

The four columns JAK1, JAK2, JAK3, and TYK2 contain the corresponding docking-derived binding-energy values.

Statistical analysis

Because the analysis uses small group sizes and docking-derived measurements, non-parametric pairwise comparisons were performed using the exact two-sided Wilcoxon rank-sum test (Mann--Whitney U test).

For each comparison and JAK isoform, the script reports:

Sample size for each group

Median

Median difference

Wilcoxon W statistic

IQR value

Two-sided P value

Hodges--Lehmann difference estimate p value

There are 3 predefined group comparisons across 4 JAK isoforms, resulting in 12 statistical tests.

The Holm step-down procedure is applied across all 12 P values.

The adjusted significance threshold used by the script is:

P_adjusted < 0.05

Output

The script creates a results directory containing three separate CSV result tables:

results/ ├── Table_1_Reference_JAKi_vs_Flavonoids/ │ └── Table_1_Reference_JAKi_vs_Flavonoids.csv ├── Table_2_Other_phytochemicals_vs_Reference_JAKi/ │ └── Table_2_Other_phytochemicals_vs_Reference_JAKi.csv ├── Table_3_Flavonoids_vs_Other_phytochemicals/ │ └── Table_3_Flavonoids_vs_Other_phytochemicals.csv └── input_binding_affinity_data.csv

The input-data CSV is also exported to document the exact values used by the analysis.

Software and R packages

The script uses:

dplyr

purrr

readr

The script checks whether these packages are installed and installs missing packages from CRAN when necessary.

The exact R/RStudio version used for the publication should be recorded separately in the manuscript or repository metadata if required by the journal.

Reproducibility

To reproduce the analysis:

Install R.

Download or clone this repository.

Open the R script in RStudio or another R environment.

Run the complete script.

The three comparison tables and the input-data CSV will be generated automatically in the results directory.

Interpretation

The values analyzed here are molecular docking-derived binding-energy estimates. They are computational estimates and should not be interpreted as experimentally measured binding affinities.

The statistical analysis evaluates differences between the predefined molecular groups in the available docking-derived dataset.

Repository contents

Statistical-Analysis/ ├── README.md └── RStudio_statistical_analysis.R

Citation

If the code or data from this repository are used, please cite the associated publication.

DOI
https://doi.org/10.5281/zenodo.23051817 

Authors
Kriti Gupta, Richa Ashma, Manjari Jonnalagadda
