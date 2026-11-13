# Fish pipeline comparisons

## Overview

This repository contains a set of bioinformatics workflows and analysis notebooks developed to compare alternative sequence inference and taxonomic assignment approaches for freshwater fish environmental DNA (eDNA) metabarcoding datasets. The project evaluates how different sample inference and taxonomic classification methods influence biodiversity assessments across multiple fish datasets. Four datasets are included: Marchamley, Loch Insh (labelled as Ring Trial), Windermere and simulated mock dataset. 

The analysis compares outputs from three sequence inference pipelines (MetaBEAT, DADA2 and VSEARCH) and five taxonomic assignment methods (MetaBEAT, BLAST + LCA, RDP, VSEARCH + LCA, MAPseq). Visual data from the three real ecologial datasets also enables assessment of concordance between bioinformatic approaches and known species records.

## Repository workflow

The project is organised into seven main notebooks:

- `00_reference_database_formatting`: prepares and formats reference databases for use by the different bioinformatics tools.
- `01_run_sequence_inference`: runs the sequence inference steps of the tools/pipelines included in the comparison.
- `02_run_assign_taxonomy`: runs the taxonomy assignment steps of the tools/pipelines included in the comparison.
- `03_run_groundtruth_simulated`: processes simulated ground-truth datasets.
- `04_data_formatting`: formats and standardises outputs for downstream analysis into standardised long-format datasets and phyloseq objects.
- `05_compare_tools`: compares results between the different bioinformatics tools.
- `06_explore_simulated`: explores results obtained from the simulated datasets.

The repository also contains `Data`, `Results`, and `Scripts` folders, which store the input data, analysis outputs, and supporting scripts used throughout the main analysis notebooks.

## Outputs

The repository automatically produces:

- ASV abundance tables
- Taxonomic assignment tables
- Long-format integrated datasets
- phyloseq objects
- Comparative analyses and visualisations (for internal EA report)
- Concordance assessments against historical fish records