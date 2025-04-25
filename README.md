# PRISTINE: PRimer-based Identification Suite Targeting Informative Nucleotide Elements

## Installation

This pipeline is distributed as a self-contained [Apptainer](https://apptainer.org/) (formerly known as Singularity) container, which includes all necessary dependencies and tools pre-installed.

To use the pipeline, you will need to have **Apptainer** installed on your system. Installation instructions are available on the [Apptainer documentation site](https://apptainer.org/docs/).

Once Apptainer is installed, no further setup is required. Simply download the pipeline container available at https://nextcloud.uni-greifswald.de/index.php/s/CFPAgy82SSLicqN and run the pipeline using:

```bash
./primer_design_tool_full.sif
```

This single .sif file encapsulates the entire environment, ensuring consistent and reproducible results across different systems.

## Input

To run the pipeline, you need:

- A configuration file named `config.yaml`
- A properly structured input directory with raw genome assemblies

### Configuration File

The file `config.yaml` must be placed in the **same directory as the container** (`primer_design_tool_full.sif`). It defines all global parameters, paths, and settings required for the analysis (for details, see below).

## Quick Start with Toy Dataset

A toy dataset is provided in the `toy_dataset/` directory. It includes:

- A small set of raw genome FASTA files organized in the correct input structure
- A ready-to-use configuration file: `config.yaml`

To quickly test the pipeline:

1. **Download the `toy_dataset/` directory**
2. **Place the Apptainer container** (`primer_design_tool_full.sif`) inside the `toy_dataset/` folder
3. From within the `toy_dataset/` directory, run:

```bash
./primer_design_tool_full.sif
```
---

## Pipeline overview

![Pipeline Diagram](pipeline.svg)

### Input Directory Structure

Currently, the pipeline supports input in the form of **raw genome FASTA files** (annotation-based support coming soon).

The input directory is organized by species. Each target species must have its own folder, named in the format Genus_species, where:

Genus is the genus name (capitalized)

species is the species name (lowercase)

Inside each species folder:

There must be one or more genome files in FASTA format

Each FASTA file must be named exactly as the folder, with a suffix _i, where i is a unique integer (e.g., Yersinia_pseudotuberculosis_1.fasta)

In addition, each species folder contains a subfolder named non-targets/, which holds genome files for closely related non-target species. These files follow the same naming rule: Genus_species_i.fasta.

At the top level of the input directory, there may also be a folder named non-targets/. This contains additional non-target genome files that will be used in all comparisons, regardless of species.

An example input directory structure:

```bash
input/
├── Yersinia_pseudotuberculosis/
│   ├── Yersinia_pseudotuberculosis_1.fasta
│   ├── Yersinia_pseudotuberculosis_2.fasta
│   └── non-targets/
│       ├── Yersinia_enterocolitica_1.fasta
│       └── Yersinia_intermedia_1.fasta
│
├── Escherichia_coli/
│   ├── Escherichia_coli_1.fasta
│   ├── Escherichia_coli_2.fasta
│   └── non-targets/
│       ├── Shigella_sonnei_1.fasta
│       └── Klebsiella_pneumoniae_1.fasta
│
└── non-targets/
    ├── Salmonella_enterica_1.fasta
    ├── Vibrio_cholerae_1.fasta
    └── Campylobacter_jejuni_1.fasta
```

## Configuration File: `config.yaml`

The configuration file defines all input paths, pipeline settings, and primer design parameters. Below is a breakdown of each section:

### `input_type`
Specifies the format of the input data.

- `raw` – Input consists of raw genome FASTA files (currently fully supported)
- `prokka` – Input is in the form of Prokka-annotated files (not yet supported)
- `panaroo` – Input comes from a Panaroo output folder (not yet supported)

### `input_paths`

- `raw_dir`: Path to the directory containing raw genome input folders (one per species)
- `prokka_dir`: Placeholder for Prokka input (set to `null`)
- `panaroo_dir`: Placeholder for Panaroo input (set to `null`)

### `output_dir`

- Directory where all output results, logs, and intermediate files will be stored.

### `max_cores`

- Maximum number of CPU cores to use for multi-threaded steps (e.g. alignment, BLAT, Panaroo).

### `aligner`

- Tool used for multiple sequence alignment.
- Options:
  - `mafft`
  - `clustal`
  - `prank`

### `snp_avg_prop_threshold`

- Minimum average proportion of informative SNPs a locus must have to be considered diagnostic and retained in downstream analysis.

### `primer3_config_file`

- Path to an external Primer3 config file (not fully implemented).
- If set, it overrides the `primer3.global_params` section.
- **Keep as `null` for now.**

### `primer3.global_params`

Parameters passed directly to Primer3 for primer design. These include:

- `PRIMER_NUM_RETURN`: Number of primer pairs Primer3 should return (default = 1)
- `PRIMER_MAX_NS_ACCEPTED`: Max number of ambiguous bases allowed in the template
- `PRIMER_LIBERAL_BASE`: Whether to use relaxed base recognition
- `PRIMER_MIN_SIZE` / `PRIMER_MAX_SIZE`: Minimum and maximum allowed primer lengths
- `PRIMER_MIN_TM` / `PRIMER_MAX_TM`: Minimum and maximum melting temperatures (°C)
- `PRIMER_PRODUCT_SIZE_RANGE`: Desired size range for PCR products (e.g. `"500-1800"`)

### `snp_primer_design`

Settings for SNP-aware primer design strategy:

- `snp_window_size`: Size of the window (in bp) used to scan for SNP-rich regions
- `snp_top_n`: Number of top-ranked loci to attempt primer design on
- `min_snps`: Minimum number of informative SNPs required within a window to proceed with primer design

### `validation` (not fully implemented yet)

Optional post-design validation settings using pBLAT:

- `perform`: Set to `yes` or `no` to enable/disable in silico validation
- `database`: Path to the reference `.2bit` database used for pBLAT validation
- `pblat_min_identity`: Minimum identity threshold (%) for considering alignments
- `match_median_filter_tolerance`: Allowed deviation in median match length between targets and non-targets

## Output

The pipeline generates an organized output directory containing results for each species. The top-level `output/` folder includes the following subdirectories:

### `prokka_output/`
Contains Prokka annotation results for each genome. Automatically generated from raw input if `input_type` is set to `raw`.

### `panaroo_output/`
Stores results from Panaroo pan-genome analysis, used internally for identifying shared and variable loci.

### `informative_loci/` **(Key output for users)**

This folder contains loci identified as having high diagnostic potential based on SNP profiles:

- **FASTA files** for each high-value locus
- `heatmap_informativeness.png`: visual overview of SNP informativeness across loci and species
- `snp_summary.csv`: tabular summary of all loci, including SNP positions, proportions, and informativeness scores

The file `snp_summary.csv` provides a detailed summary of all analyzed loci with respect to their SNP-based diagnostic potential. For each locus, it includes:

- The number and proportion of **informative SNPs** (per non-target species and averaged)
- The **positions** of informative SNPs
- The **median length** of aligned sequences in target and non-target genomes
- The **difference in median lengths**, which may reflect indel-based divergence
- A list of informative SNP **positions** and the **position range** that captures the majority of SNPs (±2 SD)

This file serves as the analytical foundation for selecting loci for downstream primer design.

#### Subfolders:

- **`consensus_sequences/`**: consensus FASTA sequences of informative loci, used for primer design
- **`snp_density_plots/`**: line plots showing SNP distributions across loci for visualization and interpretation

### `primers/` **(Key output for users)**

Contains designed primer pairs targeting the most informative SNP-rich regions:

- `primer_design_summary.csv`: detailed summary of primer sequences, SNP coverage, Tm, GC content, and amplicon size
- One FASTA file **per locus** containing the left and right primer sequences (`Locus_primers.fasta`)

All outputs are grouped by species to keep analyses modular and easily navigable.

