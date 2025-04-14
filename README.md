# Marker Loci Identification and Primer Design Tool

## Installation

This pipeline is distributed as a self-contained [Apptainer](https://apptainer.org/) (formerly known as Singularity) container, which includes all necessary dependencies and tools pre-installed.

To use the pipeline, you will need to have **Apptainer** installed on your system. Installation instructions are available on the [Apptainer documentation site](https://apptainer.org/docs/).

Once Apptainer is installed, no further setup is required. Simply download or clone the project and run the pipeline container using:

```bash
./primer_design_tool.sif
```

This single .sif file encapsulates the entire environment, ensuring consistent and reproducible results across different systems.

## Input

To run the pipeline, you need:

- A configuration file named `config.yaml`
- A properly structured input directory with raw genome assemblies

### Configuration File

The file `config.yaml` must be placed in the **same directory as the container** (`primer_design_tool.sif`). It defines all global parameters, paths, and settings required for the analysis.

---

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
