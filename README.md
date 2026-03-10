# TreeMS2

## Description

TreeMS2 is a Python-based tool for efficient comparison of large-scale tandem mass spectrometry (MS/MS) datasets. It
calculates pairwise distances between sets of MS/MS spectra—without requiring peptide or protein identification.

TreeMS2 is designed to support applications such as:

- Species identification
- Molecular phylogenetics
- Food and feed authentication
- Quality control in proteomics workflows

TreeMS2 requires the following input files:

- a set of MGF files: One or more Mascot Generic Format (MGF) files containing MS/MS spectra.
- A metadata file (CSV or TSV): A table that maps each MGF file to a user-defined label representing a set of spectra. A
  set of spectra in TreeMS2 can thus consist of spectra from multiple MGF files. The file paths listed in the metadata
  file must be relative to the location of the metadata file.

Each set of spectra is then compared against all others to compute a distance matrix.

An example for a valid TSV file might be:

```
file        group
file1.mgf   Gorilla
file2.mgf   Gorilla
file3.mgf   Chimpanzee
file4.mgf   Macaque
```

## Usage

### Installing the Conda environment

To get started with TreeMS2, it's recommended to use the provided Conda environment to ensure all dependencies are
correctly installed.

Clone the repository:

```bash
git clone https://github.com/bittremieuxlab/TreeMS2.git
cd TreeMS2
```

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate TreeMS2
```

### Configuration and Environment variables

To request an overview of the available parameters:

```bash
python3 -m TreeMS2 --help
```

Additionally, some environment variables specific to TreeMS2 can be set:

- `TREEMS2_NUM_CPUS`: Specifies how many CPU threads TreeMS2 should use.
- `TREEMS2_MEM_PER_CPU`: Memory in GB per core/cpu that TreeMS2 is allowed to use.

Some other environment variables, used by dependencies of TreeMS2:

Lance:

- `LANCE_IO_THREADS`
- `LANCE_CPU_THREADS`

NumPy, SciPy, BLAS, etc.:

- `MKL_NUM_THREADS`
- `NUMEXPR_NUM_THREADS`
- `BLIS_NUM_THREADS`
- `OPENBLAS_NUM_THREADS`
- `NUMBA_NUM_THREADS`

In general it is not necessary to set the environment variables explicitly, the application will determine appropriate
settings automatically.

### Running TreeMS2

Run the program by either providing a configuration file or by providing the required arguments:

```bash
python3 -m TreeMS2 -c config.ini
```

### Output

TreeMS2 generates several output files during execution, but the most important result is the pairwise distance matrix.
The distance matrix is always saved to `<work_directory>/results/distance_matrix.meg`. This file contains the computed
pairwise distances between all sets of spectra. It is formatted for compatibility with the MEGA (Molecular Evolutionary
Genetics Analysis) software (version 11), which can be used to construct
a rooted tree using the UPGMA algorithm. (Phylogeny → Construct/Test UPGMA Tree → Pairwise Distance → Lower Left Matrix)

## Considerations and known issues

- Additional testing should be done. Testing during development was mainly done on a
  CPU. Additional testing on a GPU should be done, especially for larger datasets

## Acknowledgements

TreeMS2 is inspired by and builds upon the distance metric introduced in compareMS2, developed by Marissen et al. (
2022).
