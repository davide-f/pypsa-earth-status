# PyPSA-Earth-Status: A package for validating any PyPSA network on Earth :D

## Development Status: **under development**

[![Test workflows](https://github.com/pypsa-meets-earth/pypsa-earth-status/actions/workflows/test.yaml/badge.svg)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/test.yaml)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth-status)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
<!-- [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth/main) -->
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing)


**PyPSA-Earth-Status: A package for validating any PyPSA network on Earth**

PyPSA-Earth-Status is a github repository designed to validate PyPSA networks against real-world data. It provides an automated procedure to compare PyPSA networks provided by the user with state-of-the-art statistics from available databases, such as IRENA, IEA, among others. PyPSA-Earth-Status aims to facilitate life of modelers in facilitating validation and quality of modelling results.

The model contains automated procedures for the validation of:
- **Installed capacities** of existing generators and storage units, namely the "p_nom" of generators and storage units
- **Optimal capacities**: the "p_nom_opt" of generators and storage units
- **Line capacities**: the "s_nom" of lines
- **Demand**

## Features
- Automated validation of PyPSA networks against real-world data
- Creation of reference statistics:
  - demand data from Our World in Data
  - installed capacities from IRENA and IEA
  - cross-border line capacities from [Global Transmission Database](https://zenodo.org/records/15527469)
- Comparison of network data with reference statistics
- Generation of tables and visualizations for easy interpretation of results

## Installation

To install PyPSA-Earth-Status, you can simply clone the repository:

```bash
    git clone https://github.com/pypsa-meets-earth/pypsa-earth-status --recurse-submodules
```

The python dependencies needed to execute the script are a proven subset of those of PyPSA-Earth and we recommend using the PyPSA-Earth environment. PyPSA-Earth provides OS-specific lock files for conda. For Linux, you can create the environment as follows:

```bash
    cd pypsa-earth-status
    conda env create -f workflows/pypsa-earth/envs/linux-64.lock.yaml
```

For other operating systems, you can use the same command with a different path in agreement to your OS and platform. You can find the available lock files in `workflows/pypsa-earth/envs/`, e.g. `win-64.lock.yaml` for Windows.

## Tutorial

To validate a first PyPSA network, you can run the following command in your terminal from the `pypsa-earth-status` folder:

```bash
    conda activate pypsa-earth
    snakemake -j 1 validate_data
```

This command will create the sample network scigrid_de from PyPSA and save it as `resources/example_DE.nc` with minimal changes and execute the validation procedure. The results will be saved in the folder `results/`, including tables and visualizations.

## Execute your first custom validation

If you want to validate your own PyPSA network, you can:

1. Open the file `config.yaml` in the `pypsa-earth-status` folder
2. Change the path of the network to your own network in the field `network_path` under `network_validation`
3. Adapt the list of countries you want to validate in the field `countries` under `network_validation` using 2-letter code naming convention; please keep at least two neighbouring countries in the list, e.g. ['DE', 'IT'] for the tutorial case for Germany.
4. Execute:
    ```bash
        snakemake -j 1 validate_data
    ```
5. Check the results in the `results/` folder
