"""
This script generates the configuration file for the simulation.

Starting by a base configuration file, it will create changes according to the provided yaml file.
"""

import pandas as pd
from helpers import (
    three_2_two_digits_country,
    read_csv_nafix,
    to_csv_nafix,
    configure_logging,
    country_name_2_two_digits,
)
import os

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("config_generator")
    
    configure_logging(snakemake)