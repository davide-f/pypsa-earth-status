"""

This script reads a PyPSA network and builds reference statistics to be used for comparison.

"""
import pandas as pd
import shutil
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

        snakemake = mock_snakemake("build_network_statistics")