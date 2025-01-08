"""

This script collects clean statistics data and merge the datasets to create reference statistic to be used to validate energy systems.

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

        snakemake = mock_snakemake("build_reference_statistics")
    
    shutil.copy(snakemake.input["demand_owid"], snakemake.output["demand"])
    shutil.copy(snakemake.input["cap_irena"], snakemake.output["installed_capacity"])
