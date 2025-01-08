# -*- coding: utf-8 -*-
"""

This script collects clean statistics data and merge the datasets to create reference statistic to be used to validate energy systems.

"""
import os
import shutil

import pandas as pd
from helpers import (
    configure_logging,
    country_name_2_two_digits,
    read_csv_nafix,
    three_2_two_digits_country,
    to_csv_nafix,
)

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("build_reference_statistics")

    shutil.copy(snakemake.input["demand_owid"], snakemake.output["demand"])
    shutil.copy(snakemake.input["cap_irena"], snakemake.output["installed_capacity"])
