# -*- coding: utf-8 -*-
"""

This script collects clean statistics data and merges the datasets to create reference statistics to be used to validate energy systems.

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

def filter_data_by_config(df, column, valid_values):
    """
    Filters the dataframe based on the provided column and valid values from config.yaml.
    """
    return df[df[column].isin(valid_values)]

def process_reference_statistics(inputs, outputs, config):
    """
    Processes demand and installed capacity data based on the specified year and countries.
    """
    year = config["network_validation"]["year"][0]  # Extracting the year from config
    countries = config["network_validation"]["countries"]  # Extracting the list of countries
    
    # Process demand data
    df_demand = read_csv_nafix(inputs["demand_owid"])
    df_demand = filter_data_by_config(df_demand, "iso_code_2", countries)
    df_demand = df_demand[df_demand["year"] == year]
    df_demand = df_demand[["iso_code_2","electricity_demand"]].rename(columns={"electricity_demand": "demand","iso_code_2": "country"}).set_index("country")
    to_csv_nafix(df_demand, outputs["demand"])
    
    # Process installed capacity data
    df_capacity = read_csv_nafix(inputs["cap_irena"])
    df_capacity = filter_data_by_config(df_capacity, "alpha2", countries)
    df_capacity = df_capacity[df_capacity["Year"] == year]
    df_capacity = df_capacity.rename(columns={"Technology": "carrier","alpha2": "country"})
    df_capacity = df_capacity[["country","carrier","p_nom"]].set_index("country")
    to_csv_nafix(df_capacity, outputs["installed_capacity"])

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake
        snakemake = mock_snakemake("build_reference_statistics")
    
    configure_logging(snakemake)
    
    process_reference_statistics(snakemake.input, snakemake.output, snakemake.config)
