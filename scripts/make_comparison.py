# -*- coding: utf-8 -*-
"""
This script compares the reference and network statistics by searching unique
region and carrier combinations for capacities and unique countries for demand.
"""
import os
import shutil
import pandas as pd
from helpers import (
    configure_logging,
    country_name_2_two_digits,
    read_csv_nafix,
    to_csv_nafix,
)

def compare_capacity_statistics(reference_df, network_df):
    """
    Compare the installed and optimal capacities between the reference and network data.
    Calculate the difference and ratio for both capacities with respect to reference capacities.
    """
    # Prepare empty list to store comparison results
    comparison_results = []
    
    # Find unique region and carrier combinations in the reference data
    unique_combinations = reference_df[['region', 'carrier']].drop_duplicates()

    for index, row in unique_combinations.iterrows():
        region = row['region']
        carrier = row['carrier']
        
        # Find the matching rows for each combination
        reference_row = reference_df[(reference_df['region'] == region) & 
                                     (reference_df['carrier'] == carrier)]
        network_row = network_df[(network_df['region'] == region) & 
                                 (network_df['carrier'] == carrier)]

        if not reference_row.empty and not network_row.empty:
            comparison_results.append({
                    'region': region,
                    'carrier': carrier,
                    'network_capacity': network_row['p_nom'].values[0],
                    'reference_capacity': reference_row['p_nom'].values[0]
                })
    
    # Convert the comparison results into a dataframe
    comparison_df = pd.DataFrame(comparison_results)
    comparison_df = comparison_df.set_index('region')
    
    return comparison_df

def compare_demand_statistics(reference_df, network_df):
    """
    Compare the demand between the reference and network data.
    Calculate the difference and ratio with respect to reference demand.
    """
    # Prepare empty list to store comparison results
    comparison_results = []
    
    # Find unique countries in the reference data
    unique_countries = reference_df['region'].drop_duplicates()

    for region in unique_countries:
        # Find the matching rows for each region
        reference_row = reference_df[reference_df['region'] == region]
        network_row = network_df[network_df['region'] == region]

        if not reference_row.empty and not network_row.empty:
            comparison_results.append({
                'region': region,
                'network_demand': network_row['demand'].values[0],
                'reference_demand': reference_row['demand'].values[0]
            })

    # Convert the comparison results into a dataframe
    comparison_df = pd.DataFrame(comparison_results)
    comparison_df = comparison_df.set_index('region')
    
    return comparison_df

def make_comparison(inputs, outputs):
    """
    This function loads the reference and network statistics, performs comparisons,
    and saves the results to CSV files.
    """
    # Load reference statistics
    df_reference_installed_capacity = read_csv_nafix(inputs["installed_capacity_reference"])
    df_reference_optimal_capacity = read_csv_nafix(inputs["installed_capacity_reference"])  # Assuming the same reference file
    df_reference_demand = read_csv_nafix(inputs["demand_reference"])

    # Load network statistics
    df_network_installed_capacity = read_csv_nafix(inputs["installed_capacity_network"])
    df_network_optimal_capacity = read_csv_nafix(inputs["optimal_capacity_network"])
    df_network_demand = read_csv_nafix(inputs["demand_network"])

    # Perform comparison for installed capacity
    installed_capacity_comparison = compare_capacity_statistics(
        df_reference_installed_capacity, df_network_installed_capacity
    )

    # Perform comparison for optimal capacity
    optimal_capacity_comparison = compare_capacity_statistics(
        df_reference_optimal_capacity, df_network_optimal_capacity
    )
    
    # Perform comparison for demand
    demand_comparison = compare_demand_statistics(
        df_reference_demand, df_network_demand
    )

    # Save the results to the output files
    to_csv_nafix(installed_capacity_comparison, outputs["installed_capacity_comparison"])
    to_csv_nafix(optimal_capacity_comparison, outputs["optimal_capacity_comparison"])
    to_csv_nafix(demand_comparison, outputs["demand_comparison"])

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("make_comparison")

    configure_logging(snakemake)

    make_comparison(snakemake.input, snakemake.output)