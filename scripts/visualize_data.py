# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from helpers import (
    configure_logging,
    read_csv_nafix,
)
colors = [
    "#b80404",  
    "#0c6013",
    "#707070",
    "#ba91b1",
    "#6895dd",
    "#262626",
    "#235ebc",
    "#4adbc8",
    "#f9d002"
]

def plot_demand_comparison(demand_df, output_path):
    """
    Plot a side-by-side bar graph comparing electricity demand for each country (reference vs network).
    """

    # Set up the plot
    plt.figure(figsize=(10, 6))

    # Plot in reverse order: reference first, then network
    demand_df[['reference_demand', 'network_demand']].fillna(0).plot(kind='bar', stacked=False, color=colors[:2], zorder=3)

    plt.title("Electricity Demand Comparison (Reference vs Network)")
    plt.ylabel('Demand (TWh)')
    plt.xlabel('Country')
    plt.xticks(ticks=range(len(demand_df)), labels=demand_df['country'])

    # Remove plot box (spines)
    ax = plt.gca() 
    ax.spines['top'].set_visible(False)  
    ax.spines['right'].set_visible(False) 
    ax.spines['left'].set_visible(False)  

    # Add horizontal grid lines
    plt.grid(True, axis='y', zorder=0)

    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path)
    plt.close()


def plot_carrier_capacity_comparison(installed_capacity_df, optimal_capacity_df, output_path, carrier='coal', normalize='False'):
    """
    Plot a side-by-side bar graph comparing installed, optimal, and reference capacities for a given carrier (default: coal).
    """

    # Filter for the chosen carrier
    installed_capacity_df = installed_capacity_df[installed_capacity_df['carrier'] == carrier]
    optimal_capacity_df = optimal_capacity_df[optimal_capacity_df['carrier'] == carrier]
    reference_capacity_df = optimal_capacity_df[optimal_capacity_df['carrier'] == carrier]

    # Merge the dataframes
    capacity_df = pd.merge(installed_capacity_df[['country', 'network_capacity']], 
                           optimal_capacity_df[['country', 'network_capacity']], 
                           on=['country'], suffixes=('_network', '_optimal'))
    
    capacity_df = pd.merge(capacity_df, reference_capacity_df[['country', 'reference_capacity']], on='country', how='left')

    # Rename the columns to correct names
    capacity_df = capacity_df.rename(columns={'network_capacity_network': 'network_capacity', 'network_capacity_optimal': 'optimal_capacity'})

    if normalize == True:
        # Normalize data with respect to reference data (element-wise division)
        capacity_df['network_capacity'] = capacity_df['network_capacity'] / capacity_df['reference_capacity']
        capacity_df['optimal_capacity'] = capacity_df['optimal_capacity'] / capacity_df['reference_capacity']
        capacity_df['reference_capacity'] = capacity_df['reference_capacity'] / capacity_df['reference_capacity']
    
    # Set up the plot
    plt.figure(figsize=(10, 6))

    # Plot in reverse order: reference_capacity, network_capacity, optimal_capacity
    capacity_df[['reference_capacity', 'network_capacity', 'optimal_capacity']].fillna(0).plot(kind='bar', stacked=False, color=colors[:3], zorder=3)

    plt.title(f"{'Normalized ' if normalize else ''}Capacity Comparison for {carrier.capitalize()} per Country")
    plt.ylabel(f"Capacity {'(Ratio to Reference)' if normalize else '(MW)'}")
    plt.xlabel('Country')
    plt.xticks(ticks=range(len(capacity_df)), labels=capacity_df['country'])

    # Remove plot box (spines)
    ax = plt.gca() 
    ax.spines['top'].set_visible(False)  
    ax.spines['right'].set_visible(False) 
    ax.spines['left'].set_visible(False)  

    # Add horizontal grid lines
    plt.grid(True, axis='y', zorder=0)  # Enable grid lines along the y-axis
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path)
    plt.close()


def plot_stack_carrier_capacity_comparison(installed_capacity_df, optimal_capacity_df, output_path, stack_percent=False):
    """
    Plot a single stacked bar graph comparing network, optimal, and reference capacity mix for each country.
    Each country will have 3 bars side by side: one for reference capacity, one for network capacity, and one for optimal capacity.
    The bars will be stacked with different carriers' capacities.
    """

    # Merge the dataframes on country and carrier
    merged_df = pd.merge(installed_capacity_df[['carrier', 'country', 'network_capacity']], 
                         optimal_capacity_df[['carrier', 'country', 'network_capacity']], 
                         on=['carrier', 'country'], suffixes=('_network', '_optimal'))
    
    merged_df = pd.merge(merged_df, installed_capacity_df[['carrier', 'country', 'reference_capacity']], 
                         on=['carrier', 'country'], how='left')

    # Rename columns for clarity
    merged_df = merged_df.rename(columns={'network_capacity_network': 'network_capacity', 
                                          'network_capacity_optimal': 'optimal_capacity'})

    if stack_percent:
        # Normalize per country’s total capacity mix
        merged_df['network_capacity'] = merged_df.groupby('country')['network_capacity'].transform(lambda x: x*100 / x.sum())
        merged_df['optimal_capacity'] = merged_df.groupby('country')['optimal_capacity'].transform(lambda x: x*100 / x.sum())
        merged_df['reference_capacity'] = merged_df.groupby('country')['reference_capacity'].transform(lambda x: x*100 / x.sum())
    
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create a pivot table for each type of capacity
    network_pivot = merged_df.pivot_table(index='country', columns='carrier', values='network_capacity', aggfunc='sum')
    optimal_pivot = merged_df.pivot_table(index='country', columns='carrier', values='optimal_capacity', aggfunc='sum')
    reference_pivot = merged_df.pivot_table(index='country', columns='carrier', values='reference_capacity', aggfunc='sum')

    # Plot the stacked bars for each capacity type (reference, network, optimal)
    width = 0.25  # Width of the bars
    x = np.arange(len(network_pivot))  # X positions for the countries
    offsets = [-width, 0, width]  # Adjusted offsets

    # Stacked bar plots without duplicate legends
    reference_pivot.fillna(0).plot(kind='bar', stacked=True, ax=ax, width=width, position=1, color=colors, zorder=3)
    network_pivot.fillna(0).plot(kind='bar', stacked=True, ax=ax, width=width, position=0, color=colors, zorder=3)
    optimal_pivot.fillna(0).plot(kind='bar', stacked=True, ax=ax, width=width, position=-1, color=colors, zorder=3)

    # Get legend handles from just one of the plots to avoid duplicates
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:len(reference_pivot.columns)], labels[:len(reference_pivot.columns)], 
          title="Carriers", loc="center left", bbox_to_anchor=(1.02, 0.5), 
          ncol=1, frameon=False)

    # Create formatted x-axis labels (each country appears 3 times)
    xtick_labels = []
    xtick_positions = []
    for i, country in enumerate(network_pivot.index):
        xtick_labels.extend([
            f"{country}, REFR",
            f"{country}, BASE",
            f"{country}, OPTI"
        ])
        xtick_positions.extend([x[i] + offsets[0], x[i] + offsets[1], x[i] + offsets[2]])

    # Customize the plot
    ax.set_title("Capacity Mix Comparison per Country")
    ax.set_xlabel('Country')
    ax.set_ylabel(f"Capacity Mix {'(%)' if stack_percent else '(MW)'}")
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels)

    # Add grid lines and clean up the plot
    ax.grid(True, axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Adjust the layout for better spacing
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("visualize_data")

    configure_logging(snakemake)

    # Load comparison data
    demand_comparison = read_csv_nafix(snakemake.input["demand_comparison"])
    installed_capacity_comparison = read_csv_nafix(snakemake.input["installed_capacity_comparison"])
    optimal_capacity_comparison = read_csv_nafix(snakemake.input["optimal_capacity_comparison"])

    plot_demand_comparison(demand_comparison, snakemake.output['plot_demand'])

    # Compares capacities per country one carrier at a time
    # Select carrier value: ['solar' 'onwind' 'offwind-dc' 'coal' 'CCGT' 'ror' 'biomass' 'oil' 'geothermal']
    plot_carrier_capacity_comparison(installed_capacity_comparison, optimal_capacity_comparison, snakemake.output['plot_installed_capacity'], carrier='solar', normalize=True)

    # Compares network capacity mix per country with respect to reference with a stacked bargraph
    plot_stack_carrier_capacity_comparison(installed_capacity_comparison, optimal_capacity_comparison, snakemake.output['plot_capacity_mix'], stack_percent=True)