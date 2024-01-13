import matplotlib.pyplot as plt
import plotly.graph_objects as go

import seaborn as sns
import pandas as pd
import numpy as np
import string
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import ListedColormap


def plot_coloc(data, biobank=False):
    # Filter data if biobank is True
    if biobank:
        data = data[data['studyId'].str.startswith('GCST')]

    # Trim 'traitReported' if they have more than 35 characters
    data.loc[:, 'traitReported'] = data['traitReported'].apply(lambda x: x[:35] if len(x) > 35 else x)

    # Remove duplicated rows for 'studyId', 'variant_id', 'gene_id' combination
    data = data.drop_duplicates(subset=['studyId', 'variant_id', 'gene_id'])

    # Group by 'gene_id' and 'variant_id', and sort by 'log2h4h3'
    data_sorted = data.groupby(['gene_id', 'variant_id']).apply(lambda x: x.nlargest(3, 'log2h4h3')).reset_index(drop=True)

    # Select rows with 'log2h4h3' > 7
    data_sorted = data_sorted[data_sorted['log2h4h3'] > 7]

    # Create a categorical color palette based on unique 'studyId' values
    color_palette = sns.color_palette("husl", len(data_sorted['studyId'].unique()))
    color_mapping = {study: color for study, color in zip(data_sorted['studyId'].unique(), color_palette)}
    data_sorted['color'] = data_sorted['studyId'].map(color_mapping)

    # Create the barplot
    plt.figure(figsize=(10, 6))
    bars = plt.barh(data_sorted['traitReported'], data_sorted['log2h4h3'], color=data_sorted['color'], alpha=0.7)

    # Add 'rsId' as text on each corresponding bar
    for bar, rsId in zip(bars, data_sorted['rsId']):
        plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2, rsId, ha='left', va='center')

    plt.xlabel('log2(h4/h3)')
    plt.ylabel('Trait Reported')
    plt.title('Colocalization Analysis')
    plt.grid(axis='x', linestyle='--', alpha=0.6)

    plt.show()

# Example usage
# plot_coloc(res, biobank=True)



def plot_l2g(data, disease_efo=None):
    # Exclude irrelevant trait categories
    exclude = ["phenotype", "measurement", "Uncategorised", "biological process"]
    data = data[~data['traitCategory'].isin(exclude)]

    #Select relevant columns and rename them
    df = data[['L2G', 'Distance', 'Interaction', 'mQTL', 'Pathogenicity', 'gene_symbol',
              'traitReported', 'traitEfos',"traitReported" ,'pval']].copy()


    if disease_efo is not None:
        # Filter by disease EFO ID and select top-scoring gene for each trait
        df = df[df['traitEfos'] == disease_efo]
        df = df.groupby('gene_symbol').apply(lambda group: group[group['L2G'] == group['L2G'].max()])
        df = df.reset_index(drop=True)
        df_data = df.iloc[:, :6] # number of variable
        
        categories=list(df_data)[6:]
        N = len(categories)
        
        # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
        angles = [n / float(N) * 2 * pi for n in range(N)]
        angles += angles[:1]
        
        # Initialise the spider plot
        ax = plt.subplot(111, polar=True)
        
        # If you want the first axis to be on top:
        ax.set_theta_offset(pi / 2)
        ax.set_theta_direction(-1)
        
        # Draw one axe per variable + add labels
        plt.xticks(angles[:-1], categories)
        
        # Draw ylabels
        ax.set_rlabel_position(0)
        plt.yticks([10,20,30], ["10","20","30"], color="grey", size=7)
        plt.ylim(0,40)
        

        # ------- PART 2: Add plots
        
        # Plot each individual = each line of the data
        # I don't make a loop, because plotting more than 3 groups makes the chart unreadable
        
        # Ind1
        values=df.loc[0].drop('group').values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=1, linestyle='solid', label="group A")
        ax.fill(angles, values, 'b', alpha=0.1)
        
        # Ind2
        values=df.loc[1].drop('group').values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=1, linestyle='solid', label="group B")
        ax.fill(angles, values, 'r', alpha=0.1)
        
        # Add legend
        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

        # Show the graph
        plt.show()
                

      
# Example usage:
# plot = plot_l2g(data, disease_efo="EFO_0004339")
# print(plot)

def plot_manhattan(data):
    # Extract relevant columns and rename them
    gwas_results = data[['variant.position', 'variant.chromosome', 'pval', 'variant.rsId']].copy()
    gwas_results = gwas_results.drop_duplicates()

    # Simulate DataFrame
    df = pd.DataFrame({
        'rsid': gwas_results["variant.rsId"],
        'chrom': gwas_results['variant.chromosome'],
        'pos': gwas_results['variant.position'],
        'pval': gwas_results['pval']
    })
    df['-logp'] = -np.log10(df.pval)
    df = df.sort_values(['chrom', 'pos'])
    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    # Calculate significance cutoff for each chromosome
    df['p_cutoff'] = df.groupby('chrom')['pval'].transform(lambda x: x.quantile(0.05))  # Adjust quantile as needed

    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    plot = sns.relplot(data=df, x='i', y='-logp', aspect=3.7,
            hue='chrom', palette='bright', legend=None)
    chrom_df = df.groupby('chrom')['i'].median()
    plot.ax.set_xlabel('chrom')
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.fig.suptitle('Manhattan plot')

    # Add SNP labels for outliers
    outliers = df[df['pval'] < df['p_cutoff']]
    
    for index, row in outliers.iterrows():
        plot.ax.annotate(row['rsid'], (row['i'], row['-logp']),
                textcoords="offset points", xytext=(0, 10), ha='center')

    plt.show()
# Example usage:
# plot_manhattan(data)
def plot_phewas(data, disease=True, source=["GCST", "FINNGEN", "NEALE", "SAIGE"]):
    # Prepare the data
    dt0 = data.copy()
    dt0['traitCategory'] = dt0['study.traitCategory'].str.lower()
    dt0['traitReported_trimmed'] = dt0['study.traitReported'].str.replace(r'[:punct:]|[:symbol:]', '', regex=True)
    dt0['traitReported_trimmed'] = dt0['traitReported_trimmed'].str.slice(stop=35)

    # Filter and categorize the data
    if disease:
        dt2 = dt0[
            (dt0['study.source'].isin(source)) &
            (~dt0['traitCategory'].isin(["measurement", "phenotype", "biological process", "uncategorized"]))
        ].copy()
        dt2['beta_shape'] = np.where(dt2['beta'] > 0, "positive", "negative")
    else:
        dt2 = dt0[
            (dt0['study.source'].isin(source)) &
            (dt0['traitCategory'].isin(["measurement", "phenotype", "biological process", "uncategorized"]))
        ].copy()
        dt2['beta_shape'] = np.where(dt2['beta'] > 0, "positive", "negative")

    # Create a dictionary to map sources to colors
    source_colors = {"NEALE": "red", "SAIGE": "blue", "FINNGEN": "green", "GCST": "purple"}

    # Map source names to colors
    dt2['source_color'] = dt2['study.source'].map(source_colors)

    # Create the plot
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(
        x=dt2['traitCategory'],
        y=-np.log10(dt2['pval']),
        c=dt2['source_color'],  # Use the mapped source colors
        s=50,
        alpha=0.5
    )

    handles = [mpatches.Patch(color=color, label=source) for source, color in source_colors.items()]

    plt.legend(handles=handles, title="Data source", loc="upper right", bbox_to_anchor=(1.0, 1.0))
    plt.axhline(y=5, color="grey", linestyle='--')
    plt.xlabel("")
    plt.ylabel("-log10(pval)")
    plt.xticks(rotation=45, ha="right")
    plt.title("Phewas Plot")
    plt.tight_layout()
    plt.show()


# Example usage:
# plot_phewas(data)