import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import string
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import ListedColormap
from fuzzywuzzy import fuzz


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

def plot_l2g(data, efo_id=None):
    # Exclude irrelevant trait categories
    
    exclude = ["phenotype", "measurement", "Uncategorised", "biological process"]
    data = data[~data['traitCategory'].isin(exclude)]
    data_f = data[['symbol','traitReported',"traitEfos",'yProbaModel', 'yProbaDistance', 'yProbaInteraction',
       'yProbaMolecularQTL', 'yProbaPathogenicity', 'pval']]
    
    data_f.columns = ['symbol','traitReported',"traitEfos",'L2G', 'Distance', 'Interaction',
       'MolecularQTL', 'Pathogenicity', 'pval']
    melted_df = pd.melt(data_f, id_vars=['symbol', 'traitReported', 'traitEfos', 'pval'], 
                    value_vars=['L2G', 'Distance', 'Interaction', 'MolecularQTL', 'Pathogenicity'], 
                    var_name='L2GMODEL', value_name='MODEL_SCORE')
    # cleaning the trait terms
    split_chars = r"[(,:|;]"

    melted_df["traitReported_trimmed"] = melted_df["traitReported"].str.split(split_chars).str[0]
    melted_df["traitReported_trimmed"] = melted_df["traitReported_trimmed"].str.lower().str.strip()



# Function to find the best match for a given string in a list of strings
    def find_best_match(string, string_list):
        best_match = None
        best_ratio = -1
        for s in string_list:
            ratio = fuzz.ratio(string, s)
            if ratio > best_ratio:
                best_ratio = ratio
                best_match = s
        return best_match

# Create a copy of the dataframe
    collapsed_df = melted_df.copy()

# Iterate over each row in the dataframe
    for index, row in collapsed_df.iterrows():
    # Get the current disease name
        current_disease = row['traitReported_trimmed']
    
    # Find the best match for the current disease in the remaining rows
        remaining_rows = collapsed_df.loc[index+1:, 'traitReported_trimmed']
        best_match = find_best_match(current_disease, remaining_rows)
    
    # Replace the current disease with the best match
        collapsed_df.at[index, 'traitReported_trimmed'] = best_match
    
    # Replace the longer version with the shorter version
        if best_match and len(best_match) < len(current_disease):
            collapsed_df['traitReported_trimmed'] = collapsed_df['traitReported_trimmed'].replace(current_disease, best_match)



    
    
    
    if efo_id is None:

        
        df["rank"] = collapsed_df.groupby(["L2GMODEL", "traitReported_trimmed"])["MODEL_SCORE"].rank(ascending=False)
        df_ft = df[df["rank"] == 1]
        
        collapsed_df = collapsed_df[collapsed_df["L2GMODEL"].isin(["L2G","Interaction","MolecularQTL"])]
        p = sns.catplot(x="MODEL_SCORE", y="traitReported_trimmed", data=collapsed_df,
                                kind="box", col="symbol", hue="L2GMODEL", col_wrap=2,
                                palette="Set3", legend=True)
        
    else:
        collapsed_df = collapsed_df[collapsed_df["traitEfos"].isin(efo_id)]
        p = sns.catplot(x="MODEL_SCORE", y="L2GMODEL", data=collapsed_df,
                                kind="box", col="symbol", hue="L2GMODEL", col_wrap=2,
                                palette="Set3", legend=True)
        
    return p

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