import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import string
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import ListedColormap


def plot_coloc(data, biobank=False):
    # Data preparation
    data['study'] = data['study'].str.lower()
    data['study_trimmed'] = data['study'].str.replace(r'[:punct:]|[:symbol:]', '', regex=True)
    data['study_trimmed'] = data['study_trimmed'].str.slice(0, 35)

    # Removing duplicates
    data = data.drop_duplicates(subset=['Source', 'leftVariant', 'gene_id'])

    # Selecting the top records
    data = data.sort_values(by='log2_H4_H3', ascending=False).groupby(['gene_id', 'leftVariant']).head(1)
    data = data[data['log2_H4_H3'] > 7]

    # Filtering for biobank studies if required
    if biobank:
        data = data[data['Source'].str.startswith('GCST')]

    # Plotting
    plt.figure(figsize=(10, 6))
    sns.barplot(data=data, y='study_trimmed', x='log2_H4_H3', hue='Source')
    plt.xlabel('')
    plt.ylabel('')
    plt.title('Colocalisation Plot')
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Adding labels
    for index, row in data.iterrows():
        plt.text(row['log2_H4_H3'], index, row['leftVariant'], color='black', ha="right")

    plt.show()

# Example usage
# plot_coloc(res, biobank=True)



def plot_l2g(data, disease_efo=None):
    # Exclude irrelevant trait categories
    exclude = ["phenotype", "measurement", "Uncategorised", "biological process"]
    data = data[~data['Trait_category'].isin(exclude)]

    # Select relevant columns and rename them
    df = data[['L2G_score', 'Distance', 'Interaction', 'mQTL', 'Pathogenicity', 'Gene_name',
               'Traits', 'EFO_ID', 'Trait_category', 'pval']].copy()
    df.columns = ['L2G_score', 'Distance', 'Interaction', 'mQTL', 'Pathogenicity', 'Gene_name',
                  'Traits', 'EFO_ID', 'Trait_category', 'pval']

    if disease_efo is not None:
        # Filter by disease EFO ID and select top-scoring gene for each trait
        df = df[df['EFO_ID'] == disease_efo]
        df = df.groupby('Gene_name').apply(lambda group: group[group['L2G_score'] == group['L2G_score'].max()])
        df = df.reset_index(drop=True)
        df_data = df.iloc[:, :6]

        # Generate radar plot with title based on the first trait
        plot = pn.ggplot(df_data, pn.aes(color='Gene_name')) + \
            pn.geom_radar(rescale=False, use_label=True, alpha=0.12, size=2) + \
            pn.labs(title=df.loc[0, 'Traits'])
    else:
        # Select top-scoring genes for each trait and plot in separate panels
        df = df.groupby('Traits').apply(lambda group: group.nlargest(3, 'L2G_score'))
        df = df.reset_index(drop=True)
        df_data = df.iloc[:, :7]

        # Generate radar plot with facetting by traits
        plot = pn.ggplot(df_data, pn.aes(color='Gene_name')) + \
            pn.geom_radar(rescale=False, use_label=True, size=2, alpha=0.12) + \
            pn.facet_wrap('~Traits') + \
            pn.theme(legend_position='right')

    return plot

# Example usage:
# plot = plot_l2g(data, disease_efo="EFO_0004339")
# print(plot)
def plot_manhattan(data):
    # Extract relevant columns and rename them
    gwas_results = data[['variant.position', 'variant.chromosome', 'pval', 'variant.id', 'bestLocus2Genes_score', 'bestLocus2Genes_gene_symbol']].copy()
    gwas_results = gwas_results.rename(columns={
    'variant.position': 'BP',
    'variant.chromosome': 'CHR',
    'pval': 'P',
    'variant.id': 'SNP'
})

# ... (the rest of your code remains the same)

    # Convert CHR column to integer
    gwas_results['CHR'] = gwas_results['CHR'].astype(int)

    # Calculate cumulative positions for each chromosome
    chr_len = gwas_results.groupby('CHR')['BP'].max().reset_index()
    chr_len['tot'] = chr_len.groupby('CHR')['BP'].cumsum() - chr_len['BP']
    gwas_results = gwas_results.merge(chr_len, on='CHR')
    gwas_results['BPcum'] = gwas_results['BP'] + gwas_results['tot']

    # Calculate center positions for each chromosome for axis labels
    df_axis = gwas_results.groupby('CHR').agg(center=('BPcum', lambda x: (max(x) + min(x)) / 2)).reset_index()

    # Select top genes based on bestLocus2Genes_score
    l2g_annot = gwas_results.groupby('CHR').apply(lambda x: x.nlargest(3, 'bestLocus2Genes_score')).reset_index(drop=True)

    # Set p-value cutoff
    p_cutoff = 1e-8

    # Make the plot
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.set(style="whitegrid")

    cmap = ListedColormap(['darkblue', 'darkgreen'])

    ax.scatter(gwas_results['BPcum'], -np.log10(gwas_results['P']), c=gwas_results['CHR'],
               cmap=cmap, alpha=0.85, s=50, linewidths=0.5, edgecolors='k')

    ax.set_xscale('symlog')
    ax.set_yscale('log')
    ax.set_xlabel('Chromosome', fontsize=14)
    ax.set_ylabel('-log10(P)', fontsize=14)

    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_major_locator(plt.FixedLocator(df_axis['center']))
    ax.set_xticklabels(df_axis['CHR'], rotation=45, fontsize=12)

    for i, row in l2g_annot.iterrows():
        if -np.log10(row['P']) > -np.log10(p_cutoff):
            ax.annotate(row['bestLocus2Genes_gene_symbol'], (row['BPcum'], -np.log10(row['P'])),
                        fontsize=10, textcoords="offset points", xytext=(0, 5), ha='center', va='bottom')

    plt.grid(True, which="both", ls="--", c='gray')
    plt.title("Manhattan Plot", fontsize=16)
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