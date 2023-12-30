import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
