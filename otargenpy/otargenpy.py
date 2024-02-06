"""Main module."""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
#from fuzzywuzzy import fuzz
import matplotlib.pyplot as plt
import seaborn as sns

def plot_coloc(data, biobank=False, plot_type="talk"):
    """
    Plot colocalization analysis based on the provided data.

    Parameters
    ----------
    data : DataFrame
    A pandas DataFrame containing colocalization analysis data.
    biobank : bool, optional
        If True, filters the data to include only studies starting with 'GCST',
        by default False.
    plot_type : str, optional
    Seaborn plot style to use. Options include "paper", "notebook", "talk",
    "poster". Default is "talk."

    Notes
    -----
    The function requires 'matplotlib' and 'seaborn' libraries for plotting.
    """
    if biobank:
        data = data[data["studyId"].str.startswith("GCST")]

    data["traitReported"] = data["traitReported"].apply(
        lambda x: x[:35] if len(x) > 35 else x
    )
    data.drop_duplicates(
        subset=["studyId", "variant_id", "gene_id"], inplace=True
    )

    data_sorted = (
        data.groupby(["gene_id", "variant_id"])
        .apply(lambda x: x.nlargest(3, "log2h4h3"))
        .reset_index(drop=True)
    )
    data_sorted = data_sorted[data_sorted["log2h4h3"] > 7]

    color_palette = sns.color_palette(
        "husl", len(data_sorted["studyId"].unique())
    )
    color_mapping = {
        study: color
        for study, color in zip(data_sorted["studyId"].unique(), color_palette)
    }
    data_sorted["color"] = data_sorted["studyId"].map(color_mapping)

    plt.figure(figsize=(10, 6))
    bars = plt.barh(
        data_sorted["traitReported"],
        data_sorted["log2h4h3"],
        color=data_sorted["color"],
        alpha=0.7,
    )

    for bar, rsId in zip(bars, data_sorted["rsId"]):
        plt.text(
            bar.get_width() + 0.1,
            bar.get_y() + bar.get_height() / 2,
            rsId,
            ha="left",
            va="center",
        )

    plt.xlabel("log2(h4/h3)")
    plt.ylabel("Trait Reported")
    plt.title("Colocalization Analysis")
    plt.grid(axis="x", linestyle="--", alpha=0.6)
    plt.show()


def plot_l2g(data, efo_id=None, plot_type="talk"):
    """
    Plot Locus-to-Gene (L2G) modeling results based on the provided data.

    Parameters
    ----------
    data : DataFrame
    A pandas DataFrame containing L2G modeling data.
    efo_id : str or list, optional
    EFO (Experimental Factor Ontology) ID(s) to filter the data,
    default is None.
    plot_type : str, optional
        Seaborn plot style to use. Options include "talk," "poster," etc.
        Default is "talk."

    Returns
    -------
    Seaborn axis-level plot
        A seaborn plot object representing the L2G data.

    Notes
    -----
    The function modifies column names and performs data transformation for
    plotting.
    """
    # Exclude irrelevant trait categories
    exclude = [
        "phenotype",
        "measurement",
        "Uncategorised",
        "biological process",
    ]
    data = data[~data["traitCategory"].isin(exclude)]
    data_f = data[
        [
            "symbol",
            "traitReported",
            "traitEfos",
            "yProbaModel",
            "yProbaDistance",
            "yProbaInteraction",
            "yProbaMolecularQTL",
            "yProbaPathogenicity",
            "pval",
        ]
    ]

    data_f.columns = [
        "symbol",
        "traitReported",
        "traitEfos",
        "L2G",
        "Distance",
        "Interaction",
        "MolecularQTL",
        "Pathogenicity",
        "pval",
    ]
    melted_df = pd.melt(
        data_f,
        id_vars=["symbol", "traitReported", "traitEfos", "pval"],
        value_vars=[
            "L2G",
            "Distance",
            "Interaction",
            "MolecularQTL",
            "Pathogenicity",
        ],
        var_name="L2GMODEL",
        value_name="MODEL_SCORE",
    )
    # cleaning the trait terms
    split_chars = r"[(,:|;]"

    melted_df["traitReported_trimmed"] = (
        melted_df["traitReported"].str.split(split_chars).str[0]
    )
    melted_df["traitReported_trimmed"] = (
        melted_df["traitReported_trimmed"].str.lower().str.strip()
    )

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
        current_disease = row["traitReported_trimmed"]

        # Find the best match for the current disease in the remaining rows
        remaining_rows = collapsed_df.loc[index + 1 :, "traitReported_trimmed"]
        best_match = find_best_match(current_disease, remaining_rows)

        # Replace the current disease with the best match
        collapsed_df.at[index, "traitReported_trimmed"] = best_match

        # Replace the longer version with the shorter version
        if best_match and len(best_match) < len(current_disease):
            collapsed_df["traitReported_trimmed"] = collapsed_df[
                "traitReported_trimmed"
            ].replace(current_disease, best_match)

    if efo_id is None:

        collapsed_df["rank"] = collapsed_df.groupby(
            ["L2GMODEL", "traitReported_trimmed"]
        )["MODEL_SCORE"].rank(ascending=False)

        collapsed_df = collapsed_df[
            collapsed_df["L2GMODEL"].isin(
                ["L2G", "Interaction", "MolecularQTL"]
            )
        ]
        sns.catplot(
            x="MODEL_SCORE",
            y="traitReported_trimmed",
            data=collapsed_df,
            kind="box",
            col="symbol",
            hue="L2GMODEL",
            col_wrap=2,
            palette="Set3",
            legend=True,
        )

    else:
        collapsed_df = collapsed_df[collapsed_df["traitEfos"].isin(efo_id)]
        sns.catplot(
            x="MODEL_SCORE",
            y="L2GMODEL",
            data=collapsed_df,
            kind="box",
            col="symbol",
            hue="L2GMODEL",
            col_wrap=2,
            palette="Set3",
            legend=True,
        )

    plt.show()


def plot_manhattan(data, plot_type="talk"):
    """
    Generate a Manhattan plot from GWAS results.

    Parameters
    ----------
    data : DataFrame
    A pandas DataFrame containing GWAS results with columns for variant
    position, chromosome, p-value, and rsID.
    plot_type : str, optional
    Seaborn plot style to use. Options include "paper", "notebook",
    "talk", "poster". Default is "talk."

    Notes
    -----
    This function requires 'pandas', 'numpy', 'matplotlib',
    and 'seaborn' libraries. It creates a Manhattan plot to visualize
    GWAS results, highlighting significant SNPs.
    """
    # Extract relevant columns and rename them
    gwas_results = data[
        ["variant.position", "variant.chromosome", "pval", "variant.rsId"]
    ].copy()
    gwas_results = gwas_results.drop_duplicates()

    # Simulate DataFrame
    df = pd.DataFrame(
        {
            "rsid": data["variant.rsId"],
            "chrom": data["variant.chromosome"],
            "pos": data["variant.position"],
            "pval": data["pval"],
        }
    ).drop_duplicates()
    # Calculate -log10 of p-values
    df["-logp"] = -np.log10(df.pval)
    df.sort_values(["chrom", "pos"], inplace=True)
    df.reset_index(inplace=True, drop=True)
    df["i"] = df.index
    # Calculate significance cutoff for each chromosome
    df["p_cutoff"] = df.groupby("chrom")["pval"].transform(
        lambda x: x.quantile(0.05)
    )  # Adjust quantile as needed

    # Set Seaborn style
    sns.set(style="whitegrid", context=plot_type, palette="bright")

    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    plot = sns.relplot(
        data=df,
        x="i",
        y="-logp",
        aspect=3.7,
        hue="chrom",
        palette="bright",
        legend=None,
    )
    chrom_df = df.groupby("chrom")["i"].median()
    plot.ax.set_xlabel("chrom")
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.fig.suptitle("Manhattan plot")

    # Add SNP labels for outliers
    outliers = df[df["pval"] < df["p_cutoff"]]

    for index, row in outliers.iterrows():
        plot.ax.annotate(
            row["rsid"],
            (row["i"], row["-logp"]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
        )

    plt.show()


def plot_phewas(data, disease=True, source=["GCST", "FINNGEN", "NEALE", "SAIGE"], plot_type="talk"
):
    """
    Generate a PheWAS (Phenome-Wide Association Study) plot from the
    given dataset.

    Parameters
    ----------
    data : DataFrame
        A pandas DataFrame containing PheWAS results with columns for study
        traits, beta values, p-values, and sources.
    disease : bool, optional
        If True, filters the data to include only disease-related traits;
        otherwise, includes other traits. Default is True.
    source : list of str, optional
        List of sources to include in the plot. Default is
        ["GCST", "FINNGEN", "NEALE", "SAIGE"].
    plot_type : str, optional
        Seaborn plot style to use. Options include "paper", "notebook", "talk",
        "poster"
        Default is "talk."

    Notes
    -----
    The function requires 'pandas', 'numpy', 'matplotlib', and
    'seaborn' libraries.
    It creates a scatter plot to visualize PheWAS results.
    """
    # Prepare the data
    dt0 = data.copy()
    dt0["traitCategory"] = dt0["study.traitCategory"].str.lower()
    dt0["traitReported_trimmed"] = dt0["study.traitReported"].str.replace(
        r"[:punct:]|[:symbol:]", "", regex=True
    )
    dt0["traitReported_trimmed"] = dt0["traitReported_trimmed"].str.slice(
        stop=35
    )

    # Filter and categorize the data
    if disease:
        dt2 = dt0[
            (dt0["study.source"].isin(source))
            & (
                ~dt0["traitCategory"].isin(
                    [
                        "measurement",
                        "phenotype",
                        "biological process",
                        "uncategorized",
                    ]
                )
            )
        ].copy()
        dt2["beta_direction"] = dt2["beta"].apply(
            lambda x: (
                "positive"
                if x > 0
                else ("negative" if x < 0 else "no information")
            )
        )
    else:
        dt2 = dt0[
            (dt0["study.source"].isin(source))
            & (
                dt0["traitCategory"].isin(
                    [
                        "measurement",
                        "phenotype",
                        "biological process",
                        "uncategorized",
                    ]
                )
            )
        ].copy()
        dt2["beta_direction"] = dt2["beta"].apply(
            lambda x: (
                "positive"
                if x > 0
                else ("negative" if x < 0 else "no information")
            )
        )

    # Create a dictionary to map sources to colors
    source_colors = {
        "NEALE": "red",
        "SAIGE": "blue",
        "FINNGEN": "green",
        "GCST": "purple",
    }
    dt2["source_color"] = dt2["study.source"].map(source_colors)

    # Define custom markers
    custom_markers = {
        "positive": "^",  # Triangle up for positive
        "negative": "v",  # Triangle down for negative
        "no information": "o",  # Circle for no information
    }

    # Set Seaborn style for a modern look
    sns.set(style="whitegrid", context=plot_type, palette="muted")

    # Create the plot
    plt.figure(figsize=(14, 8))
    sns.scatterplot(
        x="traitCategory",
        y=-np.log10(dt2["pval"]),
        hue="beta_direction",
        style="beta_direction",
        data=dt2,
        markers=custom_markers,
        s=100,  # Size of the dots
        alpha=0.7,  # Transparency of the dots
        edgecolor="w",  # White edges for the markers
        linewidth=0.5,  # Edge line width
    )

    # Calculate the 75th percentile of -log10(pval) for each category
    percentile_75 = dt2.groupby("traitCategory")["pval"].apply(
        lambda x: np.percentile(-np.log10(x), 75)
    )

    # Add text to points above the 75th percentile in each category
    for index, row in dt2.iterrows():
        category_75th_percentile = percentile_75[row["traitCategory"]]
        if -np.log10(row["pval"]) > category_75th_percentile:
            # Trim the text to 15 characters
            label = row["study.traitReported"][:30]
            plt.text(
                row["traitCategory"],
                -np.log10(row["pval"]),
                label,
                ha="center",
                va="bottom",
                fontsize=8,
            )

    # p-value threshold (e.g., p=0.05, -log10(p) = 1.3)
    pvalue_threshold = 0.05
    threshold_line = -np.log10(pvalue_threshold)
    plt.axhline(y=threshold_line, color="grey", linestyle="--", linewidth=1.5)

    # Improve the aesthetics
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Trait Category")
    plt.ylabel("-log10(pval)")
    plt.title("PheWAS Plot")
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.legend(
        title="Beta Direction", bbox_to_anchor=(1.05, 1), loc="upper left"
    )
    plt.tight_layout()

    plt.show()

