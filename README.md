<!-- badges -->
[![PyPI version](https://img.shields.io/pypi/v/otargenpy)](https://pypi.org/project/otargenpy/)
[![GitHub stars](https://img.shields.io/github/stars/amirfeizi/otargenpy?style=flat)](https://github.com/amirfeizi/otargenpy/stargazers)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/amirfeizi/otargenpy/blob/main/examples/quickstart.ipynb)

# otargenpy

**Tidy Python interface to the [Open Targets Platform](https://platform.opentargets.org) GraphQL API.**

Query genes, diseases, drugs, variants, and genetic evidence directly from Python and receive
analysis-ready pandas DataFrames — no manual JSON wrangling required.

Sister package of [otargen](https://github.com/amirfeizi/otargen) (R).

---

## Installation

**PyPI (stable)**

```bash
pip install otargenpy
```

**GitHub (recommended — latest features and bug fixes)**

```bash
pip install git+https://github.com/amirfeizi/otargenpy.git
```

---

## Quick start

Every function takes a single identifier and returns a pandas DataFrame.

### Drug queries (by ChEMBL ID)

```python
import otargenpy as ot

# Adverse events for imatinib
ae = ot.adverse_events_query("CHEMBL941")

# Mechanism of action
moa = ot.mechanisms_of_action_query("CHEMBL941")

# Drug indications with clinical stage
ind = ot.indications_query("CHEMBL941")
```

### Gene queries (by Ensembl ID)

```python
# Known drugs targeting TP53
drugs = ot.known_drugs_gene_query("ENSG00000141510")

# Cancer hallmarks for TP53
hall = ot.hallmarks_query("ENSG00000141510")

# Protein-protein interactions from IntAct
inter = ot.interactions_query("ENSG00000141510", source_database="intact", size=25)

# DepMap essentiality for EGFR
dep = ot.depmap_query("ENSG00000146648")

# Safety liabilities
safe = ot.safety_query("ENSG00000146648")
```

### Gene + disease evidence (by Ensembl ID + EFO ID)

```python
# ChEMBL evidence linking PARP1 to breast cancer
ev = ot.chembl_query("ENSG00000143799", "EFO_0000305")

# GWAS credible sets for PCSK9 and hyperlipidemia
gwas = ot.gwas_credible_sets_query("ENSG00000169174", "EFO_0004911")

# ClinVar evidence for BRCA1 and ovarian cancer
cv = ot.clinvar_query("ENSG00000012048", "EFO_0001075")
```

### Pharmacogenomics & variants

```python
# Pharmacogenomics for a drug
pgx = ot.pharmacogenomics_chembl_query("CHEMBL1016")

# Variant effect predictions
vep = ot.variant_effect_predictor_query("1_154453788_C_T")
```

### Genetics & colocalisation

```python
# Locus-to-gene predictions
l2g = ot.locus2gene_query("fa375739ca2a6b825ce5cc69d117e84b")

# GWAS colocalisation
coloc = ot.gwas_colocalisation("5a86bfd40d2ebecf6ce97bbe8a737512")
```

---

## Visualization

Built-in plotting functions turn query results into publication-ready figures.

```python
import otargenpy as ot

ae = ot.adverse_events_query("CHEMBL941")
ot.plot_adverse_events(ae)

inter = ot.interactions_query("ENSG00000141510", source_database="intact", size=25)
ot.plot_interactions(inter)

l2g = ot.locus2gene_query("fa375739ca2a6b825ce5cc69d117e84b")
ot.plot_l2g(l2g)

coloc = ot.gwas_colocalisation("5a86bfd40d2ebecf6ce97bbe8a737512")
ot.plot_colocalisation(coloc)

ind = ot.indications_query("CHEMBL941")
ot.plot_indications(ind)
```

| Function | Input | Plot type |
|---|---|---|
| `plot_adverse_events` | `adverse_events_query` | Lollipop chart with significance threshold |
| `plot_interactions` | `interactions_query` | Circular network graph |
| `plot_l2g` | `locus2gene_query` | Ranked bar chart of L2G scores |
| `plot_colocalisation` | `gwas_colocalisation` | H4 vs variant count scatter |
| `plot_indications` | `indications_query` | Faceted clinical stage chart |

---

## Available functions (40)

| Category | Functions |
|---|---|
| **Drug queries** | `adverse_events_query`, `indications_query`, `known_drugs_chembl_query`, `mechanisms_of_action_query`, `pharmacogenomics_chembl_query` |
| **Gene queries** | `comp_genomics_query`, `depmap_query`, `gene_ontology_query`, `genetic_constraint_query`, `hallmarks_query`, `interactions_query`, `known_drugs_gene_query`, `mouse_phenotypes_query`, `pathways_query`, `pharmacogenomics_gene_query`, `safety_query` |
| **Gene + disease** | `chembl_query`, `clinvar_query`, `europe_pmc_query`, `gene_burden_query`, `genomics_england_query`, `gwas_credible_sets_query`, `orphanet_query`, `uniprot_literature_query` |
| **Variant queries** | `pharmacogenomics_variant_query`, `qtl_credible_sets_query`, `uniprot_variants_query`, `variant_effect_predictor_query`, `variant_effect_query` |
| **Genetics / GWAS** | `gwas_colocalisation`, `gwas_credible_set`, `locus2gene_query`, `overlap_info_for_study`, `shared_trait_studies_query`, `variants_query` |
| **Visualization** | `plot_adverse_events`, `plot_colocalisation`, `plot_indications`, `plot_interactions`, `plot_l2g` |

---

## Citation

If you use `otargenpy` in your research, please cite:

> Feizi A, Ray D (2023). otargen: an R package for accessing
> and visualizing Open Targets Genetics data. *Bioinformatics*, 39(7).
> https://doi.org/10.1093/bioinformatics/btad441

---

## Contributing

Bug reports and feature requests: [GitHub Issues](https://github.com/amirfeizi/otargenpy/issues)

---

## License

MIT
