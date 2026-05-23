"""otargenpy — Tidy Python interface to the Open Targets Platform GraphQL API."""

__author__ = "Amir Feizi"
__email__ = "afeizi@gmail.com"
__version__ = "2.0.1"

# --- Drug queries ---
from .drug import (
    adverse_events_query,
    indications_query,
    known_drugs_chembl_query,
    mechanisms_of_action_query,
    pharmacogenomics_chembl_query,
)

# --- Gene queries ---
from .gene import (
    comp_genomics_query,
    depmap_query,
    gene_ontology_query,
    genetic_constraint_query,
    hallmarks_query,
    interactions_query,
    known_drugs_gene_query,
    mouse_phenotypes_query,
    pathways_query,
    pharmacogenomics_gene_query,
    safety_query,
)

# --- Gene + disease evidence ---
from .disease import (
    chembl_query,
    clinvar_query,
    europe_pmc_query,
    gene_burden_query,
    genomics_england_query,
    gwas_credible_sets_query,
    orphanet_query,
    uniprot_literature_query,
)

# --- Variant queries ---
from .variant import (
    pharmacogenomics_variant_query,
    qtl_credible_sets_query,
    uniprot_variants_query,
    variant_effect_predictor_query,
    variant_effect_query,
)

# --- Genetics / GWAS ---
from .genetics import (
    gwas_colocalisation,
    gwas_credible_set,
    locus2gene_query,
    overlap_info_for_study,
    shared_trait_studies_query,
    variants_query,
)

# --- Plotting ---
from .plotting import (
    plot_adverse_events,
    plot_colocalisation,
    plot_indications,
    plot_interactions,
    plot_l2g,
)
