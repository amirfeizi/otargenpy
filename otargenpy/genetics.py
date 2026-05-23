"""Genetics / GWAS queries (by study locus, study ID, disease IDs)."""

import pandas as pd
from ._client import query_api, flatten_to_df, PLATFORM_API, GENETICS_API

# ---------------------------------------------------------------------------
# gwas_colocalisation
# ---------------------------------------------------------------------------

_GWAS_COLOC_QUERY = """
query GWASColocQuery($studyLocusId: String!, $size: Int!, $index: Int!) {
  credibleSet(studyLocusId: $studyLocusId) {
    colocalisation(studyTypes: [gwas], page: { size: $size, index: $index }) {
      count
      rows {
        otherStudyLocus {
          studyLocusId
          study { id projectId traitFromSource publicationFirstAuthor }
          variant { id chromosome position referenceAllele alternateAllele }
          pValueMantissa pValueExponent
        }
        numberColocalisingVariants colocalisationMethod h3 h4 clpp betaRatioSignAverage
      }
    }
  }
}
"""


def gwas_colocalisation(study_locus_id: str, size: int = 500, index: int = 0) -> pd.DataFrame:
    """Retrieve GWAS colocalisation data for a study locus.

    Args:
        study_locus_id: Open Targets study locus ID (e.g. ``"5a86bfd40d2ebecf6ce97bbe8a737512"``).
        size: Page size.
        index: Page index.

    Returns:
        DataFrame with colocalisation results including H3, H4, and CLPP values.
    """
    data = query_api(_GWAS_COLOC_QUERY, {"studyLocusId": study_locus_id, "size": size, "index": index})
    rows = (data.get("credibleSet") or {}).get("colocalisation", {}).get("rows", [])
    if not rows:
        return pd.DataFrame()
    df = flatten_to_df(rows)
    # Rename to match R package output column names
    rename_map = {
        "otherStudyLocus.study.id": "study.studyId",
        "otherStudyLocus.study.projectId": "study.projectId",
        "otherStudyLocus.study.traitFromSource": "study.traitReported",
        "otherStudyLocus.study.publicationFirstAuthor": "study.publicationFirstAuthor",
        "otherStudyLocus.variant.id": "indexVariant.id",
        "otherStudyLocus.variant.chromosome": "indexVariant.chromosome",
        "otherStudyLocus.variant.position": "indexVariant.position",
        "otherStudyLocus.variant.referenceAllele": "indexVariant.referenceAllele",
        "otherStudyLocus.variant.alternateAllele": "indexVariant.alternateAllele",
        "otherStudyLocus.pValueMantissa": "pValueMantissa",
        "otherStudyLocus.pValueExponent": "pValueExponent",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    # Drop nested column if present
    df = df.drop(columns=["otherStudyLocus.studyLocusId"], errors="ignore")
    return df


# ---------------------------------------------------------------------------
# gwas_credible_set (legacy Genetics API)
# ---------------------------------------------------------------------------

_SEARCH_RSID_QUERY = """
query ConvertRSIDtoVID($queryString: String!) {
  search(queryString: $queryString) { totalVariants variants { id } }
}
"""

_GWAS_CREDSET_QUERY = """
query credsetQuery($studyId: String!, $variantId: String!) {
  gwasCredibleSet(studyId: $studyId, variantId: $variantId) {
    tagVariant { id rsId }
    beta postProb pval se MultisignalMethod logABF is95 is99
  }
}
"""


def gwas_credible_set(study_id: str, variant_id: str) -> pd.DataFrame:
    """Retrieve GWAS credible set for a study and lead variant (legacy Genetics API).

    Args:
        study_id: GWAS study ID (e.g. ``"GCST006614"``).
        variant_id: Variant ID or rsID (e.g. ``"rs12345"`` or ``"1_55053079_C_T"``).
    """
    # Convert rsID to variant ID if needed
    if variant_id.startswith("rs"):
        search_data = query_api(_SEARCH_RSID_QUERY, {"queryString": variant_id}, endpoint=GENETICS_API)
        variants = search_data.get("search", {}).get("variants", [])
        if not variants:
            raise ValueError(f"Could not resolve rsID: {variant_id}")
        variant_id = variants[0]["id"]

    data = query_api(_GWAS_CREDSET_QUERY, {"studyId": study_id, "variantId": variant_id}, endpoint=GENETICS_API)
    rows = data.get("gwasCredibleSet", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# locus2gene_query
# ---------------------------------------------------------------------------

_L2G_QUERY = """
query Locus2GeneQuery($studyLocusId: String!) {
  credibleSet(studyLocusId: $studyLocusId) {
    l2GPredictions {
      count
      rows {
        shapBaseValue
        features { shapValue value name }
        score
        target { id approvedSymbol }
      }
    }
  }
}
"""


def locus2gene_query(study_locus_id: str) -> pd.DataFrame:
    """Retrieve locus-to-gene (L2G) predictions for a study locus.

    Args:
        study_locus_id: Study locus ID (e.g. ``"fa375739ca2a6b825ce5cc69d117e84b"``).
    """
    data = query_api(_L2G_QUERY, {"studyLocusId": study_locus_id})
    rows = (data.get("credibleSet") or {}).get("l2GPredictions", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# overlap_info_for_study (legacy Genetics API)
# ---------------------------------------------------------------------------

_OVERLAP_QUERY = """
query overlapinfostudyquery($studyId: String!, $studyIds: [String!]!) {
  overlapInfoForStudy(studyId: $studyId, studyIds: $studyIds) {
    study { studyId traitReported traitCategory }
    overlappedVariantsForStudies {
      overlaps { variantIdA variantIdB overlapAB distinctA distinctB }
      study { studyId traitReported traitCategory }
    }
    variantIntersectionSet
  }
}
"""


def overlap_info_for_study(study_id: str, study_ids: list = None) -> dict:
    """Retrieve variant overlap info between studies (legacy Genetics API).

    Args:
        study_id: Primary GWAS study ID.
        study_ids: List of study IDs to compare against.

    Returns:
        Raw dict with ``study``, ``overlappedVariantsForStudies``, and ``variantIntersectionSet``.
    """
    if study_ids is None:
        study_ids = []
    data = query_api(_OVERLAP_QUERY, {"studyId": study_id, "studyIds": study_ids}, endpoint=GENETICS_API)
    return data.get("overlapInfoForStudy", {})


# ---------------------------------------------------------------------------
# shared_trait_studies_query
# ---------------------------------------------------------------------------

_SHARED_TRAIT_QUERY = """
query SharedTraitStudiesQuery($diseaseIds: [String!]!, $size: Int!, $index: Int!) {
  studies(diseaseIds: $diseaseIds, page: { size: $size, index: $index }) {
    count
    rows {
      id traitFromSource projectId
      diseases { id name }
      publicationFirstAuthor publicationDate publicationJournal nSamples cohorts
      ldPopulationStructure { ldPopulation relativeSampleSize }
      pubmedId
    }
  }
}
"""


def shared_trait_studies_query(disease_ids: list, size: int = 500, index: int = 0) -> pd.DataFrame:
    """Retrieve shared-trait GWAS studies for disease IDs.

    Args:
        disease_ids: List of EFO disease IDs.
        size: Page size.
        index: Page index.
    """
    data = query_api(_SHARED_TRAIT_QUERY, {"diseaseIds": disease_ids, "size": size, "index": index})
    rows = (data.get("studies") or {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# variants_query
# ---------------------------------------------------------------------------

_VARIANTS_QUERY = """
query VariantsQuery($studyLocusId: String!, $size: Int!, $index: Int!) {
  credibleSet(studyLocusId: $studyLocusId) {
    studyLocusId
    locus(page: { size: $size, index: $index }) {
      count
      rows {
        logBF posteriorProbability
        variant {
          id chromosome position referenceAllele alternateAllele
          mostSevereConsequence { id label }
        }
        pValueMantissa pValueExponent beta standardError r2Overall
      }
    }
  }
}
"""


def variants_query(study_locus_id: str, size: int = 500, index: int = 0) -> pd.DataFrame:
    """Retrieve variant-level locus data for a study locus.

    Args:
        study_locus_id: Study locus ID.
        size: Page size.
        index: Page index.
    """
    data = query_api(_VARIANTS_QUERY, {"studyLocusId": study_locus_id, "size": size, "index": index})
    rows = (data.get("credibleSet") or {}).get("locus", {}).get("rows", [])
    return flatten_to_df(rows)
