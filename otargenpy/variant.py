"""Variant-level queries (by variant ID)."""

import pandas as pd
from ._client import query_api, flatten_to_df

# ---------------------------------------------------------------------------
# pharmacogenomics_variant_query
# ---------------------------------------------------------------------------

_PHARMACOGENOMICS_VARIANT_QUERY = """
query PharmacogenomicsQuery($variantId: String!) {
  variant(variantId: $variantId) {
    id referenceAllele alternateAllele
    pharmacogenomics {
      genotypeId isDirectTarget
      target { id approvedSymbol }
      drugs { drugFromSource drugId }
      phenotypeFromSourceId genotypeAnnotationText phenotypeText
      pgxCategory evidenceLevel studyId literature
    }
  }
}
"""


def pharmacogenomics_variant_query(variant_id: str) -> pd.DataFrame:
    """Retrieve pharmacogenomics data for a variant.

    Args:
        variant_id: Variant ID (e.g. ``"4_1804392_G_A"``).
    """
    data = query_api(_PHARMACOGENOMICS_VARIANT_QUERY, {"variantId": variant_id})
    rows = (data.get("variant") or {}).get("pharmacogenomics", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# qtl_credible_sets_query
# ---------------------------------------------------------------------------

_QTL_CREDIBLE_SETS_QUERY = """
query QTLCredibleSetsQuery($variantId: String!, $size: Int!, $index: Int!) {
  variant(variantId: $variantId) {
    id referenceAllele alternateAllele
    qtlCredibleSets: credibleSets(
      studyTypes: [scsqtl, sceqtl, scpqtl, sctuqtl, sqtl, eqtl, pqtl, tuqtl]
      page: { size: $size, index: $index }
    ) {
      count
      rows {
        studyLocusId pValueMantissa pValueExponent beta finemappingMethod confidence isTransQtl
        variant { id chromosome position referenceAllele alternateAllele }
        study {
          id studyType condition
          target { id approvedSymbol }
          biosample { biosampleId biosampleName }
        }
        locus(variantIds: [$variantId]) { rows { posteriorProbability } }
        locusSize: locus { count }
      }
    }
  }
}
"""


def qtl_credible_sets_query(variant_id: str, size: int = 500, index: int = 0) -> pd.DataFrame:
    """Retrieve QTL credible sets for a variant.

    Args:
        variant_id: Variant ID (e.g. ``"1_154453788_C_T"``).
        size: Page size.
        index: Page index.
    """
    data = query_api(_QTL_CREDIBLE_SETS_QUERY, {"variantId": variant_id, "size": size, "index": index})
    rows = (data.get("variant") or {}).get("qtlCredibleSets", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# uniprot_variants_query
# ---------------------------------------------------------------------------

_UNIPROT_VARIANTS_QUERY = """
query UniProtVariantsQuery($variantId: String!) {
  variant(variantId: $variantId) {
    id referenceAllele alternateAllele
    evidences(datasourceIds: ["uniprot_variants"]) {
      count
      rows {
        targetFromSourceId confidence diseaseFromSource
        disease { id name }
        literature
      }
    }
  }
}
"""


def uniprot_variants_query(variant_id: str) -> pd.DataFrame:
    """Retrieve UniProt variant evidence.

    Args:
        variant_id: Variant ID (e.g. ``"4_1804392_G_A"``).
    """
    data = query_api(_UNIPROT_VARIANTS_QUERY, {"variantId": variant_id})
    rows = (data.get("variant") or {}).get("evidences", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# variant_effect_predictor_query
# ---------------------------------------------------------------------------

_VEP_QUERY = """
query VariantEffectPredictorQuery($variantId: String!) {
  variant(variantId: $variantId) {
    id
    transcriptConsequences {
      variantConsequences { id label }
      aminoAcidChange uniprotAccessions codons
      distanceFromFootprint distanceFromTss
      target { id approvedSymbol biotype }
      impact consequenceScore transcriptIndex transcriptId
      lofteePrediction siftPrediction polyphenPrediction
    }
    referenceAllele alternateAllele
  }
}
"""


def variant_effect_predictor_query(variant_id: str) -> pd.DataFrame:
    """Retrieve variant effect predictor (VEP) transcript consequences.

    Args:
        variant_id: Variant ID (e.g. ``"1_154453788_C_T"``).
    """
    data = query_api(_VEP_QUERY, {"variantId": variant_id})
    rows = (data.get("variant") or {}).get("transcriptConsequences", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# variant_effect_query
# ---------------------------------------------------------------------------

_VARIANT_EFFECT_QUERY = """
query VariantEffectQuery($variantId: String!) {
  variant(variantId: $variantId) {
    id
    variantEffect { method assessment score assessmentFlag normalisedScore }
    referenceAllele alternateAllele
  }
}
"""


def variant_effect_query(variant_id: str) -> pd.DataFrame:
    """Retrieve variant effect scores (CADD, etc.).

    Args:
        variant_id: Variant ID (e.g. ``"1_154453788_C_T"``).
    """
    data = query_api(_VARIANT_EFFECT_QUERY, {"variantId": variant_id})
    rows = (data.get("variant") or {}).get("variantEffect", [])
    return flatten_to_df(rows) if rows else pd.DataFrame()
