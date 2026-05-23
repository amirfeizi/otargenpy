"""Gene + disease evidence queries (by Ensembl ID + EFO ID)."""

import pandas as pd
from ._client import query_api, flatten_to_df


def _evidence_query(ensembl_id, efo_id, datasource_id, query_str, size, cursor=None):
    """Shared helper for disease-evidence queries with optional cursor pagination."""
    variables = {"ensemblId": ensembl_id, "efoId": efo_id, "size": size}
    if cursor is not None:
        variables["cursor"] = cursor
        query_str = query_str.replace("{optional_cursor}", ", $cursor: String")
        query_str = query_str.replace("{cursor_param}", "cursor: $cursor")
    else:
        query_str = query_str.replace("{optional_cursor}", "")
        query_str = query_str.replace("{cursor_param}", "")
    return query_api(query_str, variables)


# ---------------------------------------------------------------------------
# chembl_query
# ---------------------------------------------------------------------------

_CHEMBL_QUERY = """
query ChemblQuery($ensemblId: String!, $efoId: String!, $size: Int!{optional_cursor}) {
  disease(efoId: $efoId) {
    id
    chembl: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["chembl"]
      size: $size {cursor_param}
    ) {
      cursor count
      rows {
        disease { id name }
        target { id approvedSymbol }
        drug {
          id name drugType
          mechanismsOfAction { rows { mechanismOfAction targets { id approvedSymbol } } }
        }
        directionOnTrait targetFromSourceId clinicalStage studyStartDate
        trialWhyStopped trialStopReasonCategories cohortPhenotypes
        urls { niceName url }
      }
    }
  }
}
"""


def chembl_query(ensembl_id: str, efo_id: str, size: int = 10, cursor: str = None) -> pd.DataFrame:
    """Retrieve ChEMBL evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID (e.g. ``"ENSG00000143799"``).
        efo_id: EFO/MONDO disease ID (e.g. ``"EFO_0000305"``).
        size: Number of records.
        cursor: Pagination cursor.
    """
    data = _evidence_query(ensembl_id, efo_id, "chembl", _CHEMBL_QUERY, size, cursor)
    rows = (data.get("disease") or {}).get("chembl", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# clinvar_query
# ---------------------------------------------------------------------------

_CLINVAR_QUERY = """
query ClinvarQuery($ensemblId: String!, $efoId: String!, $size: Int!{optional_cursor}) {
  disease(efoId: $efoId) {
    id name
    eva: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["eva"]
      size: $size {cursor_param}
    ) {
      cursor count
      rows {
        disease { id name }
        variant { id hgvsId referenceAllele alternateAllele }
        directionOnTrait diseaseFromSource variantRsId studyId
        variantFunctionalConsequence { id label }
        clinicalSignificances allelicRequirements alleleOrigins
        confidence literature cohortPhenotypes
      }
    }
  }
}
"""


def clinvar_query(ensembl_id: str, efo_id: str, size: int = 10, cursor: str = None) -> pd.DataFrame:
    """Retrieve ClinVar evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID (e.g. ``"ENSG00000012048"``).
        efo_id: EFO disease ID (e.g. ``"EFO_0001075"``).
        size: Number of records.
        cursor: Pagination cursor.
    """
    data = _evidence_query(ensembl_id, efo_id, "eva", _CLINVAR_QUERY, size, cursor)
    rows = (data.get("disease") or {}).get("eva", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# europe_pmc_query
# ---------------------------------------------------------------------------

_EUROPE_PMC_QUERY = """
query EuropePMCQuery($ensemblId: String!, $efoId: String!, $size: Int!{optional_cursor}) {
  disease(efoId: $efoId) {
    id
    europePmc: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true size: $size datasourceIds: ["europepmc"]
      {cursor_param}
    ) {
      count cursor
      rows {
        disease { name id }
        target { approvedSymbol id }
        literature
        textMiningSentences { tStart tEnd dStart dEnd section text }
        resourceScore
      }
    }
  }
}
"""


def europe_pmc_query(ensembl_id: str, efo_id: str, size: int = 50, cursor: str = None) -> pd.DataFrame:
    """Retrieve Europe PMC literature evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID.
        efo_id: EFO disease ID.
        size: Number of records.
        cursor: Pagination cursor.
    """
    data = _evidence_query(ensembl_id, efo_id, "europepmc", _EUROPE_PMC_QUERY, size, cursor)
    rows = (data.get("disease") or {}).get("europePmc", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# gene_burden_query
# ---------------------------------------------------------------------------

_GENE_BURDEN_QUERY = """
query GeneBurdenQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    geneBurdenSummary: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["gene_burden"] size: $size
    ) {
      count
      rows {
        disease { id name }
        diseaseFromSource target { id approvedSymbol }
        releaseVersion targetFromSourceId urls { url }
        directionOnTrait allelicRequirements studyId ancestry ancestryId
        resourceScore cohortId projectId statisticalMethod statisticalMethodOverview
        studyCases studyCasesWithQualifyingVariants studySampleSize
        oddsRatio oddsRatioConfidenceIntervalLower oddsRatioConfidenceIntervalUpper
        beta betaConfidenceIntervalLower betaConfidenceIntervalUpper
        pValueMantissa pValueExponent literature
      }
    }
  }
}
"""


def gene_burden_query(ensembl_id: str, efo_id: str, size: int = 3500) -> pd.DataFrame:
    """Retrieve gene burden evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID.
        efo_id: EFO disease ID.
        size: Number of records.
    """
    data = query_api(_GENE_BURDEN_QUERY, {"ensemblId": ensembl_id, "efoId": efo_id, "size": size})
    rows = (data.get("disease") or {}).get("geneBurdenSummary", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# genomics_england_query
# ---------------------------------------------------------------------------

_GENOMICS_ENGLAND_QUERY = """
query GenomicsEnglandQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name
    genomicsEngland: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["genomics_england"] size: $size
    ) {
      count
      rows {
        disease { id name }
        target { approvedSymbol }
        diseaseFromSource cohortPhenotypes confidence allelicRequirements
        studyOverview studyId literature
      }
    }
  }
}
"""


def genomics_england_query(ensembl_id: str, efo_id: str, size: int = 3500) -> pd.DataFrame:
    """Retrieve Genomics England evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID.
        efo_id: EFO disease ID.
        size: Number of records.
    """
    data = query_api(_GENOMICS_ENGLAND_QUERY, {"ensemblId": ensembl_id, "efoId": efo_id, "size": size})
    rows = (data.get("disease") or {}).get("genomicsEngland", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# gwas_credible_sets_query
# ---------------------------------------------------------------------------

_GWAS_CREDIBLE_SETS_QUERY = """
query GwasCredibleSetsQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) { approvedSymbol }
  disease(efoId: $efoId) {
    id name
    gwasCredibleSets: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["gwas_credible_sets"] size: $size
    ) {
      count
      rows {
        disease { id name }
        credibleSet {
          studyLocusId
          study { traitFromSource id projectId publicationFirstAuthor publicationDate pubmedId nSamples }
          variant { id chromosome position referenceAllele alternateAllele }
          pValueMantissa pValueExponent beta finemappingMethod confidence
        }
        score
      }
    }
  }
}
"""


def gwas_credible_sets_query(ensembl_id: str, efo_id: str, size: int = 500) -> pd.DataFrame:
    """Retrieve GWAS credible sets evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID (e.g. ``"ENSG00000169174"``).
        efo_id: EFO disease ID (e.g. ``"EFO_0004911"``).
        size: Number of records.
    """
    data = query_api(_GWAS_CREDIBLE_SETS_QUERY, {"ensemblId": ensembl_id, "efoId": efo_id, "size": size})
    disease = data.get("disease") or {}
    rows = disease.get("gwasCredibleSets", {}).get("rows", [])
    if not rows:
        return pd.DataFrame()
    df = flatten_to_df(rows)
    df["targetEnsemblId"] = ensembl_id
    df["targetSymbol"] = (data.get("target") or {}).get("approvedSymbol")
    df["diseaseId"] = disease.get("id")
    df["diseaseName"] = disease.get("name")
    return df


# ---------------------------------------------------------------------------
# orphanet_query
# ---------------------------------------------------------------------------

_ORPHANET_QUERY = """
query OrphanetQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    orphanetSummary: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["orphanet"] size: $size
    ) {
      count
      rows {
        target { id approvedSymbol }
        disease { id name }
        directionOnTrait diseaseFromSource diseaseFromSourceId diseaseFromSourceMappedId
        targetFromSource targetFromSourceId alleleOrigins confidence literature
        variantFunctionalConsequence { id label }
      }
    }
  }
}
"""


def orphanet_query(ensembl_id: str, efo_id: str, size: int = 3500) -> pd.DataFrame:
    """Retrieve Orphanet evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID.
        efo_id: EFO disease ID.
        size: Number of records.
    """
    data = query_api(_ORPHANET_QUERY, {"ensemblId": ensembl_id, "efoId": efo_id, "size": size})
    rows = (data.get("disease") or {}).get("orphanetSummary", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# uniprot_literature_query
# ---------------------------------------------------------------------------

_UNIPROT_LIT_QUERY = """
query UniprotLiteratureQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    uniprotLiteratureSummary: evidences(
      ensemblIds: [$ensemblId] enableIndirect: true datasourceIds: ["uniprot_literature"] size: $size
    ) {
      count
      rows {
        disease { id name }
        diseaseFromSource targetFromSourceId studyId literature confidence
      }
    }
  }
}
"""


def uniprot_literature_query(ensembl_id: str, efo_id: str, size: int = 3500) -> pd.DataFrame:
    """Retrieve UniProt literature evidence for a gene–disease pair.

    Args:
        ensembl_id: Ensembl gene ID.
        efo_id: EFO disease ID.
        size: Number of records.
    """
    data = query_api(_UNIPROT_LIT_QUERY, {"ensemblId": ensembl_id, "efoId": efo_id, "size": size})
    rows = (data.get("disease") or {}).get("uniprotLiteratureSummary", {}).get("rows", [])
    return flatten_to_df(rows)
