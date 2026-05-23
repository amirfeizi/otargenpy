"""Gene/target-level queries (by Ensembl ID)."""

import pandas as pd
from ._client import query_api, flatten_to_df

# ---------------------------------------------------------------------------
# comp_genomics_query
# ---------------------------------------------------------------------------

_COMP_GENOMICS_QUERY = """
query CompGenomics($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    homologues {
      speciesId speciesName homologyType isHighConfidence
      targetGeneId targetGeneSymbol
      queryPercentageIdentity targetPercentageIdentity
    }
  }
}
"""


def comp_genomics_query(ensembl_id: str) -> pd.DataFrame:
    """Retrieve comparative genomics (homologues) data for a gene.

    Args:
        ensembl_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_COMP_GENOMICS_QUERY, {"ensemblId": ensembl_id})
    rows = (data.get("target") or {}).get("homologues", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# depmap_query
# ---------------------------------------------------------------------------

_DEPMAP_QUERY = """
query DepmapQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    depMapEssentiality {
      tissueName
      screens { depmapId cellLineName diseaseFromSource geneEffect expression }
    }
  }
}
"""


def depmap_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve DepMap essentiality data for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000146648"``).
    """
    data = query_api(_DEPMAP_QUERY, {"ensgId": ensg_id})
    tissues = (data.get("target") or {}).get("depMapEssentiality", [])
    if not tissues:
        return pd.DataFrame()
    # Unnest screens per tissue
    records = []
    for t in tissues:
        for s in t.get("screens", []):
            records.append({**s, "tissueName": t["tissueName"]})
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# gene_ontology_query
# ---------------------------------------------------------------------------

_GENE_ONTOLOGY_QUERY = """
query GeneOntologyQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    geneOntology { term { id label } aspect evidence geneProduct source }
  }
}
"""


def gene_ontology_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve Gene Ontology annotations for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_GENE_ONTOLOGY_QUERY, {"ensgId": ensg_id})
    rows = (data.get("target") or {}).get("geneOntology", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# genetic_constraint_query
# ---------------------------------------------------------------------------

_GENETIC_CONSTRAINT_QUERY = """
query GeneticConstraintQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id approvedSymbol
    geneticConstraint { constraintType score upperBin upperBin6 }
  }
}
"""


def genetic_constraint_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve genetic constraint metrics (pLI, LOEUF, etc.) for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_GENETIC_CONSTRAINT_QUERY, {"ensgId": ensg_id})
    target = data.get("target") or {}
    rows = target.get("geneticConstraint", [])
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["geneId"] = target.get("id")
    df["approvedSymbol"] = target.get("approvedSymbol")
    return df


# ---------------------------------------------------------------------------
# hallmarks_query
# ---------------------------------------------------------------------------

_HALLMARKS_QUERY = """
query HallmarksQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    hallmarks {
      attributes { name pmid description }
      cancerHallmarks { pmid impact description label }
    }
  }
}
"""


def hallmarks_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve cancer hallmarks data for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_HALLMARKS_QUERY, {"ensgId": ensg_id})
    target = data.get("target") or {}
    hallmarks = target.get("hallmarks") or {}
    records = []
    for a in hallmarks.get("attributes", []):
        records.append({**a, "type": "attributes", "geneId": target.get("id")})
    for c in hallmarks.get("cancerHallmarks", []):
        records.append({**c, "type": "cancerHallmarks", "geneId": target.get("id")})
    return pd.DataFrame(records) if records else pd.DataFrame()


# ---------------------------------------------------------------------------
# interactions_query
# ---------------------------------------------------------------------------

_INTERACTIONS_QUERY = """
query InteractionsSectionQuery($ensgId: String!, $sourceDatabase: InteractionSourceEnum, $index: Int = 0, $size: Int = 10) {
  target(ensemblId: $ensgId) {
    id approvedName approvedSymbol
    interactions(sourceDatabase: $sourceDatabase, page: { index: $index, size: $size }) {
      count
      rows {
        intA intABiologicalRole
        targetA { id approvedSymbol }
        speciesA { mnemonic }
        intB intBBiologicalRole
        targetB { id approvedSymbol }
        speciesB { mnemonic }
        score count sourceDatabase
        evidences {
          evidenceScore hostOrganismScientificName
          interactionDetectionMethodMiIdentifier interactionDetectionMethodShortName
          interactionIdentifier interactionTypeShortName
          participantDetectionMethodA { miIdentifier shortName }
          participantDetectionMethodB { miIdentifier shortName }
          expansionMethodShortName pubmedId
        }
      }
    }
  }
}
"""


def interactions_query(ensg_id: str, source_database: str = None, index: int = 0, size: int = 10) -> pd.DataFrame:
    """Retrieve molecular interaction data for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
        source_database: Optional filter (e.g. ``"intact"``).
        index: Page index.
        size: Page size.
    """
    variables = {"ensgId": ensg_id, "index": index, "size": size}
    if source_database is not None:
        variables["sourceDatabase"] = source_database
    data = query_api(_INTERACTIONS_QUERY, variables)
    rows = (data.get("target") or {}).get("interactions", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# known_drugs_gene_query
# ---------------------------------------------------------------------------

_KNOWN_DRUGS_GENE_QUERY = """
query KnownDrugsQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug { id name drugType mechanismsOfAction { rows { actionType targets { id } } } }
        diseases { diseaseFromSource disease { id name } }
        clinicalReports { id source url clinicalStage trialOverallStatus }
      }
    }
  }
}
"""


def known_drugs_gene_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve drug and clinical candidate data for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_KNOWN_DRUGS_GENE_QUERY, {"ensgId": ensg_id})
    rows = (data.get("target") or {}).get("drugAndClinicalCandidates", {}).get("rows", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# mouse_phenotypes_query
# ---------------------------------------------------------------------------

_MOUSE_PHENOTYPES_QUERY = """
query MousePhenotypes($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    mousePhenotypes {
      targetInModel targetInModelMgiId modelPhenotypeId modelPhenotypeLabel
      modelPhenotypeClasses { id label }
      biologicalModels { id allelicComposition geneticBackground literature }
    }
  }
}
"""


def mouse_phenotypes_query(ensembl_id: str) -> pd.DataFrame:
    """Retrieve mouse phenotype data for a gene.

    Args:
        ensembl_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_MOUSE_PHENOTYPES_QUERY, {"ensemblId": ensembl_id})
    rows = (data.get("target") or {}).get("mousePhenotypes", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# pathways_query
# ---------------------------------------------------------------------------

_PATHWAYS_QUERY = """
query PathwaysQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id approvedSymbol
    pathways { pathwayId pathway topLevelTerm }
  }
}
"""


def pathways_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve pathway data for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_PATHWAYS_QUERY, {"ensgId": ensg_id})
    rows = (data.get("target") or {}).get("pathways", [])
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ---------------------------------------------------------------------------
# pharmacogenomics_gene_query
# ---------------------------------------------------------------------------

_PHARMACOGENOMICS_GENE_QUERY = """
query PharmacogenomicsQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    pharmacogenomics {
      variantRsId genotypeId
      variantFunctionalConsequence { id label }
      haplotypeId haplotypeFromSourceId isDirectTarget
      drugs { drugId drugFromSource }
      phenotypeFromSourceId genotypeAnnotationText phenotypeText
      pgxCategory evidenceLevel studyId literature
    }
  }
}
"""


def pharmacogenomics_gene_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve pharmacogenomics data for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000141510"``).
    """
    data = query_api(_PHARMACOGENOMICS_GENE_QUERY, {"ensgId": ensg_id})
    rows = (data.get("target") or {}).get("pharmacogenomics", [])
    return flatten_to_df(rows)


# ---------------------------------------------------------------------------
# safety_query
# ---------------------------------------------------------------------------

_SAFETY_QUERY = """
query SafetyQuery($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    safetyLiabilities {
      event eventId
      biosamples { cellFormat cellLabel tissueLabel tissueId }
      effects { dosing direction }
      studies { name type description }
      datasource literature url
    }
  }
}
"""


def safety_query(ensg_id: str) -> pd.DataFrame:
    """Retrieve safety liabilities for a gene.

    Args:
        ensg_id: Ensembl gene ID (e.g. ``"ENSG00000146648"``).
    """
    data = query_api(_SAFETY_QUERY, {"ensgId": ensg_id})
    rows = (data.get("target") or {}).get("safetyLiabilities", [])
    return flatten_to_df(rows)
