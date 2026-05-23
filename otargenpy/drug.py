"""Drug-level queries (by ChEMBL ID)."""

import pandas as pd
from ._client import query_api, flatten_to_df

# ---------------------------------------------------------------------------
# adverse_events_query
# ---------------------------------------------------------------------------

_ADVERSE_EVENTS_QUERY = """
query AdverseEventsQuery($chemblId: String!, $index: Int = 0, $size: Int = 10) {
  drug(chemblId: $chemblId) {
    id
    maxLlr: adverseEvents(page: { index: 0, size: 1 }) {
      rows { logLR }
    }
    adverseEvents(page: { index: $index, size: $size }) {
      criticalValue
      count
      rows { name count logLR meddraCode }
    }
  }
}
"""


def adverse_events_query(chembl_id: str, index: int = 0, size: int = 10) -> pd.DataFrame:
    """Retrieve adverse events for a drug.

    Args:
        chembl_id: ChEMBL drug ID (e.g. ``"CHEMBL941"``).
        index: Page index for pagination.
        size: Number of records per page.

    Returns:
        DataFrame with columns ``name``, ``count``, ``logLR``, ``meddraCode``,
        ``drugId``, ``criticalValue``.
    """
    data = query_api(_ADVERSE_EVENTS_QUERY, {"chemblId": chembl_id, "index": index, "size": size})
    drug = data.get("drug")
    if not drug:
        return pd.DataFrame()
    rows = drug.get("adverseEvents", {}).get("rows", [])
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["drugId"] = drug["id"]
    df["criticalValue"] = drug["adverseEvents"]["criticalValue"]
    return df


# ---------------------------------------------------------------------------
# indications_query
# ---------------------------------------------------------------------------

_INDICATIONS_QUERY = """
query IndicationsQuery($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    indications {
      rows {
        maxClinicalStage
        disease { id name therapeuticAreas { id name } }
        clinicalReports { id source url }
      }
      count
    }
  }
}
"""


def indications_query(chembl_id: str) -> pd.DataFrame:
    """Retrieve indications for a drug.

    Args:
        chembl_id: ChEMBL drug ID (e.g. ``"CHEMBL941"``).

    Returns:
        DataFrame with indication details including clinical stage and disease info.
    """
    data = query_api(_INDICATIONS_QUERY, {"chemblId": chembl_id})
    drug = data.get("drug")
    if not drug:
        return pd.DataFrame()
    rows = drug.get("indications", {}).get("rows", [])
    if not rows:
        return pd.DataFrame()
    df = flatten_to_df(rows)
    df["drugId"] = drug["id"]
    df["indicationsCount"] = drug["indications"]["count"]
    return df


# ---------------------------------------------------------------------------
# known_drugs_chembl_query
# ---------------------------------------------------------------------------

_KNOWN_DRUGS_CHEMBL_QUERY = """
query KnownDrugsQuery($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    indications {
      count
      rows {
        maxClinicalStage
        disease { id name }
        clinicalReports { id source url clinicalStage trialOverallStatus }
      }
    }
  }
}
"""


def known_drugs_chembl_query(chembl_id: str) -> pd.DataFrame:
    """Retrieve indications and clinical reports for a drug.

    Args:
        chembl_id: ChEMBL drug ID (e.g. ``"CHEMBL1016"``).

    Returns:
        DataFrame with indication rows including clinical report details.
    """
    data = query_api(_KNOWN_DRUGS_CHEMBL_QUERY, {"chemblId": chembl_id})
    drug = data.get("drug")
    if not drug:
        return pd.DataFrame()
    rows = drug.get("indications", {}).get("rows", [])
    if not rows:
        return pd.DataFrame()
    df = flatten_to_df(rows)
    df["drugId"] = drug["id"]
    return df


# ---------------------------------------------------------------------------
# mechanisms_of_action_query
# ---------------------------------------------------------------------------

_MECHANISMS_QUERY = """
query MechanismsOfActionSectionQuery($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    mechanismsOfAction {
      rows {
        mechanismOfAction
        targetName
        targets { id approvedSymbol }
        references { source urls }
      }
      uniqueActionTypes
      uniqueTargetTypes
    }
    parentMolecule { id name }
    childMolecules { id name }
  }
}
"""


def mechanisms_of_action_query(chembl_id: str) -> pd.DataFrame:
    """Retrieve mechanisms of action for a drug.

    Args:
        chembl_id: ChEMBL drug ID (e.g. ``"CHEMBL941"``).

    Returns:
        DataFrame with mechanism of action rows.
    """
    data = query_api(_MECHANISMS_QUERY, {"chemblId": chembl_id})
    drug = data.get("drug")
    if not drug:
        return pd.DataFrame()
    rows = drug.get("mechanismsOfAction", {}).get("rows", [])
    if not rows:
        return pd.DataFrame()
    df = flatten_to_df(rows)
    df["drugId"] = drug["id"]
    return df


# ---------------------------------------------------------------------------
# pharmacogenomics_chembl_query
# ---------------------------------------------------------------------------

_PHARMACOGENOMICS_CHEMBL_QUERY = """
query PharmacogenomicsQuery($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    pharmacogenomics {
      variantRsId genotypeId
      variantFunctionalConsequence { id label }
      target { id approvedSymbol }
      haplotypeId haplotypeFromSourceId isDirectTarget
      phenotypeFromSourceId genotypeAnnotationText phenotypeText
      pgxCategory evidenceLevel studyId literature
    }
  }
}
"""


def pharmacogenomics_chembl_query(chembl_id: str) -> pd.DataFrame:
    """Retrieve pharmacogenomics data for a drug.

    Args:
        chembl_id: ChEMBL drug ID (e.g. ``"CHEMBL1016"``).

    Returns:
        DataFrame with pharmacogenomics records.
    """
    data = query_api(_PHARMACOGENOMICS_CHEMBL_QUERY, {"chemblId": chembl_id})
    drug = data.get("drug")
    if not drug:
        return pd.DataFrame()
    rows = drug.get("pharmacogenomics", [])
    if not rows:
        return pd.DataFrame()
    df = flatten_to_df(rows)
    df["drugId"] = drug["id"]
    return df
