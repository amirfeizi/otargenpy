"""Shared GraphQL client for the Open Targets Platform API."""

import requests
import pandas as pd

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"
GENETICS_API = "https://api.genetics.opentargets.org/graphql"


def query_api(query: str, variables: dict, endpoint: str = PLATFORM_API) -> dict:
    """Execute a GraphQL query and return the parsed JSON data.

    Args:
        query: GraphQL query string.
        variables: Dict of query variables.
        endpoint: GraphQL API URL.

    Returns:
        The ``data`` field from the JSON response.

    Raises:
        RuntimeError: On HTTP errors or GraphQL-level errors.
    """
    resp = requests.post(
        endpoint,
        json={"query": query, "variables": variables},
        headers={"Content-Type": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    body = resp.json()
    if "errors" in body:
        msgs = "; ".join(e.get("message", str(e)) for e in body["errors"])
        raise RuntimeError(f"GraphQL error: {msgs}")
    return body.get("data", {})


def flatten_to_df(records) -> pd.DataFrame:
    """Normalize a list of dicts (possibly with nested dicts) into a flat DataFrame."""
    if not records:
        return pd.DataFrame()
    return pd.json_normalize(records, sep=".")
