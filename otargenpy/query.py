#!/usr/bin/env python3
import requests
import json
import pandas as pd
import pandas as pd
import requests
import json
import re


def convert_symbol2ensembl(gene_ids):
    """Convert gene symbols to Ensembl gene IDs using the Open Targets Genetics
    GraphQL API.

    This function takes a list of gene identifiers (gene symbols or Ensembl
    gene IDs) and returns a list of Ensembl gene IDs. If a gene symbol is
    provided, it queries the Open Targets Genetics API to retrieve the
    corresponding Ensembl ID.  If an Ensembl gene ID is already provided, it is
    returned as is.

    Parameters
    ----------
    gene_ids : list of str
        A list of gene identifiers. Each identifier can be a gene symbol or an
        Ensembl gene ID (formatted as 'ENSG' followed by 11 digits).

    Returns
    -------
    list of str
        A list containing the Ensembl gene IDs corresponding to the input gene
        identifiers.

    Raises
    ------
    ValueError
        If no valid Ensembl gene IDs or gene symbols are provided.

    Examples
    --------
    >>> convert_symbol2ensembl(["BRAF", "GPR35"])
    ['ENSG00000157764', 'ENSG00000157764']

    Notes
    -----
    This function requires an internet connection to access the
    Open Targets Genetics
    GraphQL API.
    It also depends on the 'requests' and 'pandas' libraries.

    The function checks if the input identifiers are already in the
    Ensembl gene ID
    format. If not, it performs an API call for each gene symbol to fetch the
    corresponding Ensembl gene ID. If no match is found for a given gene symbol,
    or if an invalid identifier is provided, a ValueError is raised.
    """

    # Define the GraphQL API endpoint
    api_url = "https://api.genetics.opentargets.org/graphql"

    # Initialize a GraphQL query
    query_search = """
    query gene2ensembl($queryString: String!) {
        search(queryString: $queryString) {
            genes {
                id
                symbol
            }
        }
    }
    """

    # Check format
    match_result = [bool(re.match(r"ENSG\d{11}", gene)) for gene in gene_ids]
    df_id = pd.DataFrame()

    if all(match_result) is False:
        for g in gene_ids:
            variables = {"queryString": g}
            response = requests.post(
                api_url, json={"query": query_search, "variables": variables}
            )
            data = response.json()

            if "data" in data and "search" in data["data"]:
                genes_data = data["data"]["search"]["genes"]
                id_result = pd.DataFrame(genes_data)

                if not id_result.empty:
                    name_match = id_result[id_result["symbol"] == g]

                    if not name_match.empty:
                        ensembl_ids = name_match["id"].tolist()
                        df_id = pd.concat(
                            [df_id,
                             pd.DataFrame({"ensembl_ids": ensembl_ids})],
                            ignore_index=True,
                        )

        if df_id.empty:
            raise ValueError("\nPlease provide Ensemble gene ID or gene name")
        else:
            ensembl_ids = df_id["ensembl_ids"].tolist()
    else:
        ensembl_ids = gene_ids

    return ensembl_ids

def fetch_gene_colocs(genes):
    """This function performs the same query as `colocalisationsForGene` from
    the OTG's GraphQL schema. It retrieve colocalisation data for a list of
    genes using the Open Targets Genetics GraphQL API.

    This function fetches colocalisation information for each gene in the input
    list. It returns a DataFrame containing various details about
    colocalisations such as variant ID, study ID, publication details, and gene
    information.

    Parameters
    ----------
    genes : list of str     A list of gene identifiers
    (gene symbols or Ensembl gene IDs).

    Returns
    -------
    DataFrame:A pandas DataFrame containing colocalisation
    data and associated gene information. Columns include     'variant_id',
    'rsId', 'studyId', 'traitReported', 'pubJournal', 'pubTitle', 'pubAuthor',
    'hasSumstats',     'nInitial', 'nReplication', 'nCases', 'numAssocLoci',
    'pubDate', 'pmid', 'tissue_name', 'phenotypeId',     'h3', 'h4',
    'log2h4h3', 'qtlStudyId', 'gene_id', 'symbol', 'description', 'chromosome',
    'start', 'end'.

    Raises
    ------
    HTTPError     If the request to the Open Targets Genetics API
    fails.

    Notes
    -----
    The function requires an active internet connection to access
    the Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries, and also utilizes the `convert_symbol2ensembl` function for
    converting gene symbols to Ensembl gene IDs.
    """
    # Define the GraphQL API endpoint
    api_url = "https://api.genetics.opentargets.org/graphql"

    # GraphQL query for colocalizations
    query_coloc = """
    query geneandcolocal($gene: String!) {
        geneInfo(geneId: $gene) {
            id
            symbol
            description
            chromosome
            start
            end
    }

    colocalisationsForGene(geneId: $gene) {
        leftVariant {
            id
            rsId
        }
        study {
            studyId
            traitReported
            pubJournal
            pubTitle
            pubAuthor
            hasSumstats
            nInitial
            nReplication
            nCases
            numAssocLoci
            pubDate
            pmid
        }
        tissue {
            name
        }
            phenotypeId
            h3
            h4
            log2h4h3
            qtlStudyId
        }
    }
    """

    coloc_final = pd.DataFrame()
    ensembl_ids = convert_symbol2ensembl(genes)
    for input_gene in ensembl_ids:
        variables = {"gene": input_gene}
        r = requests.post(
            api_url,
            json={
                "query": query_coloc,
                "variables": variables})

        query_output = json.loads(r.text)

        if query_output.get("data"):
            gene_info = query_output["data"]["geneInfo"]
            gene_info_df = pd.DataFrame([gene_info])
            coloc_dt = query_output["data"]["colocalisationsForGene"]
            coloc_keys = coloc_dt[0].keys()

            coloc_df = pd.DataFrame()
            for el in coloc_dt:
                coloc_row = pd.DataFrame()
                for k in coloc_keys:
                    if isinstance(el[k], dict):
                        el_dt = el[k]  # dict
                        k_df = pd.DataFrame([el_dt])
                        coloc_row = pd.concat([coloc_row, k_df], axis=1)
                    else:
                        k_df = pd.DataFrame([{k: el[k]}])
                        # Wrap non-dictionary data in a list to create a
                        # DataFrame
                        coloc_row = pd.concat([coloc_row, k_df], axis=1)

                coloc_df = pd.concat([coloc_df, coloc_row], ignore_index=True)

            # Now coloc_df contains the flattened data in the desired format
            g_info = pd.concat(
                [gene_info_df] * len(coloc_df),
                ignore_index=True)

            # Reset the index of coloc_df and g_info before concatenation
            coloc_df.reset_index(drop=True, inplace=True)
            g_info.reset_index(drop=True, inplace=True)

            # Concatenate coloc_df and g_info column-wise
            coloc_all = pd.concat([coloc_df, g_info], axis=1)

            # Reset the index of coloc_final before concatenation
            coloc_final.reset_index(drop=True, inplace=True)

            # Concatenate coloc_all with coloc_final
            coloc_final = pd.concat(
                [coloc_final, coloc_all], ignore_index=True)

    # Now coloc_final should contain the final concatenated dataframe
    coloc_final.columns = [
        "variant_id",
        "rsId",
        "studyId",
        "traitReported",
        "pubJournal",
        "pubTitle",
        "pubAuthor",
        "hasSumstats",
        "nInitial",
        "nReplication",
        "nCases",
        "numAssocLoci",
        "pubDate",
        "pmid",
        "tissue_name",
        "phenotypeId",
        "h3",
        "h4",
        "log2h4h3",
        "qtlStudyId",
        "gene_id",
        "symbol",
        "description",
        "chromosome",
        "start",
        "end",
    ]

    return coloc_final

def fetch_gene_info(gene):
    """Retrieve detailed information for a specific gene from the Open Targets
    Genetics GraphQL API.

    This function performs the same query as `geneInfo` from the OTG's
    GraphQL schema.
    It queries detailed information about a gene, given a gene symbol
    or an Ensembl gene ID.

    Parameters
    ----------
    gene : str
        A gene symbol or an Ensembl gene ID.

    Returns
    -------
    DataFrame
        A pandas DataFrame with detailed information about the gene.

    Raises
    ------
    ValueError
        If the gene symbol is invalid or the gene ID is not found.

    Examples
    --------
    >>> result = fetch_gene_info(gene="ENSG00000169174")
    >>> print(result)
    # This will print the DataFrame containing detailed information for
    # the specified gene.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API.
    It depends on the 'requests' and 'pandas' libraries.
    """
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"
    query_search = """
    query convertnametoid($queryString:String!) {
        search(queryString:$queryString){
            genes{
                id
                symbol
            }
        }
    }"""
    query_gene_info = """
    query geneInfoQuery($geneId: String!){
        geneInfo(geneId:$geneId){
            id
            symbol
            bioType
            description
            chromosome
            tss
            start
            end
            fwdStrand
            exons
        }
    }"""
    if not gene.startswith("ENSG"):
        response = requests.post(
            base_url, json={"query": query_search,
                            "variables": {"queryString": gene}}
        )
        data = response.json()
        gene_info = data["data"]["search"]["genes"]
        if gene_info:
            matched_genes = [g for g in gene_info if g["symbol"] == gene]
            gene_input = matched_genes[0]["id"] if matched_genes else None
        else:
            raise ValueError(
                "Please provide Ensemble gene ID or     gene name")
    else:
        gene_input = gene
    if not gene_input:
        raise ValueError("Gene not found")

    response = requests.post(
        base_url, json={"query": query_gene_info,
                        "variables": {"geneId": gene_input}}
    )
    if response.status_code == 200:
        gene_info = response.json()["data"]["geneInfo"]
        output = pd.DataFrame([gene_info]) if gene_info else pd.DataFrame()
        return output
    else:
        raise ValueError(
            "Error retrieving gene information. Status Code: {}".format(
                response.status_code
            )
        )

def convert_variant_id(rs_id):
    """Convert an rsID to a variant ID using the Open Targets Genetics GraphQL
    API.

    Parameters ---------- rs_id : str     The rsID to be converted.

    Returns ------- str     The converted variant ID.
    """
    query_searchid = """
        query ConvertRSIDtoVID($queryString:String!) {
            search(queryString:$queryString){
                totalVariants
                variants{
                    id
                }
            }
        }"""

    base_url = "https://api.genetics.opentargets.org/graphql"
    response = requests.post(
        base_url, json={"query": query_searchid,
                        "variables": {"queryString": rs_id}}
    )
    if response.status_code != 200:
        raise ValueError("Error converting rsID to variant ID.")
    data = response.json()
    variants = data.get("data", {}).get("search", {}).get("variants", [])
    if not variants:
        raise ValueError("No variants found for the given rsID.")

    return variants[0]["id"]

def map_variant2genes(variant_id):
    """Retrieve gene information associated with a specific variant from the
    Open Targets Genetics GraphQL API.

    This function performs the same query as `genesForVariant`
    from the OTG's GraphQL
    schema. It retrieves detailed
    information about genes associated with a given variant,
    either specified by an
    rsID or a variant ID.

    Parameters
    ----------
    variant_id : str
        A variant identifier, either an rsID or a variant ID.

    Returns
    -------
    DataFrame
        A pandas DataFrame with detailed information about the genes
        associated with the specified variant.

    Raises
    ------
    ValueError
        If the variant ID is invalid or not found.

    Examples
    --------
    >>> result = map_variant2genes(variant_id="rs28362263")
    >>> print(result)
    # This will print the DataFrame containing gene information for the
    # specified rsID or variant ID.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It relies on the
    'requests' and 'pandas' libraries for API communication and data processing.
    """
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    elif "_" not in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    query = """
    query v2gquery($variantId: String!){
        genesForVariant(variantId: $variantId) {
            gene{
                id
                symbol
            }
            variant
            overallScore
            # Other fields as required...
        }
    }"""

    response = requests.post(
        base_url, json={"query": query, "variables": {"variantId": variant_id}}
    )
    if response.status_code != 200:
        raise ValueError("Error retrieving variant-to-gene information.")

    data = response.json().get("data", {}).get("genesForVariant", [])
    return pd.json_normalize(data)

def fetch_gwas_coloc(study_id, variant_id):
    """Retrieve GWAS colocalisation data from the Open Targets Genetics GraphQL
    API.

    This function performs the same query as `gwasColocalisation` from the OTG's
    GraphQL schema. It queries GWAS colocalisation data based on a given
    study ID and variant ID.

    Parameters
    ----------
    study_id : str
        The study identifier.
    variant_id : str
        The variant identifier, either an rsID or a variant ID.

    Returns
    -------
    DataFrame
        A pandas DataFrame with GWAS colocalisation data.

    Raises
    ------
    ValueError
        If the variant ID is invalid, not found, or could not be converted.

    Examples
    --------
    >>> result = fetch_gwas_coloc(study_id="GCST90025951",
    variant_id="9_90790168_A_G")
    >>> print(result)
    # This will print the DataFrame containing GWAS colocalisation data for
    # the specified study ID and variant ID.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"
    # Check variant id format
    if "rs" in variant_id:
        converted_id = convert_variant_id(variant_id)
        if converted_id:
            variant_id = converted_id
        else:
            raise ValueError("Variant ID could not be converted.")
    elif "_" not in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    query = """
    query gwascolquery($studyId: String!, $variantId: String!) {
        gwasColocalisation(studyId: $studyId, variantId: $variantId) {
            indexVariant {
                id
                position
                chromosome
                rsId
            }
            study {
                studyId
                traitReported
                traitCategory
            }
            beta
            h3
            h4
            log2h4h3
        }
    }"""

    response = requests.post(
        base_url,
        json={
            "query": query,
            "variables": {"studyId": study_id, "variantId": variant_id},
        },
    )
    response_json = response.json()

    if "errors" in response_json:
        raise ValueError("Errors in API response:", response_json["errors"])

    data = response_json.get("data", {}).get("gwasColocalisation", [])
    return pd.json_normalize(data) if data else pd.DataFrame()

def fetch_gwas_coloc_region(chromosome, start, end):
    """Retrieve GWAS colocalisation data for a specific genomic region from the
    Open Targets Genetics GraphQL API.

    This function performs the same query as `gwasColocalisationForRegion`
    from the OTG's GraphQL schema. It queries GWAS colocalisation data
    based on a specified chromosome region defined by start and end positions.

    Parameters
    ----------
    chromosome : str
        The chromosome identifier.
    start : int
        The start position on the chromosome.
    end : int
        The end position on the chromosome.

    Returns
    -------
    DataFrame
        A pandas DataFrame with GWAS colocalisation data for the specified region.

    Raises
    ------
    ValueError
        If the input parameters are invalid or missing.

    Examples
    --------
    >>> result = fetch_gwas_coloc_region(chromosome="1", start=1000000, end=2000000)
    >>> print(result)
    # This will print the DataFrame containing GWAS colocalisation data for the specified region.

    Notes
    -----
    The function requires an active internet connection to access the Open Targets Genetics API.
    It depends on the 'requests' and 'pandas' libraries.
    """
    if not chromosome or not start or not end:
        raise ValueError(
            "Please provide values for all the arguments: chromosome, start, and end."
        )

    print("Connecting to the Open Targets Genetics GraphQL API...")
    base_url = "https://api.genetics.opentargets.org/graphql"

    query = """
    query gwasColForReg_query($chromosome: String!, $start: Long!, $end: Long!) {
        gwasColocalisationForRegion(chromosome: $chromosome, start: $start, end: $end) {
            leftVariant {
                id
                position
                chromosome
                rsId
            }
            leftStudy {
                studyId
                traitReported
                traitCategory
            }
            rightVariant {
                id
                position
                chromosome
                rsId
            }
            rightStudy {
                studyId
                traitReported
                traitCategory
            }
            h3
            h4
            log2h4h3
        }
    }"""

    variables = {"chromosome": chromosome, "start": start, "end": end}
    response = requests.post(
        base_url,
        json={
            "query": query,
            "variables": variables})
    response_json = response.json()

    if (
        "data" in response_json
        and "gwasColocalisationForRegion" in response_json["data"]
    ):
        data = response_json["data"]["gwasColocalisationForRegion"]
        return pd.json_normalize(data)
    else:
        print("No data found or error in response")
        return pd.DataFrame()

def fetch_gwas_credset(study_id, variant_id):
    """Retrieve GWAS credible set data from the Open Targets Genetics GraphQL
    API.

    This function queries GWAS credible set data based on a given study ID and
    variant ID.

    Parameters
    ----------
    study_id : str
        The study identifier.
    variant_id : str
        The variant identifier, either an rsID or a variant ID.

    Returns
    -------
    DataFrame
        A pandas DataFrame with GWAS credible set data.

    Raises
    ------
    ValueError
        If the variant ID is invalid, not found, or could not be converted.

    Examples
    --------
    >>> result = fetch_gwas_credset(
    ...     study_id="GCST90002401", variant_id="9_90797195_C_A"
    ... )
    >>> print(result)
    # This will print the DataFrame containing GWAS credible set data
    # for the specified study ID and variant ID.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API.
    It depends on the 'requests' and 'pandas' libraries.
    """
    print("Connecting to the Open Targets Genetics GraphQL API...")
    base_url = "https://api.genetics.opentargets.org/graphql"

    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    if "_" not in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    query = """
    query credsetQuery($studyId: String!, $variantId: String!) {
        gwasCredibleSet(studyId: $studyId, variantId: $variantId) {
            tagVariant {
                id
                rsId
            }
            beta
            postProb
            pval
            se
            MultisignalMethod
            logABF
            is95
            is99
        }
    }"""

    variables = {"studyId": study_id, "variantId": variant_id}
    response = requests.post(
        base_url,
        json={
            "query": query,
            "variables": variables})
    response_json = response.json()

    if "data" in response_json and "gwasCredibleSet" in response_json["data"]:
        data = response_json["data"]["gwasCredibleSet"]
        return pd.json_normalize(data)
    else:
        print("No data found or error in response")
        return pd.DataFrame()

def fetch_gwas_regional(study_id, chromosome, start, end):
    """Retrieve GWAS regional data from the Open Targets Genetics GraphQL API.

    This function queries GWAS regional data based on a study ID and a specific
    chromosome region defined by start and end positions.

    Parameters
    ----------
    study_id : str
        The study identifier.
    chromosome : str
        The chromosome identifier.
    start : int
        The start position on the chromosome.
    end : int
        The end position on the chromosome.

    Returns
    -------
    DataFrame
        A pandas DataFrame with GWAS regional data.

    Raises
    ------
    ValueError
        If the input parameters are invalid or missing.

    Examples
    --------
    >>> result = fetch_gwas_regional(
    ...     study_id="GCST90002357", chromosome="1", start=153992685,
    end=154155116
    ... )
    >>> print(result)
    # This will print the DataFrame containing GWAS regional data for the
    # specified parameters.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    print("Connecting to the Open Targets Genetics GraphQL API...")
    base_url = "https://api.genetics.opentargets.org/graphql"

    query = """
    query gwasregionalquery($studyId: String!, $chromosome: String!,
    $start: Long!, $end: Long!) {
        gwasRegional(studyId: $studyId, chromosome: $chromosome,
        start: $start, end: $end) {
            variant {
                id
                chromosome
                position
            }
            pval
        }
    }"""

    variables = {
        "studyId": study_id,
        "chromosome": chromosome,
        "start": start,
        "end": end,
    }
    response = requests.post(base_url,
                             json={"query": query, "variables": variables})
    response_json = response.json()

    if "data" in response_json and "gwasRegional" in response_json["data"]:
        data = response_json["data"]["gwasRegional"]
        output = pd.json_normalize(data)
        output = output.rename(
            columns={
                "variant.id": "variant_id",
                "variant.chromosome": "chromosome",
                "variant.position": "position",
            }
        )
        return output
    else:
        print("No data found or error in response")
        return pd.DataFrame()

def fetch_tagvariant_assoc(variant_id, page_index=0, page_size=20):
    """Retrieve index variants and associated studies for a given tag variant
    from the Open Targets Genetics GraphQL API.

    This function performs the same query as
    `indexVariantsAndStudiesForTagVariant`
    from the OTG's GraphQL schema. It queries the associations
    between a tag variant and index variants, including details of
    related studies.

    Parameters
    ----------
    variant_id : str
        The tag variant identifier, either an rsID or a variant ID.
    page_index : int, optional
        Page index for the query results, by default 0.
    page_size : int, optional
        Number of results per page, by default 20.

    Returns
    -------
    DataFrame
        A pandas DataFrame with index variants and associated studies.

    Raises
    ------
    ValueError
        If the variant ID is invalid, not found, or could not be converted.

    Examples
    --------
    >>> result = fetch_tagvariant_assoc(
    ...     variant_id="rs12740374", page_index=1, page_size=50
    ... )
    >>> print(result)
    # This will print the DataFrame containing index variants and associated
    # studies for the given tag variant.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API.It depends on the 'requests' and 'pandas'
    libraries.
    """
    if not variant_id:
        raise ValueError("Please provide a value for the variant ID argument.")

    print("Connecting to the Open Targets Genetics GraphQL API...")
    base_url = "https://api.genetics.opentargets.org/graphql"

    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    if "_" not in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    query = """
    query indexvariantsandstudiesquery($variantId: String!, $pageIndex: Int!, $pageSize: Int!) {
        indexVariantsAndStudiesForTagVariant(variantId: $variantId, pageIndex: $pageIndex, pageSize: $pageSize) {
            associations {
                indexVariant {
                    id
                    rsId
                }
                study {
                    studyId
                    traitReported
                    traitCategory
                }
                # ... other fields ...
            }
        }
    }"""

    variables = {
        "variantId": variant_id,
        "pageIndex": page_index,
        "pageSize": page_size,
    }
    response = requests.post(
        base_url,
        json={
            "query": query,
            "variables": variables})
    response_json = response.json()

    if (
        "data" in response_json
        and "indexVariantsAndStudiesForTagVariant" in response_json["data"]
    ):
        data = response_json["data"]["indexVariantsAndStudiesForTagVariant"][
            "associations"
        ]
        return pd.json_normalize(data)
    else:
        print("No data found or error in response")
        return pd.DataFrame()

def extract_score_from_list(lst):
    if lst:
        return lst[0]["score"]
    else:
        return None

def fetch_phewas(variant_id):
    """Perform a Phenome-Wide Association Study (PheWAS) for a given variant ID
    using the Open Targets Genetics GraphQL API.

    Parameters ---------- variant_id : str     The variant identifier, either
    in rsID format (e.g., "rs12345") or genomic location format (e.g.,
    "1_55053079_C_T").

    Returns ------- DataFrame     A pandas DataFrame with PheWAS results
    including p-values, beta values, Raises ------ ValueError     If the
    variant ID format is invalid or if the request to the API fails.

    Notes ----- The function requires 'requests' and 'pandas' libraries, and an
    active internet connection to access the Open Targets Genetics API.
    """
    url = "https://api.genetics.opentargets.org/graphql"

    # Check variant ID format and convert rsID to variant ID if necessary
    if re.match(r"rs\d+", variant_id):
        query_searchid = """
        query ConvertRSIDtoVID($queryString: String!) {
            search(queryString: $queryString) {
                totalVariants
                variants {
                    id
                }
            }
        }"""
        variables = {"queryString": variant_id}
        response = requests.post(
            url, json={"query": query_searchid, "variables": variables}
        )
        data = response.json()
        input_variant_id = data["data"]["search"]["variants"][0]["id"]
    elif re.match(r"\d+_\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id):
        input_variant_id = variant_id
    else:
        raise ValueError("Please provide a valid variant ID")

    # PheWAS query
    query = """
    query search($variantId: String!) {
        pheWAS(variantId: $variantId) {
            totalGWASStudies
            associations {
                pval
                beta
                oddsRatio
            study {
                studyId
                source
                pmid
                pubDate
                traitReported
                traitCategory
                }
        nTotal
            }
        }
    }"""
    variables = {"variantId": input_variant_id}
    response = requests.post(
        url,
        json={
            "query": query,
            "variables": variables})

    if response.status_code == 200:
        data = response.json()["data"]["pheWAS"]
        if len(data["associations"]) != 0:
            result_df = pd.json_normalize(data, "associations")
            return result_df
        else:
            result_df = pd.DataFrame()
            return result_df
    else:
        raise ValueError(f"Error in API request: {response.status_code}")

def fetch_manhattan_data(study_id, page_index=0, page_size=100):
    """Retrieve Manhattan plot data from the Open Targets Genetics GraphQL API.

    This function performs the same query as `manhattan` from the OTG's
    GraphQL schema. It queries Manhattan plot data based on a given study ID,
    with pagination options.

    Parameters
    ----------
    study_id : str
        The study identifier.
    page_index : int, optional
        Page index for the query results, by default 0.
    page_size : int, optional
        Number of results per page, by default 100.

    Returns
    -------
    DataFrame
        A pandas DataFrame with Manhattan plot data.

    Raises
    ------
    ValueError
        If the study ID is invalid or not found.

    Examples
    --------
    >>> result = fetch_manhattan_data(
    ...     study_id="GCST90002357", page_index=2, page_size=50
    ... )
    >>> print(result)
    # This will print the DataFrame containing Manhattan plot data for the
    # specified study ID.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API.
    It depends on the 'requests' and 'pandas' libraries.
    """
    base_url = "https://api.genetics.opentargets.org/graphql"

    query = """
    query manhattanquery($studyId: String!, $pageIndex: Int!, $pageSize: Int!) {
        manhattan(studyId: $studyId,
        pageIndex: $pageIndex,
        pageSize: $pageSize) {
            associations {
                pvalMantissa
                pvalExponent
                credibleSetSize
                ldSetSize
                totalSetSize
                variant {
                    id
                    position
                    chromosome
                    rsId
                }
                pval
                oddsRatio
                oddsRatioCILower
                oddsRatioCIUpper
                beta
                betaCILower
                betaCIUpper
                direction
                bestGenes {
                    gene {
                        id
                        symbol
                    }
                    score
                }
                bestColocGenes {
                    gene {
                        id
                        symbol
                    }
                    score
                }
                bestLocus2Genes {
                    gene {
                        id
                        symbol
                    }
                    score
                }
            }
        }
    }"""
    variables = {"studyId": study_id,
                 "pageIndex": page_index,
                 "pageSize": page_size}
    response = requests.post(base_url,
                             json={"query": query, "variables": variables})
    response_json = response.json()

    if "data" in response_json and "manhattan" in response_json["data"]:
        data = response_json["data"]["manhattan"]["associations"]
        output = pd.json_normalize(data)
        output["pvalMantissa"] *= 10 ** output["pvalExponent"].abs()
        output.drop(columns=["pvalExponent"], inplace=True)

        def extract_gene_info(genes):
            if not genes:
                return None, None, None
            gene = genes[0]["gene"]
            return gene["id"], gene["symbol"], genes[0]["score"]

        for column in ["bestGenes", "bestColocGenes", "bestLocus2Genes"]:
            output[
                [f"{column}_gene_id",
                 f"{column}_gene_symbol", f"{column}_score"]
            ] = (output[column].apply(extract_gene_info).apply(pd.Series))
            output.drop(columns=[column], inplace=True)

        return output
    else:
        print("No data found or error in response")
        return pd.DataFrame()

def fetch_overlap_info_btw_studies(study_id, study_ids=None):
    """Retrieve overlap information for a given study compared with a list of
    other studies from the Open Targets Genetics GraphQL API.

This function performs the same query as `overlapInfoForStudy` from the OTG's
GraphQL schema. It queries overlap data between a specified study and a list of
other studies.

Parameters
----------
study_id : str
    The primary study identifier.
study_ids : list of str, optional
    A list of study identifiers to compare with the primary study,
    default is None.

Returns
-------
dict
    A dictionary with two keys: 'overlap_info' containing a DataFrame of
    overlap details, and 'variant_intersection_set' containing a list of
    intersected variants.

Raises
------
ConnectionError
    If there is a connection error to the Open Targets Genetics API.
Exception
    For other types of errors.

Examples
--------
    >>> result = fetch_overlap_info_for_study(
    ...     study_id="GCST90002357", study_ids=["GCST90025975", "GCST90025962"]
    ... )
    >>> print(result["overlap_info"])
    >>> print(result["variant_intersection_set"])
    """
    query = """
    query overlapinfostudyquery($studyId: String!, $studyIds: [String!]!) {
        overlapInfoForStudy(studyId: $studyId, studyIds: $studyIds) {
            study {
                studyId
                traitReported
                traitCategory
            }
            overlappedVariantsForStudies {
                overlaps {
                    variantIdA
                    variantIdB
                    overlapAB
                    distinctA
                    distinctB
                }
                study {
                    studyId
                traitReported
                traitCategory
            }
        }
        variantIntersectionSet
    }
}
"""
    url = "https://api.genetics.opentargets.org/graphql"
    variables = {"studyId": study_id, "studyIds": study_ids}
    
    try:
        response = requests.post(url,
            json={"query": query, "variables": variables})

        if response.status_code == 200:
            data = response.json()["data"]["overlapInfoForStudy"]

            # Extracting overlap info
            overlap_data = data["overlappedVariantsForStudies"]
            df_overlap = pd.json_normalize(overlap_data, "overlaps", ["study"])
            variant_overlap = data["variantIntersectionSet"]
            study = data["study"]
            return {
                "overlap_info": df_overlap,
                "variant_intersection_set": variant_overlap,
                "study": study,
            }

        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Connection error: {str(e)}")

def fetch_qtl_coloc(study_id, variant_id):
    """Retrieve QTL colocalisation data for a specific study and variant from
    the Open Targets Genetics GraphQL API.

    This function performs the same query as `qtlColocalisationVariantQuery`
    from the OTG's GraphQL schema. It queries QTL colocalisation data based on
    a study ID and variant ID.

    Parameters
    ----------
    study_id : str
        The study identifier.
    variant_id : str
        The variant identifier.

    Returns
    -------
    DataFrame
        A pandas DataFrame with QTL colocalisation data.

    Raises
    ------
    ValueError
    If the response status code is not 200 (indicating an unsuccessful request).

    Examples
    --------
    >>> study_id = "GCST90002357"
    >>> variant_id = "1_154119580_C_A"
    >>> df = fetch_qtl_coloc(study_id, variant_id)
    >>> print(df)
    # This will print the DataFrame containing QTL colocalisation
    # data for the specified study and variant IDs.

    Notes
    -----
    The function requires an active internet connection to access
    the Open Targets Genetics API.It depends on the 'requests' and 'pandas'
    libraries.
    """
    query = """ query qtlColocalisationVariantQuery($studyId: String!,
        $variantId: String!) {
        qtlColocalisation(studyId: $studyId, variantId: $variantId) {
        qtlStudyName
        phenotypeId
        gene {
            id
            symbol
        }
        tissue {
            name
        }
        indexVariant {
            id
            rsId
        }
        beta
        h4
        h3
        log2h4h3
    }
    } """

    url = "https://api.genetics.opentargets.org/graphql"
    response = requests.post(
        url,
        json={
            "query": query,
            "variables": {"studyId": study_id, "variantId": variant_id},
        },
    )

    if response.status_code == 200:
        data = response.json()["data"]["qtlColocalisation"]
        return pd.DataFrame(data)
    else:
        raise ValueError(f"Error: {response.status_code}")

def fetch_qtl_cred_set(study_id, variant_id, gene, biofeature):
    """Retrieve QTL credible set data for specific study, variant, gene, and
    biofeature from the Open Targets Genetics GraphQL API.

    This function performs the same query as `qtlCredibleSet` from the OTG's
    GraphQL schema. It queries QTL credible set data based on a study ID,
    variant ID, gene ID, and biofeature.

    Parameters
    ----------
    study_id : str
        The study identifier.
    variant_id : str
        The variant identifier.
    gene : str
        The gene identifier.
    biofeature : str
        The biofeature.

    Returns
    -------
    DataFrame
        A pandas DataFrame with QTL credible set data.

    Raises
    ------
    ValueError
    If any of the input arguments are missing or if the request is unsuccessful.

    Examples
    --------
    >>> study_id = "Braineac2"
    >>> variant_id = "1_55053079_C_T"
    >>> gene = "ENSG00000169174"
    >>> biofeature = "SUBSTANTIA_NIGRA"
    >>> df = fetch_qtl_cred_set(study_id, variant_id, gene, biofeature)
    >>> print(df)
    # This will print the DataFrame containing QTL credible set data for the
    # specified parameters.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    if not all([study_id, variant_id, gene, biofeature]):
        raise ValueError(
            "Please provide values for all the arguments: study_id, variant_id,\
            gene, and biofeature.")

    query = """
    query qtlcredsetquery($studyId: String!,
        $variantId: String!,
        $geneId: String!,
        $bioFeature: String!) {
            qtlCredibleSet(studyId: $studyId,
                variantId: $variantId,
                geneId: $geneId,
                bioFeature: $bioFeature) {
            tagVariant {
                id
                rsId
            }
                pval
                se
                beta
                postProb
                MultisignalMethod
                logABF
                is95
                is99
    }
    }"""

    url = "https://api.genetics.opentargets.org/graphql"
    response = requests.post(
        url,
        json={
            "query": query,
            "variables": {
                "studyId": study_id,
                "variantId": variant_id,
                "geneId": gene,
                "bioFeature": biofeature,
            },
        },
    )

    if response.status_code == 200:
        data = response.json()["data"]["qtlCredibleSet"]
        df_qtl_cred = pd.DataFrame(data)
        df_qtl_cred.rename(
            columns=lambda x: x.replace(
                "tagVariant.", ""), inplace=True)
        return df_qtl_cred
    else:
        raise ValueError(f"Error: {response.status_code}")

def run_custom_query(variable_list, query, query_name):
    """Send a custom GraphQL query to the Open Targets Genetics GraphQL API.

    This function is a generic interface to run custom queries against the OTG's
    GraphQL API.

    Parameters
    ----------
    variable_list : dict
    A dictionary of variables required by the GraphQL query.
    query : str
    The GraphQL query string.
    query_name : str
    A name for the query (for readability and debugging purposes).

    Returns
    -------
    dict or None
    The query result as a dictionary if the query is successful;
    otherwise, None.

    Raises
    ------
    ValueError
    If the request status code is not 200 (indicating an unsuccessful request).

    Examples
    --------
    >>> variable_list = {"key1": "value1", "key2": "value2"}
    >>> query = "Your GraphQL query here"
    >>> query_name = "YourQueryName"
    >>> result = run_custom_query(variable_list, query, query_name)
    >>> print(result)
    # This will print the result of the custom GraphQL query.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' library.
    """
    url = "https://api.genetics.opentargets.org/graphql"
    headers = {"Content-Type": "application/json"}
    payload = {"query": query, "variables": variable_list}

    try:
        response = requests.post(url, json=payload, headers=headers)

        if response.status_code == 200:
            return response.json()["data"]
        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        error_message = (
            "Connection timeout" if "Timeout was reached" in str(e) else str(e)
        )
        raise ConnectionError(f"Connection error: {error_message}")

def fetch_l2g_model_data(genes, l2g=None, pvalue=None, vtype=None):
    """Retrieve studies and lead variants for a list of genes based on the
    Likelihood2Genes (L2G) score from the Open Targets Genetics GraphQL API.

    This function performs the same query as
    `studiesAndLeadVariantsForGeneByL2G`
    from the OTG's GraphQL schema. It queries studies and lead variants
    for genes based on L2G scores.

    Parameters
    ----------
    genes : list of str
        A list of gene identifiers (gene symbols or Ensembl gene IDs).
    l2g : float, optional
        The Likelihood2Genes score threshold, by default None.
    pvalue : float, optional
        The p-value threshold, by default None.
    vtype : str, optional
        The variant type, by default None.

    Returns
    -------
    DataFrame
    A pandas DataFrame with studies and lead variants for the specified genes.

    Raises
    ------
    ValueError
    If any of the input arguments are invalid or if the request is unsuccessful.

    Examples
    --------
    >>> genes = ["ENSG00000163946", "ENSG00000169174", "ENSG00000143001"]
    >>> result = fetch_l2g_model_data(genes, l2g=0.7)
    >>> print(result)
    # This will print the DataFrame containing studies and lead variants
    # for the specified genes.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    if not genes:
        raise ValueError("Please provide a list of gene IDs.")

    api_url = "https://api.genetics.opentargets.org/graphql"

    l2g_query = """
         query geneandl2g($gene_id: String!) {
            geneInfo(geneId: $gene_id ) {
                id
                symbol
                description
                chromosome
                start
                end
            }
            studiesAndLeadVariantsForGeneByL2G(geneId: $gene_id) {
                yProbaModel
                yProbaDistance
                yProbaInteraction
                yProbaMolecularQTL
                yProbaPathogenicity
                pval
                beta {
                    direction
                    betaCI
                    betaCILower
                    betaCIUpper
                }
                odds {
                    oddsCI
                    oddsCILower
                    oddsCIUpper
                }
                study {
                    studyId
                    traitReported
                    traitCategory
                    pubDate
                    pubTitle
                    pubAuthor
                    pubJournal
                    pmid
                    hasSumstats
                    nCases
                    numAssocLoci
                    nTotal
                    traitEfos
                }
                variant {
                    id
                    rsId
                    chromosome
                    position
                    refAllele
                    altAllele
                    nearestGene {
                        id
                    }
                    nearestCodingGene {
                        id
                    }
                    nearestCodingGeneDistance
                    nearestGeneDistance
                    mostSevereConsequence
                }
            }
        }
        """

    # Check if 'gene' argument is a list
    gene_list = genes if isinstance(genes, list) else [genes]

    l2g_final = pd.DataFrame()
    ensembl_ids = convert_symbol2ensembl(gene_list)

    # print(ensembl_ids)
    for input_gene in ensembl_ids:
        variables = {"gene_id": input_gene}
        r = requests.post(
            api_url,
            json={
                "query": l2g_query,
                "variables": variables})

        query_output = json.loads(r.text)

        if query_output.get("data"):
            gene_info = query_output["data"]["geneInfo"]
            gene_info_df = pd.DataFrame([gene_info])
            # print(gene_info_df)
            l2g_dt = query_output["data"]["studiesAndLeadVariantsForGeneByL2G"]
            l2g_keys = l2g_dt[0].keys()
            # print(l2g_keys)

            l2g_df = pd.DataFrame()
            for el in l2g_dt:
                l2g_row = pd.DataFrame()
                for k in l2g_keys:
                    if isinstance(el[k], dict):
                        el_dt = el[k]  # dict
                        k_df = pd.DataFrame([el_dt])
                        l2g_row = pd.concat([l2g_row, k_df], axis=1)
                    else:
                        k_df = pd.DataFrame([{k: el[k]}])
                        # Wrap non-dictionary data in a list to create a
                        # DataFrame
                        l2g_row = pd.concat([l2g_row, k_df], axis=1)

                l2g_df = pd.concat([l2g_df, l2g_row], ignore_index=True)
                # print(l2g_df.head())

            # Now coloc_df contains the flattened data in the desired format
            g_info = pd.concat([gene_info_df] * len(l2g_df), ignore_index=True)

            # Reset the index of coloc_df and g_info before concatenation
            l2g_df.reset_index(drop=True, inplace=True)
            g_info.reset_index(drop=True, inplace=True)

            # Concatenate coloc_df and g_info column-wise
            l2g_all = pd.concat([l2g_df, g_info], axis=1)

            # Reset the index of coloc_final before concatenation
            l2g_all.reset_index(drop=True, inplace=True)

            # Concatenate coloc_all with coloc_final
    l2g_final = pd.concat([l2g_final, l2g_all], ignore_index=True)
    l2g_final.columns = [
        "L2G",
        "Distance",
        "Interaction",
        "MolecularQTL",
        "Pathogenicity",
        "pval",
        "direction",
        "betaCI",
        "betaCILower",
        "betaCIUpper",
        "oddsCI",
        "oddsCILower",
        "oddsCIUpper",
        "studyId",
        "traitReported",
        "traitCategory",
        "pubDate",
        "pubTitle",
        "pubAuthor",
        "pubJournal",
        "pmid",
        "hasSumstats",
        "nCases",
        "numAssocLoci",
        "nTotal",
        "traitEfos",
        "variant_id",
        "rsId",
        "chromosome",
        "position",
        "refAllele",
        "altAllele",
        "nearestGene",
        "nearestCodingGene",
        "nearestCodingGeneDistance",
        "nearestGeneDistance",
        "mostSevereConsequence",
        "ensemble_id",
        "gene_symbol",
        "description",
        "chromosome",
        "start",
        "end",
    ]

    return l2g_final

def fetch_study_info(study_id):
    """Retrieve detailed information for a specific study from the Open Targets
    Genetics GraphQL API.

    This function performs the same query as `studyInfo`
    from the OTG's GraphQL schema. It queries detailed
    information about a study based on the study ID.

    Parameters
    ----------
    study_id : str
        The study identifier.

    Returns
    -------
    DataFrame
        A pandas DataFrame with detailed information about the specified study.

    Raises
    ------
    ValueError
        If the study ID is invalid or if the request is unsuccessful.

    Examples
    --------
    >>> study_id = "GCST90002357"
    >>> result = fetch_study_info(study_id)
    >>> print(result)
    # This will print the DataFrame containing detailed information for
    # the specified study ID.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API.
    It depends on the 'requests' and 'pandas' libraries.
    """
    url = "https://api.genetics.opentargets.org/graphql"
    headers = {"Content-Type": "application/json"}

    if not study_id:
        raise ValueError("Please provide a value for the study ID argument.")

    query = """
    query get_study_info($study_id:String!) {
        studyInfo(studyId: $study_id) {
            studyId
            traitReported
            source
            traitEfos
            pmid
            pubDate
            pubJournal
            pubTitle
            pubAuthor
            hasSumstats
            ancestryInitial
            ancestryReplication
            nInitial
            nReplication
            nCases
            traitCategory
            numAssocLoci
            nTotal
        }
    }
    """

    try:
        response = requests.post(
            url,
            json={"query": query, "variables": {"study_id": study_id}},
            headers=headers,
        )

        if response.status_code == 200:
            data = response.json()["data"]["studyInfo"]
            df = pd.DataFrame([data])  # Convert to DataFrame
            df = df.applymap(lambda x: None if x == "NULL" else x)
            return df
        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        error_message = (
            "Connection timeout" if "Timeout was reached" in str(e) else str(e)
        )
        raise ConnectionError(f"Connection error: {error_message}")

def fetch_study_locus2gene(study_id, variant_id):
    """Retrieve Locus2Gene (L2G) table data for a specific study and variant
    from the Open Targets Genetics GraphQL API.

    This function performs the same query as `studyLocus2GeneTable` from
    the OTG's GraphQL schema. It queries L2G data based on a study ID and
    variant ID.

    Parameters
    ----------
    study_id : str
        The study identifier.
    variant_id : str
    The variant identifier
    (either CHRPOSITION_REFALLELE_ALTALLELE format or rsID).

    Returns
    -------
    DataFrame
        A pandas DataFrame with Locus2Gene table data.

    Raises
    ------
    ValueError
        If the variant ID is invalid or if the request is unsuccessful.

    Examples
    --------
    >>> study_id = "GCST90002357"
    >>> variant_id = "1_154119580_C_A"
    >>> result = fetch_study_locus2gene(study_id, variant_id)
    >>> print(result)
    # This will print the DataFrame containing Locus2Gene table data
    # for the specified study and variant IDs.

    Notes
    -----
    The function requires an active internet connection to access
    the Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL query
    query = """
    query study_locus2gene($study_id:String!, $variant_id:String!) {
        studyLocus2GeneTable(studyId: $study_id, variantId: $variant_id) {
            study {
                studyId
                traitReported
            }
            variant {
                id
                rsId
            }
            rows {
                gene {
                    id
                    symbol
                }
                yProbaDistance
                yProbaModel
                yProbaMolecularQTL
                yProbaPathogenicity
                yProbaInteraction
                hasColoc
                distanceToLocus
            }
        }
    }
    """
    variables = {"study_id": study_id, "variant_id": variant_id}

    try:
        response = requests.post(
            url,
            json={
                "query": query,
                "variables": variables})

        if response.status_code == 200:
            data = response.json().get(
                "data", {}).get(
                "studyLocus2GeneTable", {})
            rows_data = data.get("rows", [])
            flat_rows_data = []
            for row in rows_data:
                row_dict = {
                    "gene_id": row["gene"]["id"],
                    "gene_symbol": row["gene"]["symbol"],
                    "yProbaDistance": row["yProbaDistance"],
                    "yProbaModel": row["yProbaModel"],
                    "yProbaMolecularQTL": row["yProbaMolecularQTL"],
                    "yProbaPathogenicity": row["yProbaPathogenicity"],
                    "yProbaInteraction": row["yProbaInteraction"],
                    "hasColoc": row["hasColoc"],
                    "distanceToLocus": row["distanceToLocus"],
                }
                flat_rows_data.append(row_dict)

            df_l2g = pd.DataFrame(flat_rows_data)
            return df_l2g
        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Connection error: {str(e)}")

def fetch_variants(study_id):
    """Retrieve variant associations for a specific study from the Open Targets
    Genetics GraphQL API.

    This function queries variant association data based on a given study ID.

    Parameters
    ----------
    study_id : str
        The study identifier.

    Returns
    -------
    DataFrame
        A pandas DataFrame with variant associations for the specified study.

    Raises
    ------
    ValueError
        If the request is unsuccessful.

    Examples
    --------
    >>> study_id = "GCST90002357"
    >>> result = studyVariants(study_id)
    >>> print(result)
    # This will print the DataFrame containing variant associations
    # for the specified study ID.

    Notes
    -----
    The function requires an active internet connection to access
    the Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    url = "https://api.genetics.opentargets.org/graphql"
    headers = {"Content-Type": "application/json"}

    query = """
    query study2variants($study_id:String!) {
        manhattan(studyId: $study_id) {
            associations {
                variant {
                    id
                    rsId
                    chromosome
                    position
                    nearestCodingGene {
                        id
                        symbol
                    }
                    nearestCodingGeneDistance
                }
                pval
                credibleSetSize
                ldSetSize
                oddsRatio
                beta
            }
        }
    }
    """

    variables = {"study_id": study_id}
    try:
        response = requests.post(
            url, json={"query": query, "variables": variables}, headers=headers
        )

        if response.status_code == 200:
            data = response.json()["data"]["manhattan"]["associations"]

            variants_df = pd.DataFrame(data)
            variants_df["variant_id"] = variants_df["variant"].apply(
                lambda x: x.get("id")
            )
            variants_df["variant_rsId"] = variants_df["variant"].apply(
                lambda x: x.get("rsId")
            )
            variants_df.drop(columns=["variant"], inplace=True)

            variants_df.rename(
                columns={
                    "chromosome": "variant_chromosome",
                    "position": "variant_position",
                    "nearestCodingGene.id": "nearest_coding_gene_id",
                    "nearestCodingGene.symbol": "nearest_coding_gene_symbol",
                    "nearestCodingGeneDistance": "variant_nearest_coding_gene_distance",
                    "credibleSetSize": "credible_set_size",
                    "ldSetSize": "ld_set_size",
                    "oddsRatio": "odds_ratio",
                    "beta": "beta",
                },
                inplace=True,
            )

            genes_df = pd.DataFrame(
                [assoc["variant"]["nearestCodingGene"] for assoc in data]
            )
            genes_df.rename(
                columns={
                    "id": "gene_id",
                    "symbol": "gene_symbol"},
                inplace=True)

            result_df = pd.concat([variants_df, genes_df], axis=1)
            return result_df
        else:
            raise ValueError(f"Error: {response.status_code}")
    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Connection error: {str(e)}")

def fetch_tag_variants_studies(variant_id, page_index=0, page_size=20):
    """Retrieve tag variants and associated studies for a given index variant
    from the Open Targets Genetics GraphQL API.

    This function queries tag variants and their associated studies based on
    an index variant ID.

    Parameters
    ----------
    variant_id : str
        The index variant identifier (either CHRPOSITION_REFALLELE_ALTALLELE
        format or rsID).
    page_index : int, optional
        Page index for the query results, by default 0.
    page_size : int, optional
        Number of results per page, by default 20.

    Returns
    -------
    DataFrame
        A pandas DataFrame with tag variants and associated studies.

    Raises
    ------
    ValueError
        If the variant ID is invalid or if the request is unsuccessful.

    Examples
    --------
    >>> variant_id = "1_154119580_C_A"
    >>> result = fetch_tag_variants_studies(variant_id,
    page_index=0, page_size=20)
    >>> print(result)
    # This will print the DataFrame containing tag variants and
    # associated studies for the specified index variant.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    url = "https://api.genetics.opentargets.org/graphql"

    if variant_id.startswith("rs"):
        variant_id = convert_variant_id(variant_id)
    elif not re.match(r"^\d+_\d+_[ACGT]_[ACGT]$", variant_id):
        raise ValueError(
            "Invalid variant ID format. Please provide a valid variant ID."
        )

    query = """
    query tagvariants($variant_id:String!, $pageindex:Int, $pagesize: Int) {
        tagVariantsAndStudiesForIndexVariant(variantId: $variant_id,
        pageIndex: $pageindex, pageSize: $pagesize) {
            associations {
                tagVariant {
                    id
                    chromosome
                    rsId
                    position
                }
                study {
                    studyId
                    traitReported
                }
                pval
                pvalMantissa
                pvalExponent
                nTotal
                nCases
                overallR2
                afr1000GProp
                amr1000GProp
                eas1000GProp
                eur1000GProp
                sas1000GProp
                oddsRatio
                oddsRatioCILower
                oddsRatioCIUpper
                posteriorProbability
                beta
                betaCILower
                betaCIUpper
                direction
                log10Abf
            }
        }
    }
    """

    variables = {
        "variant_id": variant_id,
        "pageindex": page_index,
        "pagesize": page_size,
    }
    try:
        response = requests.post(url,
                                 json={"query": query, "variables": variables})

        if response.status_code == 200:
            data = response.json()["data"]["tagVariantsAndStudiesForIndexVariant"][
                "associations"
            ]
            df = pd.DataFrame(data)

            df["tagVariant.id"] = df["tagVariant"].apply(lambda x: x.get("id"))
            df["tagVariant.chromosome"] = df["tagVariant"].apply(
                lambda x: x.get("chromosome"))

            df["tagVariant.rsId"] = df["tagVariant"].apply(
                lambda x: x.get("rsId"))
            df["tagVariant.position"] = df["tagVariant"].apply(
                lambda x: x.get("position"))
            df["study.studyId"] = df["study"].apply(lambda x: x.get("studyId"))
            df["study.traitReported"] = df["study"].apply(
                lambda x: x.get("traitReported")
            )
            df.drop(columns=["tagVariant", "study"], inplace=True)

            return df
        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Connection error: {str(e)}")

def fetch_top_overlapping_studies(study_id, page_index=0, page_size=20):
    """This function retrieves the top overlapped studies for a specific study from the Open
    Targets Genetics GraphQL API.

    This function performs the same query as `topOverlappedStudies` from the
    OTG's GraphQL schema. It queries top overlapped studies based on
    a given study ID.

    Parameters
    ----------
    study_id : str
        The study identifier.
    page_index : int, optional
        Page index for the query results, by default 0.
    page_size : int, optional
        Number of results per page, by default 20.

    Returns
    -------
    DataFrame
        A pandas DataFrame with the top overlapped studies.

    Raises
    ------
    ValueError
        If the study ID is invalid or if the request is unsuccessful.

    Examples
    --------
    >>> study_id = "GCST006614_3"
    >>> result = fetch_top_overlapping_studies(study_id,
    page_index=0, page_size=20)
    >>> print(result)
    # This will print the DataFrame containing top overlapped
    # studies for the specified study ID.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API.
    It depends on the 'requests' and 'pandas' libraries.
    """
    url = "https://api.genetics.opentargets.org/graphql"
    headers = {"Content-Type": "application/json"}

    if not study_id:
        raise ValueError("Please provide a value for the study ID argument.")

    query = """
    query overlapstudies($study_id:String!, $pageindex:Int, $pagesize: Int){
        topOverlappedStudies(studyId: $study_id,
        pageIndex: $pageindex, pageSize: $pagesize) {
            study {
                studyId
                traitReported
                traitCategory
            }
            topStudiesByLociOverlap {
                studyId
                study {
                    studyId
                    traitReported
                    traitCategory
                }
                numOverlapLoci
            }
        }
    }
    """

    variables = {"study_id": study_id, "pageindex": page_index,
                 "pagesize": page_size}
    try:
        response = requests.post(
            url, json={"query": query, "variables": variables},
            headers=headers
        )

        if response.status_code == 200:
            data = response.json()["data"]["topOverlappedStudies"]

            if not data:
                print("No data found for the given study ID.")
                return None

            df = pd.json_normalize(
                data,
                record_path=["topStudiesByLociOverlap"],
                meta=[
                    ["study", "studyId"],
                    ["study", "traitReported"],
                    ["study", "traitCategory"],
                ],
                errors="ignore",
            )

            return df

        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Connection error: {str(e)}")

def fetch_variant_info(variant_id):
    """Retrieve detailed information for a specific variant from the Open
    Targets Genetics GraphQL API.

    This function retrieves variant information for a
    given variant ID.

    Parameters
    ----------
    variant_id : str
    The variant identifier (either CHRPOSITION_REFALLELE_ALTALLELE format
        or rsID).

    Returns
    -------
    DataFrame
    A pandas DataFrame with detailed information about the specified variant.

    Raises
    ------
    ValueError
        If the variant ID is invalid or if the request is unsuccessful.

    Examples
    --------
    >>> variant_id = "1_109274968_G_T"
    >>> result = fetch_variant_info(variant_id)
    >>> print(result)
    # This will print the DataFrame containing detailed information for
    # the specified variant.

    Notes
    -----
    The function requires an active internet connection to access the
    Open Targets Genetics API. It depends on the 'requests' and 'pandas'
    libraries.
    """
    url = "https://api.genetics.opentargets.org/graphql"
    headers = {"Content-Type": "application/json"}

    if not variant_id:
        raise ValueError("Please provide a value for the variant ID argument.")

    print("Connecting to the Open Targets Genetics GraphQL API...")
    otg_cli = requests.Session()

    if variant_id.startswith("rs"):
        variant_id = convert_variant_id(variant_id)
    elif not re.match(r"\d+_\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id):
        raise ValueError("Please provide a valid variant ID")

    if not variant_id:
        raise ValueError("There is no variant ID defined for this rsID by\
                        Open Target Genetics"
                        )

    query = """
    query getvariantinfo($variant_id:String!){
        variantInfo(variantId: $variant_id) {
            chromosome
            position
            refAllele
            altAllele
            rsId
            chromosomeB37
            positionB37
            id
            nearestGene {
                id
                symbol
            }
            nearestGeneDistance
            nearestCodingGene {
                id
                symbol
            }
            nearestCodingGeneDistance
            mostSevereConsequence
            caddRaw
            caddPhred
            gnomadAFR
            gnomadAMR
            gnomadASJ
            gnomadEAS
            gnomadFIN
            gnomadNFE
            gnomadNFEEST
            gnomadNFENWE
            gnomadNFESEU
            gnomadNFEONF
            gnomadOTH
        }
    }
    """

    variables = {"variant_id": variant_id}

    try:
        response = otg_cli.post(
            url, json={"query": query, "variables": variables}, headers=headers
        )

        if response.status_code == 200:
            data = response.json()["data"]["variantInfo"]
            df = pd.DataFrame([data])  # Convert to DataFrame
            return df
        else:
            raise ValueError(f"Error: {response.status_code}")

    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Connection error: {str(e)}")
