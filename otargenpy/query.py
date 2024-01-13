#!/usr/bin/env python3

import requests
import json
import pandas as pd
import pandas as pd
import requests
import json
import re

def convert_symbol2ensembl(gene_ids):
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
    match_result = [bool(re.match(r'ENSG\d{11}', gene)) for gene in gene_ids]
    df_id = pd.DataFrame()

    if all(match_result) is False:
        for g in gene_ids:
            variables = {"queryString": g}
            response = requests.post(api_url, json={"query": query_search, "variables": variables})
            data = response.json()

            if "data" in data and "search" in data["data"]:
                genes_data = data["data"]["search"]["genes"]
                id_result = pd.DataFrame(genes_data)

                if not id_result.empty:
                    name_match = id_result[id_result["symbol"] == g]

                    if not name_match.empty:
                        ensembl_ids = name_match["id"].tolist()
                        df_id = pd.concat([df_id, pd.DataFrame({"ensembl_ids": ensembl_ids})], ignore_index=True)

        if df_id.empty:
            raise ValueError("\nPlease provide Ensemble gene ID or gene name")
        else:
            ensembl_ids = df_id["ensembl_ids"].tolist()
    else:
        ensembl_ids = gene_ids

    return ensembl_ids


def colocalisationsForGene(genes):
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
        r = requests.post(api_url, json={"query": query_coloc, "variables": variables})

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
                        # Wrap non-dictionary data in a list to create a DataFrame
                        coloc_row = pd.concat([coloc_row, k_df], axis=1)

                coloc_df = pd.concat([coloc_df, coloc_row], ignore_index=True)

            # Now coloc_df contains the flattened data in the desired format
            g_info = pd.concat([gene_info_df] * len(coloc_df), ignore_index=True)

            # Reset the index of coloc_df and g_info before concatenation
            coloc_df.reset_index(drop=True, inplace=True)
            g_info.reset_index(drop=True, inplace=True)

            # Concatenate coloc_df and g_info column-wise
            coloc_all = pd.concat([coloc_df, g_info], axis=1)

            # Reset the index of coloc_final before concatenation
            coloc_final.reset_index(drop=True, inplace=True)

            # Concatenate coloc_all with coloc_final
            coloc_final = pd.concat([coloc_final, coloc_all], ignore_index=True)

    # Now coloc_final should contain the final concatenated dataframe
    coloc_final.columns = ['variant_id', 'rsId', 'studyId', 'traitReported', 'pubJournal', 'pubTitle',
                       'pubAuthor', 'hasSumstats', 'nInitial', 'nReplication', 'nCases',
                       'numAssocLoci', 'pubDate', 'pmid', 'tissue_name', 'phenotypeId', 'h3', 'h4',
                       'log2h4h3', 'qtlStudyId', 'gene_id', 'symbol', 'description', 'chromosome',
                       'start', 'end']

    return coloc_final

# Example usage
# genes = ["ENSG00000163946", "ENSG00000169174", "ENSG00000143001"]
# result = colocalisationsForGene(genes)
# print(result)


def geneInfo(gene):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    # GraphQL query for gene name search
    query_search = """
    query convertnametoid($queryString:String!) {
        search(queryString:$queryString){
            genes{
                id
                symbol
            }
        }
    }"""

    # GraphQL query for gene information
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

    # Set base URL of Genetics Portal GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # Check gene format and convert gene symbol to Ensembl ID if necessary
    if not gene.startswith("ENSG"):
        response = requests.post(base_url, json={'query': query_search, 'variables': {'queryString': gene}})
        data = json.loads(response.text)
        gene_info = data['data']['search']['genes']
        if gene_info:
            matched_genes = [g for g in gene_info if g['symbol'] == gene]
            gene_input = matched_genes[0]['id'] if matched_genes else None
        else:
            raise ValueError("Please provide Ensemble gene ID or gene name")
    else:
        gene_input = gene

    if not gene_input:
        raise ValueError("Gene not found")

    # Execute the gene information query
    print("Downloading data...")
    response = requests.post(base_url, json={'query': query_gene_info, 'variables': {'geneId': gene_input}})
    gene_info = json.loads(response.text)['data']['geneInfo']

    # Convert to DataFrame
    output = pd.DataFrame([gene_info]) if gene_info else pd.DataFrame()

    return output

# Example usage
# result = geneInfo(gene="ENSG00000169174")
# print(result)


def genesForVariant(variant_id):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    # GraphQL endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # Function to convert rs id to variant id, if necessary
    def convert_variant_id(rs_id):
        query_searchid = """
        query ConvertRSIDtoVID($queryString:String!) {
            search(queryString:$queryString){
                totalVariants
                variants{
                    id
                }
            }
        }"""

        response = requests.post(base_url, json={'query': query_searchid, 'variables': {'queryString': rs_id}})
        data = response.json()
        return data['data']['search']['variants'][0]['id']

    # Check variant id format
    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    elif not "_" in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    # Query for variant-to-gene information
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

    variables = {'variantId': variant_id}
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    data = response.json()['data']['genesForVariant']

    # Process response data into a DataFrame
    result_df = pd.json_normalize(data)

    # Process nested data if required...
    # Example: result_qtls = pd.json_normalize(data, record_path=['qtls'], errors='ignore')

    # Return the DataFrame or a dictionary of DataFrames
    return result_df

# Example usage
# result = genesForVariant(variant_id="1_154453788_C_T")
# print(result)


def get_loci_genes(chromosome, start, end):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    # GraphQL endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # GraphQL query for fetching gene data
    query = """
    query genesquery($chromosome: String!, $start: Long!, $end: Long!){
        genes(chromosome: $chromosome, start: $start, end: $end){
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

    # Variables for the GraphQL query
    variables = {
        'chromosome': chromosome,
        'start': start,
        'end': end
    }

    # Execute the query
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    data = response.json()['data']['genes']

    # Process response data into a DataFrame
    output = pd.json_normalize(data)

    # Check if the DataFrame is empty
    if output.empty:
        print("No data found for the specified range.")
        return None

    return output

# Example usage
# result = get_loci_genes(chromosome="1", start=1000000, end=2000000)
# print(result)


def gwasColocalisation(study_id, variant_id):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # Define the function to convert rs id to variant id
    def convert_variant_id(rs_id):
        query_searchid = """
        query ConvertRSIDtoVID($queryString: String!) {
            search(queryString: $queryString) {
                totalVariants
                variants {
                    id
                }
            }
        }"""
        response = requests.post(base_url, json={'query': query_searchid, 'variables': {'queryString': rs_id}})
        data = response.json()
        return data['data']['search']['variants'][0]['id'] if data['data']['search']['variants'] else None

    # Check variant id format
    if "rs" in variant_id:
        converted_id = convert_variant_id(variant_id)
        if converted_id:
            variant_id = converted_id
        else:
            raise ValueError("Variant ID could not be converted.")
    elif not "_" in variant_id:
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

    # Prepare the variables for the query
    variables = {'studyId': study_id, 'variantId': variant_id}

    # Execute the query
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check for errors in the response
    if 'errors' in response_json:
        print("Errors in API response:", response_json['errors'])
        return pd.DataFrame()

    # Check if data is present in the response
    if 'data' in response_json and 'gwasColocalisation' in response_json['data']:
        data = response_json['data']['gwasColocalisation']
        output = pd.json_normalize(data)
    else:
        print("No data found in response")
        print("Full API response:", response_json)
        output = pd.DataFrame()

    return output

# Example usage
# result = gwasColocalisation(study_id="GCST90025951", variant_id="9_90790168_A_G")
# print(result)


def gwasColocalisationForRegion(chromosome, start, end):
    # Validate input parameters
    if not chromosome or not start or not end:
        print("Please provide values for all the arguments: chromosome, start, and end.")
        return None

    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # GraphQL query
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

    # Variables for the query
    variables = {
        'chromosome': chromosome,
        'start': start,
        'end': end
    }

    # Execute the query
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check if data is present in the response
    if 'data' in response_json and 'gwasColocalisationForRegion' in response_json['data']:
        data = response_json['data']['gwasColocalisationForRegion']
        output = pd.json_normalize(data)
    else:
        print("No data found or error in response")
        output = pd.DataFrame()

    return output

# Example usage
# result = gwasColocalisationForRegion(chromosome="1", start=1000000, end=2000000)
# print(result)


def gwasCredibleSet(study_id, variant_id):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # Function to convert rs id to variant id
    def convert_variant_id(rs_id):
        query_searchid = """
        query ConvertRSIDtoVID($queryString: String!) {
            search(queryString: $queryString) {
                totalVariants
                variants {
                    id
                }
            }
        }"""
        response = requests.post(base_url, json={'query': query_searchid, 'variables': {'queryString': rs_id}})
        data = response.json()
        return data['data']['search']['variants'][0]['id'] if data['data']['search']['variants'] else None

    # Check variant id format
    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    elif not "_" in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    # GraphQL query
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

    # Execute the query
    variables = {'studyId': study_id, 'variantId': variant_id}
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()


    # Check if data is present in the response
    if 'data' in response_json and 'gwasCredibleSet' in response_json['data']:
        data = response_json['data']['gwasCredibleSet']
        output = pd.json_normalize(data)
    else:
        print("No data found or error in response")
        output = pd.DataFrame()

    return output

# Example usage
# result = query.gwasCredibleSet(study_id="GCST90002401", variant_id="9_90797195_C_A")
# print(result)


def gwasRegional(study_id, chromosome, start, end):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # GraphQL query
    query = """
    query gwasregionalquery($studyId: String!, $chromosome: String!, $start: Long!, $end: Long!) {
        gwasRegional(studyId: $studyId, chromosome: $chromosome, start: $start, end: $end) {
            variant {
                id
                chromosome
                position
            }
            pval
        }
    }"""

    # Variables for the query
    variables = {
        'studyId': study_id,
        'chromosome': chromosome,
        'start': start,
        'end': end
    }

    # Execute the query
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check if data is present in the response
    if 'data' in response_json and 'gwasRegional' in response_json['data']:
        data = response_json['data']['gwasRegional']
        output = pd.json_normalize(data)
        output = output.rename(columns={'variant.id': 'variant_id', 
                                        'variant.chromosome': 'chromosome', 
                                        'variant.position': 'position'})
    else:
        print("No data found or error in response")
        output = pd.DataFrame()

    return output

# Example usage
# result = query.gwasRegional(study_id="GCST90002357", chromosome="1", start= 153992685, end= 154155116)
# print(result)


def indexVariantsAndStudiesForTagVariant(variant_id, pageindex=0, pagesize=20):
    if not variant_id:
        print("Please provide a value for the variant ID argument.")
        return None

    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # Function to convert rs id to variant id, if necessary
    def convert_variant_id(rs_id):
        query_searchid = """
        query ConvertRSIDtoVID($queryString: String!) {
            search(queryString: $queryString) {
                totalVariants
                variants {
                    id
                }
            }
        }"""
        response = requests.post(base_url, json={'query': query_searchid, 'variables': {'queryString': rs_id}})
        data = response.json()
        return data['data']['search']['variants'][0]['id'] if data['data']['search']['variants'] else None

    # Check variant id format
    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    elif not "_" in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    # GraphQL query
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

    # Execute the query
    variables = {
        'variantId': variant_id,
        'pageIndex': pageindex,
        'pageSize': pagesize
    }
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check if data is present in the response
    if 'data' in response_json and 'indexVariantsAndStudiesForTagVariant' in response_json['data']:
        data = response_json['data']['indexVariantsAndStudiesForTagVariant']['associations']
        output = pd.json_normalize(data)
    else:
        print("No data found or error in response")
        output = pd.DataFrame()

    return output

# Example usage
# result = indexVariantsAndStudiesForTagVariant(variant_id="rs12740374", pageindex=1, pagesize=50)
# print(result)


def extract_score_from_list(lst):
    if lst:
        return lst[0]['score']
    else:
        return None


def manhattan(study_id, pageindex=0, pagesize=100):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # GraphQL query
    query = """
    query manhattanquery($studyId: String!, $pageIndex: Int!, $pageSize: Int!) {
        manhattan(studyId: $studyId, pageIndex: $pageIndex, pageSize: $pageSize) {
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

    # Variables for the query
    variables = {
        'studyId': study_id,
        'pageIndex': pageindex,
        'pageSize': pagesize
    }

    # Execute the query
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check if data is present in the response
    if 'data' in response_json and 'manhattan' in response_json['data']:
        data = response_json['data']['manhattan']['associations']
        output = pd.json_normalize(data)

        # Handle negative exponents in 'pvalExponent' column
        output['pvalMantissa'] = output['pvalMantissa'] * 10**output['pvalExponent'].abs()

        # Drop 'pvalExponent' column
        output.drop(columns=['pvalExponent'], inplace=True)

        # Flatten 'bestGenes' column
        output['bestGenes_gene_id'] = output['bestGenes'].apply(lambda x: x[0]['gene']['id'] if x else None)
        output['bestGenes_gene_symbol'] = output['bestGenes'].apply(lambda x: x[0]['gene']['symbol'] if x else None)
        output['bestGenes_score'] = output['bestGenes'].apply(extract_score_from_list)

        # Flatten 'bestColocGenes' column
        output['bestColocGenes_gene_id'] = output['bestColocGenes'].apply(lambda x: x[0]['gene']['id'] if x else None)
        output['bestColocGenes_gene_symbol'] = output['bestColocGenes'].apply(lambda x: x[0]['gene']['symbol'] if x else None)
        output['bestColocGenes_score'] = output['bestColocGenes'].apply(extract_score_from_list)

        # Flatten 'bestLocus2Genes' column
        output['bestLocus2Genes_gene_id'] = output['bestLocus2Genes'].apply(lambda x: x[0]['gene']['id'] if x else None)
        output['bestLocus2Genes_gene_symbol'] = output['bestLocus2Genes'].apply(lambda x: x[0]['gene']['symbol'] if x else None)
        output['bestLocus2Genes_score'] = output['bestLocus2Genes'].apply(extract_score_from_list)

        # Remove the original columns
        output = output.drop(columns=['bestGenes', 'bestColocGenes', 'bestLocus2Genes'])

    else:
        print("No data found or error in response")
        output = pd.DataFrame()

    return output

# Example usage
# result = manhattan(study_id="GCST90002357", pageindex=2, pagesize=50)
# print(result)


def overlapInfoForStudy(study_id, study_ids):
    if not study_id:
        print("Please provide a study ID.")
        return None

    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # GraphQL query
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
    }"""

    # Variables for the query
    variables = {
        'studyId': study_id,
        'studyIds': study_ids
    }

    # Execute the query
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check if data is present in the response
    if 'data' in response_json and 'overlapInfoForStudy' in response_json['data']:
        data = response_json['data']['overlapInfoForStudy']
        study_df = pd.json_normalize(data, record_path=['overlappedVariantsForStudies', 'overlaps'],
                                     meta=[['study', 'studyId'], ['study', 'traitReported'], ['study', 'traitCategory']],
                                     errors='ignore')
        variant_intersection_set = data['variantIntersectionSet']
    else:
        print("No data found or error in response")
        study_df = pd.DataFrame()
        variant_intersection_set = []

    return {'overlap_info': study_df, 'variant_intersection_set': variant_intersection_set}

# Example usage
# result = overlapInfoForStudy(study_id="GCST90002357", study_ids=["GCST90025975", "GCST90025962"])
# print(result)


def pheWAS(variant_id):
    print("Connecting to the Open Targets Genetics GraphQL API...")

    base_url = "https://api.genetics.opentargets.org/graphql"

    # Function to convert rs id to variant id, if necessary
    def convert_variant_id(rs_id):
        query_searchid = """
        query ConvertRSIDtoVID($queryString: String!) {
            search(queryString: $queryString) {
                totalVariants
                variants {
                    id
                }
            }
        }"""
        response = requests.post(base_url, json={'query': query_searchid, 'variables': {'queryString': rs_id}})
        data = response.json()
        return data['data']['search']['variants'][0]['id'] if data['data']['search']['variants'] else None

    # Check variant id format
    if "rs" in variant_id:
        variant_id = convert_variant_id(variant_id)
    elif not "_" in variant_id:
        raise ValueError("Please provide a valid variant ID.")

    # GraphQL query
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

    # Execute the query
    variables = {'variantId': variant_id}
    response = requests.post(base_url, json={'query': query, 'variables': variables})
    response_json = response.json()

    # Check if data is present in the response
    if 'data' in response_json and 'pheWAS' in response_json['data']:
        data = response_json['data']['pheWAS']['associations']
        result_df = pd.json_normalize(data)
    else:
        print("No data found or error in response")
        result_df = pd.DataFrame()

    return result_df

# Example usage
# result = pheWAS(variant_id="rs72698179")
# print(result)


def qtlColocalisationVariantQuery(study_id, variant_id):
    # Define the GraphQL query
    query = f'''
    query qtlColocalisationVariantQuery($studyId: String!, $variantId: String!) {{
      qtlColocalisation(studyId: $studyId, variantId: $variantId) {{
        qtlStudyName
        phenotypeId
        gene {{
          id
          symbol
        }}
        tissue {{
          name
        }}
        indexVariant {{
          id
          rsId
        }}
        beta
        h4
        h3
        log2h4h3
      }}
    }}
    '''

    # Define the GraphQL variables
    variables = {
        "studyId": study_id,
        "variantId": variant_id
    }

    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Send the GraphQL query
    response = requests.post(url, json={"query": query, "variables": variables})

    # Check if the request was successful
    if response.status_code == 200:
        data = response.json()["data"]["qtlColocalisation"]
        df_qtlcoloc = pd.DataFrame(data)
        return df_qtlcoloc
    else:
        # Handle errors here
        print(f"Error: {response.status_code}")
        return None

# Example usage:
# study_id = "GCST90002357"
# variant_id = "1_154119580_C_A"  # Provide the variant ID
# df = qtlColocalisationVariantQuery(study_id, variant_id)
# print(df)


def qtlCredibleSet(study_id, variant_id, gene, biofeature):
    # Check if arguments are empty
    if not all([study_id, variant_id, gene, biofeature]):
        print("Please provide values for all the arguments: study_id, variant_id, gene, and biofeature.")
        return None

    # Define the GraphQL query
    query = f'''
    query qtlcredsetquery($studyId: String!, $variantId: String!, $geneId: String!, $bioFeature: String!) {{
      qtlCredibleSet(studyId: $studyId, variantId: $variantId, geneId: $geneId, bioFeature: $bioFeature) {{
        tagVariant {{
          id
          rsId
        }}
        pval
        se
        beta
        postProb
        MultisignalMethod
        logABF
        is95
        is99
      }}
    }}
    '''

    # Define the GraphQL variables
    variables = {
        "studyId": study_id,
        "variantId": variant_id,
        "geneId": gene,
        "bioFeature": biofeature
    }

    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Send the GraphQL query
    response = requests.post(url, json={"query": query, "variables": variables})

    # Check if the request was successful
    if response.status_code == 200:
        data = response.json()["data"]["qtlCredibleSet"]
        df_qtl_cred = pd.DataFrame(data)
        df_qtl_cred.rename(columns=lambda x: x.replace("tagVariant.", ""), inplace=True)
        return df_qtl_cred
    else:
        print(f"Error: {response.status_code}")
        return None

# Example usage:
# study_id = "Braineac2"
# variant_id = "1_55053079_C_T"  # Provide the variant ID
# gene = "ENSG00000169174"  # Provide the gene ID or gene name
# biofeature = "SUBSTANTIA_NIGRA"  # Provide the biofeature
# df = qtlCredibleSet(study_id, variant_id, gene, biofeature)
# print(df)


def run_custom_query(variable_list, query, query_name):
    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    # Define the GraphQL request payload
    payload = {
        "query": query,
        "variables": variable_list
    }

    try:
        # Send the GraphQL query
        response = requests.post(url, json=payload, headers=headers)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()["data"]
            return data
        else:
            print(f"Error: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None

# Example usage:
# variable_list = {"key1": "value1", "key2": "value2"}  # Provide the variable list
# query = "Your GraphQL query here"  # Provide your GraphQL query
# query_name = "YourQueryName"  # Provide a query name
# result = run_custom_query(variable_list, query, query_name)
# print(result)

def studiesAndLeadVariantsForGeneByL2G(genes, l2g=None, pvalue=None, vtype=None):
    # Define the GraphQL endpoint URL
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
    if isinstance(genes, list):
        gene_list = genes
    else:
        gene_list = [genes]
        
    #print(gene_list)
    l2g_final = pd.DataFrame()
    ensembl_ids = convert_symbol2ensembl(gene_list)
    
    #print(ensembl_ids)
    for input_gene in ensembl_ids:
        variables = {"gene_id": input_gene}
        r = requests.post(api_url, json={"query": l2g_query, "variables": variables})

        query_output = json.loads(r.text)

        if query_output.get("data"):
            gene_info = query_output["data"]["geneInfo"]
            gene_info_df = pd.DataFrame([gene_info])
           # print(gene_info_df)
            l2g_dt = query_output["data"]["studiesAndLeadVariantsForGeneByL2G"]
            l2g_keys = l2g_dt[0].keys()
            #print(l2g_keys)

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
                        # Wrap non-dictionary data in a list to create a DataFrame
                        l2g_row = pd.concat([l2g_row, k_df], axis=1)

                l2g_df = pd.concat([l2g_df, l2g_row], ignore_index=True)
                #print(l2g_df.head())

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
            l2g_final.columns = ['L2G', 'Distance', 'Interaction',
       'MolecularQTL', 'Pathogenicity', 'pval', 'direction',
       'betaCI', 'betaCILower', 'betaCIUpper', 'oddsCI', 'oddsCILower',
       'oddsCIUpper', 'studyId', 'traitReported', 'traitCategory', 'pubDate',
       'pubTitle', 'pubAuthor', 'pubJournal', 'pmid', 'hasSumstats', 'nCases',
       'numAssocLoci', 'nTotal', 'traitEfos', 'variant_id', 'rsId', 'chromosome',
       'position', 'refAllele', 'altAllele', 'nearestGene',
       'nearestCodingGene', 'nearestCodingGeneDistance', 'nearestGeneDistance',
       'mostSevereConsequence', 'ensemble_id', 'gene_symbol', 'description', 'chromosome',
       'start', 'end']



    return l2g_final

# Example usage:
# gene_list = ["ENSG00000163946", "ENSG00000169174", "ENSG00000143001"]
# result = studiesAndLeadVariantsForGeneByL2G(gene_list, l2g=0.7)
# print(result)

def studyInfo(study_id):
    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    # Check if 'study_id' argument is empty or null
    if not study_id:
        print("Please provide a value for the study ID argument.")
        return None

    # Define the GraphQL query
    query = f"""
    query {{
        studyInfo(studyId: "{study_id}") {{
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
        }}
    }}
    """

    # Define the GraphQL request payload
    payload = {
        "query": query
    }

    try:
        # Send the GraphQL query
        response = requests.post(url, json=payload, headers=headers)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()["data"]["studyInfo"]

            # Convert data to a pandas DataFrame
            df = pd.DataFrame(data)

            # Replace "NULL" values with None (equivalent to NA in R)
            df = df.applymap(lambda x: None if x == "NULL" else x)

            return df

        else:
            print(f"Error: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None

# Example usage:
# study_id = "GCST90002357"
# result = studyInfo(study_id)
# if result is not None:
#     print(result)


def studyLocus2GeneTable(study_id, variant_id):
    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    try:
        # Check variant ID format and convert rsID to variant ID if needed
        if variant_id.startswith("rs"):
            query_searchid = f"""
            query {{
                search(queryString: "{variant_id}") {{
                    totalVariants
                    variants {{
                        id
                    }}
                }}
            }}
            """
            
            # Send the GraphQL query to convert rsID to variant ID
            response = requests.post(url, json={"query": query_searchid}, headers=headers)
            
            # Check if the request was successful
            if response.status_code == 200:
                data = response.json()["data"]["search"]["variants"]
                if data:
                    input_variant_id = data[0]["id"]
                else:
                    print(f"Variant with rsID {variant_id} not found.")
                    return None
            else:
                print(f"Error: {response.status_code}")
                return None
        elif re.match(r"^\d+_\d+_[ACGT]_[ACGT]$", variant_id):
            input_variant_id = variant_id
        else:
            print("Invalid variant ID format. Please provide a valid variant ID (CHRPOSITION_REFALLELE_ALTALLELE or rsID).")
            return None

        # Define the GraphQL query
        query = f"""
        query {{
            studyLocus2GeneTable(studyId: "{study_id}", variantId: "{input_variant_id}") {{
                study {{
                    studyId
                    traitReported
                }}
                variant {{
                    id
                    rsId
                }}
                rows {{
                    gene {{
                        id
                        symbol
                    }}
                    yProbaDistance
                    yProbaModel
                    yProbaMolecularQTL
                    yProbaPathogenicity
                    yProbaInteraction
                    hasColoc
                    distanceToLocus
                }}
            }}
        }}
        """

        # Define the GraphQL request payload
        payload = {
            "query": query
        }

        # Send the GraphQL query
        response = requests.post(url, json=payload, headers=headers)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json().get("data", {}).get("studyLocus2GeneTable", {})

            # Extract rows data and flatten it
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
                    "distanceToLocus": row["distanceToLocus"]
                }
                flat_rows_data.append(row_dict)

            # Create a DataFrame from the flattened data
            df_l2g = pd.DataFrame(flat_rows_data)

            return df_l2g

        else:
            print(f"Error: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None

# Example usage:
# study_id = "GCST90002357"
# variant_id = "1_154119580_C_A"
# result = studyLocus2GeneTable(study_id, variant_id)
# if result is not None:
#     print(result)

def studyVariants(study_id):
    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    # Define the GraphQL query
    query = f"""
    query {{
        manhattan(studyId: "{study_id}") {{
            associations {{
                variant {{
                    id
                    rsId
                    chromosome
                    position
                    nearestCodingGene {{
                        id
                        symbol
                    }}
                    nearestCodingGeneDistance
                }}
                pval
                credibleSetSize
                ldSetSize
                oddsRatio
                beta
            }}
        }}
    }}
    """

    # Define the GraphQL request payload
    payload = {
        "query": query
    }

    try:
        # Send the GraphQL query
        response = requests.post(url, json=payload, headers=headers)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()["data"]["manhattan"]["associations"]

            # Create a pandas DataFrame for variants
            variants_df = pd.DataFrame(data)
            
            # Split the "variant" column into "variant_id" and "variant_rsId"
            variants_df["variant_id"] = variants_df["variant"].apply(lambda x: x.get("id"))
            variants_df["variant_rsId"] = variants_df["variant"].apply(lambda x: x.get("rsId"))
            
            # Remove the original "variant" column
            variants_df.drop(columns=["variant"], inplace=True)
            
            # Rename columns
            variants_df.rename(columns={
                "chromosome": "variant_chromosome",
                "position": "variant_position",
                "nearestCodingGene.id": "nearest_coding_gene_id",
                "nearestCodingGene.symbol": "nearest_coding_gene_symbol",
                "nearestCodingGeneDistance": "variant_nearest_coding_gene_distance",
                "credibleSetSize": "credible_set_size",
                "ldSetSize": "ld_set_size",
                "oddsRatio": "odds_ratio",
                "beta": "beta"
            }, inplace=True)
            
            # Create a pandas DataFrame for genes
            genes_df = pd.DataFrame([assoc["variant"]["nearestCodingGene"] for assoc in data])
            
            # Rename columns
            genes_df.rename(columns={
                "id": "gene_id",
                "symbol": "gene_symbol"
            }, inplace=True)
            
            # Concatenate the two DataFrames
            result_df = pd.concat([variants_df, genes_df], axis=1)
            
            return result_df

        else:
            print(f"Error: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None


# Example usage:
# study_id = "GCST003155"
# result = studyVariants(study_id)
# if result is not None:
#     variants_df, genes_df = result
#     print("Variants DataFrame:")
#     print(variants_df)
#     print("Genes DataFrame:")
#     print(genes_df)

def tagVariantsAndStudiesForIndexVariant(variant_id, pageindex=0, pagesize=20):
    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    # Check variant ID format
    if variant_id.startswith("rs"):
        # Convert rsID to variant ID
        query_searchid = f"""
        query {{
            search(queryString: "{variant_id}") {{
                totalVariants
                variants {{
                    id
                }}
            }}
        }}
        """

        try:
            # Send the GraphQL query to convert rsID to variant ID
            response = requests.post(url, json={"query": query_searchid}, headers=headers)

            # Check if the request was successful
            if response.status_code == 200:
                data = response.json()["data"]["search"]["variants"]
                if data:
                    input_variant_id = data[0]["id"]
                else:
                    print(f"Variant with rsID {variant_id} not found.")
                    return None
            else:
                print(f"Error: {response.status_code}")
                return None

        except requests.exceptions.RequestException as e:
            # Handle connection timeout and other errors
            if "Timeout was reached" in str(e):
                print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
            else:
                print(f"Error: {str(e)}")
            return None

    elif re.match(r"^\d+_\d+_[ACGT]_[ACGT]$", variant_id):
        input_variant_id = variant_id
    else:
        print("Invalid variant ID format. Please provide a valid variant ID (CHRPOSITION_REFALLELE_ALTALLELE or rsID).")
        return None

    # Define the GraphQL query
    query = f"""
    query {{
        tagVariantsAndStudiesForIndexVariant(variantId: "{input_variant_id}", pageIndex: {pageindex}, pageSize: {pagesize}) {{
            associations {{
                tagVariant {{
                    id
                    chromosome
                    rsId
                    position
                }}
                study {{
                    studyId
                    traitReported
                }}
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
            }}
        }}
    }}
    """

    # Define the GraphQL request payload
    payload = {
        "query": query
    }

    try:
        # Send the GraphQL query
        response = requests.post(url, json=payload, headers=headers)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()["data"]["tagVariantsAndStudiesForIndexVariant"]["associations"]

            # Create a pandas DataFrame
            df = pd.DataFrame(data)

            # Extract and format information from the dictionaries in DataFrame columns
            df["tagVariant.id"] = df["tagVariant"].apply(lambda x: x.get("id"))
            df["tagVariant.chromosome"] = df["tagVariant"].apply(lambda x: x.get("chromosome"))
            df["tagVariant.rsId"] = df["tagVariant"].apply(lambda x: x.get("rsId"))
            df["tagVariant.position"] = df["tagVariant"].apply(lambda x: x.get("position"))
            df["study.studyId"] = df["study"].apply(lambda x: x.get("studyId"))
            df["study.traitReported"] = df["study"].apply(lambda x: x.get("traitReported"))

            # Drop unnecessary columns
            df = df.drop(columns=["tagVariant", "study"])

            return df

        else:
            print(f"Error: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None


# Example usage:
# variant_id = "1_109274968_G_T"
# pageindex = 0
# pagesize = 20
# result = tagVariantsAndStudiesForIndexVariant(variant_id, pageindex, pagesize)
# if result is not None:
#     print("Tag Variants and Studies DataFrame:")
#     print(result)


def topOverlappedStudies(study_id, pageindex=0, pagesize=20):
    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    # Check if the study ID argument is empty or null
    if not study_id:
        print("Please provide a value for the study ID argument.")
        return None

    # Define the GraphQL query
    query = f"""
    query {{
        topOverlappedStudies(studyId: "{study_id}", pageIndex: {pageindex}, pageSize: {pagesize}) {{
            study {{
                studyId
                traitReported
                traitCategory
            }}
            topStudiesByLociOverlap {{
                studyId
                study {{
                    studyId
                    traitReported
                    traitCategory
                }}
                numOverlapLoci
            }}
        }}
    }}
    """

    # Define the GraphQL request payload
    payload = {
        "query": query
    }

    try:
        # Send the GraphQL query
        response = requests.post(url, json=payload, headers=headers)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json().get("data", {}).get("topOverlappedStudies", [])

            if not data:
                print("No data found for the given study ID.")
                return None

            # Initialize lists to store extracted data
            study_ids = []
            trait_reported = []
            trait_category = []
            top_study_ids = []
            top_trait_reported = []
            top_trait_category = []
            num_overlap_loci = []

            # Extract data from the JSON response
            for item in data:
                if isinstance(item, dict):
                    study_info = item.get("study", {})
                    top_study_info = item.get("topStudiesByLociOverlap", {}).get("study", {})
                    study_ids.append(study_info.get("studyId"))
                    trait_reported.append(study_info.get("traitReported"))
                    trait_category.append(study_info.get("traitCategory"))
                    top_study_ids.append(top_study_info.get("studyId"))
                    top_trait_reported.append(top_study_info.get("traitReported"))
                    top_trait_category.append(top_study_info.get("traitCategory"))
                    num_overlap_loci.append(item.get("topStudiesByLociOverlap", {}).get("numOverlapLoci"))
                elif isinstance(item, str):
                    # Handle cases where `item` is a string (error message or other unexpected content)
                    print(f"Skipping unexpected item in the response: {item}")
                else:
                    print(f"Unexpected item type in the response: {type(item)}")

            # Create a pandas DataFrame from the extracted data
            df = pd.DataFrame({
                "study.studyId": study_ids,
                "study.traitReported": trait_reported,
                "study.traitCategory": trait_category,
                "topStudiesByLociOverlap.studyId": top_study_ids,
                "topStudiesByLociOverlap.study.studyId": top_study_ids,
                "topStudiesByLociOverlap.study.traitReported": top_trait_reported,
                "topStudiesByLociOverlap.study.traitCategory": top_trait_category,
                "topStudiesByLociOverlap.numOverlapLoci": num_overlap_loci
            })

            return df

        else:
            print(f"Error: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None
# Example usage:
# study_id = "GCST006614_3"
# pageindex = 0
# pagesize = 20
# result = topOverlappedStudies(study_id, pageindex, pagesize)
# if result is not None:
#     print("Top Overlapped Studies DataFrame:")
#     print(result)


def variantInfo(variant_id):
    # Check if the variant ID argument is empty or null
    if not variant_id:
        print("Please provide a value for the variant ID argument.")
        return None

    # Define the GraphQL endpoint URL
    url = "https://api.genetics.opentargets.org/graphql"

    # Define the GraphQL request headers
    headers = {
        "Content-Type": "application/json"
    }

    try:
        # Set up to query Open Targets Genetics API
        print("Connecting to the Open Targets Genetics GraphQL API...")
        otg_cli = requests.Session()

        # Check variant ID format
        if variant_id.startswith("rs"):
            # Convert rs ID to variant ID
            query_searchid = f"""
            query {{
                search(queryString: "{variant_id}") {{
                    totalVariants
                    variants {{
                        id
                    }}
                }}
            }}
            """

            variables = {}
            response = otg_cli.post(url, json={"query": query_searchid}, headers=headers)
            id_result = response.json()["data"]["search"]
            input_variant_id = id_result["variants"][0]["id"] if id_result["totalVariants"] > 0 else None
        elif re.match(r"\d+_\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id):
            input_variant_id = variant_id
        else:
            raise ValueError("Please provide a valid variant ID")

        # Check if the input_variant_id is null or empty
        if not input_variant_id:
            raise ValueError("There is no variant ID defined for this rsID by Open Target Genetics")

        # Define the GraphQL query
        query = f"""
        query {{
            variantInfo(variantId: "{input_variant_id}") {{
                chromosome
                position
                refAllele
                altAllele
                rsId
                chromosomeB37
                positionB37
                id
                nearestGene {{
                    id
                    symbol
                }}
                nearestGeneDistance
                nearestCodingGene {{
                    id
                    symbol
                }}
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
            }}
        }}
        """

        # Define the GraphQL request payload
        payload = {
            "query": query
        }

        # Send the GraphQL query
        response = otg_cli.post(url, json=payload, headers=headers)
        data = response.json()["data"]["variantInfo"]

        # Create a pandas DataFrame
        df = pd.DataFrame(data)

        return df

    except requests.exceptions.RequestException as e:
        # Handle connection timeout and other errors
        if "Timeout was reached" in str(e):
            print("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
        else:
            print(f"Error: {str(e)}")
        return None

# Example usage:
# variant_id = "1_109274968_G_T"
# result = variantInfo(variant_id)
# if result is not None:
#     print("Variant Information DataFrame:")
#     print(result)
