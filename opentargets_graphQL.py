#!/usr/bin/env python3

# Import relevant libraries for HTTP request and JSON formatting
import requests
import json
import pandas as pd

fs="\t"

# Set study_id and variant_id variables


# Build query string

variants = []
varinfo={}
path="/data/Users/francesco/Covid19_WES/5Fold_2000/"
#var=pd.read_csv(path+"Fold_1U2U3U4U5_multi_col_improved_grading_1200_updated.tsv", dtype=object, comment="#", sep="\t")
df=pd.read_csv(path+"Fold_1U2U3U4U5_multi_col_improved_grading_2000_updated_RF_XGB_stat_whole_LORs.tsv", dtype=str, comment="#", sep="\t")
for index, row in df.iterrows():
  var=row['Chr'].strip("chr")+"_"+row['Start']+"_"+row['Ref']+"_"+row['Alt']
  if var not in variants:
    variants.append(str(var))
    varinfo[var]=[row['Log Odds-Ratio'], row['RF_nonzero'], row['XGB_nonzero']]
    
###Comment the line below for real production
#variants=["9_76328028_G_C"]
# Set base URL of Genetics Portal GraphQL API endpoint
base_url = "https://api.genetics.opentargets.org/graphql"

print ("Chr_start_Ref_alt\tSource\ttraiCategory\ttraitReported\teaf\tbeta\tse\tnTotal\tnCases\tOddsRatio\tpval\tLogOddsRatio_covid\tRF_nonzero\tXGB_nonzero")
for variant_id in variants:
  #print (variant_id)
  # Set variables object of arguments to be passed to endpoint
  variables = {"myVariantId": variant_id}
  query_string = """
  query pheWAS {
  pheWAS(variantId: "%s") {
    totalGWASStudies
    associations {
      eaf
      beta
      se
      nTotal
      nCases
      oddsRatio
      study {
        source
        pmid
        pubDate
        pubJournal
        pubTitle
        pubAuthor
        hasSumstats
        nInitial
        nReplication
        nCases
        traitCategory
        numAssocLoci
        traitReported
        }
      studyId
      pval
      }
    }
  }
  """ % (variant_id)

  # Perform POST request and check status code of response
  #r=""
  #query_string = query_string.format(variant_id)
  #print (query_string)
  #r = requests.post(base_url, json={"query": query_string, "variables": variables})
  r = requests.post(base_url, json={"query": query_string})
  #print (r.json())
  #print(r.status_code)

  # Transform API response into JSON 
  #api_response_as_json = ""
  api_response_as_json = json.loads(r.text)

  # Print first element of JSON response data
  for asso in (api_response_as_json["data"]["pheWAS"]["associations"]):
    print (variant_id+fs+asso["study"]["source"]+fs+asso["study"]["traitCategory"]+fs+asso["study"]["traitReported"]+fs+str(asso["eaf"])+fs+str(asso["beta"])+fs+str(asso["se"])+fs+str(asso["nTotal"])+fs+str(asso["nCases"])+fs+str(asso["oddsRatio"])+fs+str(asso["pval"])+fs+varinfo[variant_id][0]+fs+varinfo[variant_id][1]+fs+varinfo[variant_id][2])
