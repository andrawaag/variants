__author__ = 'andra'

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
from pprint import pprint
import requests
from SPARQLWrapper import SPARQLWrapper, JSON
import ProteinBoxBot_Core.wdi_core as PBB_Core
import ProteinBoxBot_Core.PBB_Debug as PBB_Debug
import ProteinBoxBot_Core.PBB_login as PBB_login
import ProteinBoxBot_Core.PBB_settings as PBB_settings
from time import gmtime, strftime
import copy

logincreds = PBB_login.WDLogin("ProteinBoxBot", "sNxvAlNtjQ24")

# chromosomes
#Chromosomes dict
chromosomes = dict()
chromosomes['1'] = "Q430258"
chromosomes['2'] = "Q638893"
chromosomes['3'] = "Q668633"
chromosomes['4'] = "Q836605"
chromosomes['5'] = "Q840741"
chromosomes['6'] = "Q540857"
chromosomes['7'] = "Q657319"
chromosomes['8'] = "Q572848"
chromosomes['9'] = "Q840604"
chromosomes['10'] = "Q840737"
chromosomes['11'] = "Q847096"
chromosomes['12'] = "Q847102"
chromosomes['13'] = "Q840734"
chromosomes['14'] = "Q138955"
chromosomes['15'] = "Q765245"
chromosomes['16'] = "Q742870"
chromosomes['17'] = "Q220677"
chromosomes['18'] = "Q780468"
chromosomes['19'] = "Q510786"
chromosomes['20'] = "Q666752"
chromosomes['21'] = "Q753218"
chromosomes['22'] = "Q753805"
chromosomes['22'] = "Q753805"
chromosomes['X'] = "Q61333"
chromosomes['Y'] = "Q202771"
chromosomes["MT"] = "Q27075"



# Reference section
# Prepare references
variant_id = "201"
refStatedIn = PBB_Core.WDItemID(value="Q27612411", prop_nr='P248', is_reference=True)
timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
refReferenceURL = PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/variants/"+variant_id, prop_nr="P854", is_reference=True)
variant_reference = [refStatedIn, refRetrieved, refReferenceURL]

genomeBuildQualifier = PBB_Core.WDItemID(value="Q21067546", prop_nr='P659', is_qualifier=True)


r = requests.get('https://civic.genome.wustl.edu/api/variants/'+variant_id)
variant_data = r.json()
pprint(variant_data)

prep = dict()

sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")

ncbi_geneQuery = """SELECT ?item ?itemLabel
    WHERE
    {
        ?item wdt:P351 """
ncbi_geneQuery += "\""+str(variant_data["entrez_id"])
ncbi_geneQuery += """\" .
        SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
    }
    """
print(ncbi_geneQuery)
sparql.setQuery(ncbi_geneQuery)
sparql.setReturnFormat(JSON)
results = sparql.query().convert()
for result in results["results"]["bindings"]:
   prep['P361'] = [PBB_Core.WDItemID(value=result["item"]["value"].replace("http://www.wikidata.org/entity/", ""), prop_nr='P361', references=[copy.deepcopy(variant_reference)])]
   print("wd entrez:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))


# part of
partof = variant_data["entrez_id"]

# variant_id
prep['P3329'] = [PBB_Core.WDString(value=variant_id, prop_nr='P3329', references=[copy.deepcopy(variant_reference)])]


#coordinates
coordinates = variant_data["coordinates"]
if coordinates["chromosome"] != None:
    prep['P1057'] = [PBB_Core.WDItemID(value=chromosomes[coordinates["chromosome"]], prop_nr='P1057', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)])]
    if coordinates["chromosome2"] != None:
        prep['P1057'].append(PBB_Core.WDItemID(value=chromosomes[coordinates["chromosome2"]], prop_nr='P1057', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)]))

    # genomic start
    prep['P644'] = [PBB_Core.WDString(value=str(coordinates["start"]), prop_nr='P644', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)])]
    prep['P645'] = [PBB_Core.WDString(value=str(coordinates["stop"]), prop_nr='P645', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)])]

    if coordinates["start2"] != None:
        prep['P644'].append(PBB_Core.WDString(value=str(coordinates["start2"]), prop_nr='P644', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)]))
        prep['P645'].append(PBB_Core.WDString(value=str(coordinates["stop2"]), prop_nr='P645', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)]))

query = """
       SELECT DISTINCT ?item  ?itemLabel
        WHERE
        {
        	?item ?p ?o .
            ?o prov:wasDerivedFrom ?t .
            ?t ?u wd:Q7452458 .
            SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
        }
        """
sparql.setQuery(query)
sparql.setReturnFormat(JSON)
results = sparql.query().convert()
seqO = dict()
for result in results["results"]["bindings"]:
    seqO[result["itemLabel"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
    print("wd disease:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))

prep["P31"] = []
for variant_type in variant_data["variant_types"]:
    if variant_type["name"] != "N/A":
        prep['P31'].append(PBB_Core.WDItemID(value=seqO[variant_type["display_name"]], prop_nr='P31', references=[copy.deepcopy(variant_reference)]))

    print(variant_type["display_name"])

name = variant_data["variant_aliases"]
hgvs_expressions = variant_data["hgvs_expressions"]

print(name)

for evidence_item in variant_data["evidence_items"]:
    if evidence_item["evidence_level"] == "A" or evidence_item["evidence_level"] == "B" or evidence_item["evidence_level"] == "C" :
        ## Disease
        if evidence_item["disease"]["doid"] != None:
            print("DOID:"+evidence_item["disease"]["doid"])

            query = """SELECT ?item ?itemLabel
                WHERE
                {
                    ?item wdt:P699 """
            query += "\"DOID:"+evidence_item["disease"]["doid"]
            query += """\" .
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
                }
                """
            sparql.setQuery(query)
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()
            for result in results["results"]["bindings"]:
              print("wd disease:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
              disease = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
        else:
            disease = ""

        ## Drugs
        print(evidence_item["drugs"])
        drug_qualifiers = []
        wd_drugs = []
        for drug in evidence_item["drugs"]:
            query = """
            SELECT ?item ?itemLabel
            WHERE
            {
	            ?item wdt:P31 wd:Q12140 ;
                      rdfs:label ?label .
	            FILTER regex(?label,  \""""
            query += drug["name"]
            query += """\", \"i\")
            }
            """
            print(query)
            sparql.setQuery(query)
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()
            for result in results["results"]["bindings"]:
                print("wd drug:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                wd_drugs.append(result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                drug_qualifier = PBB_Core.WDItemID(value=result["item"]["value"].replace("http://www.wikidata.org/entity/", ""), prop_nr='P2176', is_qualifier=True)
                drug_qualifiers.append(drug_qualifier)

        ## Pubmed
        print("PMID: "+evidence_item["source"]["pubmed_id"])
        pubmedQuery = """SELECT ?item ?itemLabel
            WHERE
            {
                ?item wdt:P698 """
        pubmedQuery += "\""+evidence_item["source"]["pubmed_id"]
        pubmedQuery += """\" .
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
            }
            """

        print(pubmedQuery)
        sparql.setQuery(pubmedQuery)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        pubmed_references = []

        for result in results["results"]["bindings"]:
           pubmed_entry = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
           print("wd source:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
           refStatedIn = PBB_Core.WDItemID(value=result["item"]["value"].replace("http://www.wikidata.org/entity/", ""), prop_nr='P248', is_reference=True)
           diseferences.append(refStatedIn)

        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
        pubmed_references.append(refRetrieved)
        refReferenceURL = PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/variants/"+variant_id, prop_nr="P854", is_reference=True)
        pubmed_references.append(refReferenceURL)

        # Positive therapeutic predictor
        if evidence_item["evidence_type"] ==  "Predictive" and evidence_item["clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Supports":
            disease_qualifier = [PBB_Core.WDItemID(value=disease, prop_nr='P2175', is_qualifier=True)]
            for drug in wd_drugs:
                if 'P3354' not in prep.keys():
                    prep["P3354"] = []
                prep['P3354'].append(PBB_Core.WDItemID(value=drug, prop_nr='P3354', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(disease_qualifier)))

        # Positive diagnostic predictor (stated in)
        if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Positive" and evidence_item["evidence_direction"] == "Supports":
            if 'P3356' not in prep.keys():
                prep["P3356"] = []
            prep['P3356'].append(PBB_Core.WDItemID(value=disease, prop_nr='P3356', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(drug_qualifiers)))

        # Positive diagnostic predictor (disputed by)
        if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Positive" and evidence_item["evidence_direction"] == "Does Not Support":
            disputed_by_qualifier = [PBB_Core.WDItemID(value=pubmed_entry, prop_nr='P1310', is_qualifier=True)]


            if 'P3356' not in prep.keys():
                prep["P3356"] = []
            prep['P3356'].append(PBB_Core.WDItemID(value=disease, prop_nr='P3356', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(disputed_by_qualifier)))

        # Negative diagnostic predictor (stated in)
        if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Negative" and evidence_item["evidence_direction"] == "Supports":
            if 'P3357' not in prep.keys():
                prep["P3357"] = []
            prep['P3357'].append(PBB_Core.WDItemID(value=disease, prop_nr='P3357', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(drug_qualifiers)))

        # Negative diagnostic predictor (disputed by)
        if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Negative  " and evidence_item["evidence_direction"] == "Does Not Support":
            disputed_by_qualifier = [PBB_Core.WDItemID(value=pubmed_entry, prop_nr='P1310', is_qualifier=True)]


            if 'P3357' not in prep.keys():
                prep["P3357"] = []
            prep['P3357'].append(PBB_Core.WDItemID(value=disease, prop_nr='P3357', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(disputed_by_qualifier)))


        # Positive prognostic predictor
        if evidence_item["evidence_type"] ==  "Prognostic" and evidence_item["clinical_significance"] == "Better Outcome" and evidence_item["evidence_direction"] == "Supports":
            if 'P3358' not in prep.keys():
                prep["P3358"] = []
            prep['P3358'].append(PBB_Core.WDItemID(value=disease, prop_nr='P3358', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(drug_qualifiers)))


        # Negative prognostic predictor
        if evidence_item["evidence_type"] ==  "Prognostic" and evidence_item["clinical_significance"] == "Poor Outcome" and evidence_item["evidence_direction"] == "Supports":
            if 'P3359' not in prep.keys():
                prep["P3359"] = []
            prep['P3359'].append(PBB_Core.WDItemID(value=disease, prop_nr='P3359', references=[copy.deepcopy(pubmed_references)] , qualifiers=copy.deepcopy(drug_qualifiers)))

data2add = []
for key in prep.keys():
    for statement in prep[key]:
        data2add.append(statement)
        print(statement.prop_nr, statement.value)

pprint(prep)
name = variant_data["name"]
wdPage = PBB_Core.WDItemEngine( item_name=name, data=data2add, server="www.wikidata.org", domain="genes")
wdPage.set_label(name, "en")
wdPage.set_description("genetic variant", "en")

wd_json_representation = wdPage.get_wd_json_representation()

pprint(wd_json_representation)

print(wdPage.write(logincreds))
