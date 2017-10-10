__author__ = 'andra'

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
from pprint import pprint
import requests
from SPARQLWrapper import SPARQLWrapper, JSON
from wikidataintegrator import wdi_core, wdi_login, wdi_property_store
from time import gmtime, strftime
import copy
import traceback
import time

wdi_property_store.wd_properties['P3329'] = {
        'datatype': 'string',
        'name': 'CIViC Variant ID',
        'domain': ['genes'],
        'core_id': True
    }

fast_run_base_filter = {'P3329': ''}
fast_run = True

wdi_property_store.wd_properties['P31'] = {
        'datatype': 'item',
        'name': 'instance of',
        'domain': ['therapy'],
        'core_id': False
    }

logincreds = wdi_login.WDLogin("ProteinBoxBot", os.environ['wikidataApi'])

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

# CIViC evidence scores
evidence_level = dict()
evidence_level["A"] = "Q36805652"
evidence_level["B"] = "Q36806012"
evidence_level["C"] = "Q36799701"
evidence_level["D"] = "Q36806470"
evidence_level["E"] = "Q36811327"

# CIViC trust ratings
trustratings = dict()
trustratings["1"] = "Q28045396"
trustratings["2"] = "Q28045397"
trustratings["3"] = "Q28045398"
trustratings["4"] = "Q28045383"
trustratings["5"] = "Q28045399"


# r = requests.get('https://civic.genome.wustl.edu/api/variants?count=10000')

ignore_synonym_list = [
    "AMPLIFICATION",
    "EXPRESSION",
    "DELETION",
    "LOSS",
    "LOSS-OF-FUNCTION",
    "MUTATION",
    "NUCLEAR EXPRESSION",
    "OVEREXPRESSION",
    "UNDEREXPRESSION",
    "3\' UTR MUTATION",
    "BIALLELIC INACTIVATION",
    "EWSR1-FLI1",
    "EXON 12 MUTATION",
    "EXON 9 MUTATION",
    "FRAMESHIFT TRUNCATION",
    "G12",
    "G12/G13",
    "G13D",
    "METHYLATION",
    "PHOSPHORYLATION",
    "PROMOTER HYPERMETHYLATION",
    "PROMOTER METHYLATION",
    "SERUM LEVELS",
    "TMPRSS2-ERG",
    "TRUNCATING MUTATION",
]

# Reference section
# Prepare references
variant_id = "12"
refStatedIn = wdi_core.WDItemID(value="Q27612411", prop_nr='P248', is_reference=True)
timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
refRetrieved = wdi_core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
refReferenceURL = wdi_core.WDUrl("https://civic.genome.wustl.edu/links/variants/"+variant_id, prop_nr="P854", is_reference=True)
variant_reference = [refStatedIn, refRetrieved, refReferenceURL]

## prepare log file to store missing items on combination theapries. Due to the missing identifiers in CIViC
## we need do this manuaully

file = open('/tmp/quickstatements_combination_therapy.txt', 'w')


# Qualifier section
# Prepare qualifiers

genomeBuildQualifier = [wdi_core.WDItemID(value="Q21067546", prop_nr='P659', is_qualifier=True)]
r = requests.get('https://civic.genome.wustl.edu/api/variants/'+variant_id)
variant_data = r.json()
# Set core property to look at in fast_run mode
fast_run_base_filter = {'P3329': ''}
fast_run = True

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
   prep['P3433'] = [wdi_core.WDItemID(value=result["item"]["value"].replace("http://www.wikidata.org/entity/", ""), prop_nr='P3433', references=[copy.deepcopy(variant_reference)])]
   print("wd entrez:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))

# part of
partof = variant_data["entrez_id"]

# variant_id
prep['P3329'] = [wdi_core.WDString(value=variant_id, prop_nr='P3329', references=[copy.deepcopy(variant_reference)])]

# HGVS Nomenclature
prep["P3331"] = []

for hgvs in variant_data["hgvs_expressions"]:
    prep["P3331"].append(wdi_core.WDString(value=hgvs, prop_nr='P3331', references=[copy.deepcopy(variant_reference)]))


#coordinates
coordinates = variant_data["coordinates"]
if coordinates["chromosome"] != None:
    prep['P1057'] = [wdi_core.WDItemID(value=chromosomes[coordinates["chromosome"]], prop_nr='P1057', references=[copy.deepcopy(variant_reference)], qualifiers=copy.deepcopy(genomeBuildQualifier))]
    if coordinates["chromosome2"] != None:
        prep['P1057'].append(wdi_core.WDItemID(value=chromosomes[coordinates["chromosome2"]], prop_nr='P1057', references=[copy.deepcopy(variant_reference)], qualifiers=copy.deepcopy(genomeBuildQualifier)))

    # genomic start
    prep['P644'] = [wdi_core.WDString(value=str(coordinates["start"]), prop_nr='P644', references=[copy.deepcopy(variant_reference)], qualifiers=copy.deepcopy(genomeBuildQualifier))]
    prep['P645'] = [wdi_core.WDString(value=str(coordinates["stop"]), prop_nr='P645', references=[copy.deepcopy(variant_reference)], qualifiers=copy.deepcopy(genomeBuildQualifier))]

    if coordinates["start2"] != None:
        prep['P644'].append(wdi_core.WDString(value=str(coordinates["start2"]), prop_nr='P644', references=[copy.deepcopy(variant_reference)], qualifiers=copy.deepcopy(genomeBuildQualifier)))
        prep['P645'].append(wdi_core.WDString(value=str(coordinates["stop2"]), prop_nr='P645', references=[copy.deepcopy(variant_reference)], qualifiers=copy.deepcopy(genomeBuildQualifier)))

query = """
        SELECT ?item ?itemLabel WHERE {
   ?item wdt:P3986 ?seqID .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
        """
sparql.setQuery(query)
sparql.setReturnFormat(JSON)
results = sparql.query().convert()
seqO = dict()
for result in results["results"]["bindings"]:
    seqO[result["itemLabel"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
    if "alias" in result.keys():
      if result["alias"]["value"] != "":
        seqO[result["alias"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
    seqO[result["itemLabel"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
    print("wd disease:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))

prep["P31"] = []
for variant_type in variant_data["variant_types"]:
    if variant_type["name"] == "N/A":
        continue
    prep['P31'].append(wdi_core.WDItemID(value=seqO[variant_type["display_name"]], prop_nr='P31', references=[copy.deepcopy(variant_reference)]))

    print(variant_type["display_name"])

for evidence_item in variant_data["evidence_items"]:

    if evidence_item["status"] == "accepted" and evidence_item["rating"] != None:
        print(evidence_item["id"])
        evidence_qualifiers = []
        evidence_qualifiers.append(
            wdi_core.WDItemID(value=evidence_level[str(evidence_item["evidence_level"])], prop_nr="P459",
                              is_qualifier=True))
        evidence_qualifiers.append(
            wdi_core.WDItemID(value=trustratings[str(evidence_item["rating"])], prop_nr="P1552", is_qualifier=True))

        evidence_references = []
        evidence_references.append(wdi_core.WDItemID(value="Q27612411",prop_nr="P1640", is_reference=True))
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        evidence_references.append(wdi_core.WDTime(timeStringNow, prop_nr='P813', is_reference=True))
        evidence_references.append(wdi_core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+str(evidence_item["id"]), prop_nr="P854", is_reference=True))

        statement = dict()
        disease = None
        drugs = None
        ## Disease
        if evidence_item["disease"]["doid"] != None or evidence_item["disease"]["doid"] != "":
            query = """SELECT ?item ?itemLabel
                WHERE
                {
                    ?item wdt:P699 """
            query += "\"DOID:" + evidence_item["disease"]["doid"]
            query += """\" .
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
                }
                """
            sparql.setQuery(query)
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()
            statement["disease_name"] = evidence_item["disease"]["name"]
            for result in results["results"]["bindings"]:
                #print("wd disease:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                statement["wd_disease"] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
            print("=====")
        else:
            continue

        # Drugs
        if len(evidence_item["drugs"]) > 0:
            ## Drugs
            # print(len(evidence_item["drugs"]))
            # pprint(evidence_item["drugs"])
            if len(evidence_item["drugs"]) > 1:
                drugs = []
                for drug in evidence_item["drugs"]:
                    drugs.append(drug["name"])
                combination_therapy = " / ".join(drugs) + " combination therapy"
                statement["drug_label"] = combination_therapy
                query = "SELECT * WHERE { ?item wdt:P31 wd:Q1304270 ;"
                query += "rdfs:label \""+combination_therapy+"\"@en . }"

                sparql.setQuery(query)
                sparql.setReturnFormat(JSON)
                results = sparql.query().convert()
                if len(results["results"]["bindings"]) == 0:
                    file.write('CREATE\n')
                    file.write('LAST\tLen\t\"'+combination_therapy+"\"\n")
                    file.write('LAST\tDen\t\"combination therapy\"\n')
                    file.write('LAST\tP31\tQ1304270\n\n')

                for result in results["results"]["bindings"]:
                    # print(result["item"]["value"])
                    # print("wd drug:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                    statement["wd_drug"] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
            else: # len(evidence_item["drugs"]) == 1
                #pprint(evidence_item)
                statement["drug_label"] = evidence_item["drugs"][0]["name"]
                drug_query = "SELECT DISTINCT * WHERE { {{?item rdfs:label \""+evidence_item["drugs"][0]["name"]+"\"@en } UNION {?item skos:altLabel \""+evidence_item["drugs"][0]["name"]+"\"@en}}"
                drug_query += "UNION {{?item rdfs:label \""+evidence_item["drugs"][0]["name"].lower()+"\"@en } UNION {?item skos:altLabel \""+evidence_item["drugs"][0]["name"].lower()+"\"@en}} ?item rdfs:label ?label . FILTER (lang(?label) = \"en\")}"
                # print(drug_query)
                sparql.setQuery(drug_query)
                sparql.setReturnFormat(JSON)
                results = sparql.query().convert()
                for wd_drug in results["results"]["bindings"]:
                    if wd_drug["label"]["value"].lower() == evidence_item["drugs"][0]["name"].lower():
                        # print(wd_drug["item"]["value"])
                        #print("wd drug:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                        statement["wd_drug"] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")

        ## Pubmed
        print("PMID: " + evidence_item["source"]["pubmed_id"])
        pubmedQuery = """SELECT ?item ?itemLabel
            WHERE
            {
                ?item wdt:P698 """
        pubmedQuery += "\"" + evidence_item["source"]["pubmed_id"]
        pubmedQuery += """\" .
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
            }
            """

        print(pubmedQuery)
        sparql.setQuery(pubmedQuery)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        statement["pubmed"] = []
        for result in results["results"]["bindings"]:
            pubmed_entry = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
            statement["pubmed"].append(pubmed_entry)
            refStatedIn = wdi_core.WDItemID(value=pubmed_entry, prop_nr='P248', is_reference=True)
            refDisputedBy = wdi_core.WDItemID(value=pubmed_entry, prop_nr='P1310', is_qualifier=True)
            evidence_references.append(refStatedIn)

        print(statement)

        # Positive therapeutic predictor
        if evidence_item["evidence_type"] == "Predictive" and evidence_item["clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Supports":
            temp_qualifier = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P2175", is_qualifier=True)]
            for qualifier in evidence_qualifiers:
                temp_qualifier.append(qualifier)
            evidence_qualifiers = temp_qualifier
            prep["P3354"] = [wdi_core.WDItemID(value=statement["wd_drug"], prop_nr="P3354", references=[copy.deepcopy(evidence_references)], qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Evidence does not support Positive therapeutic predictor
        if evidence_item["evidence_type"] == "Predictive" and evidence_item[
            "clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Does Not Support":
            temp_qualifier = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P2175", is_qualifier=True)]
            temp_qualifier.append(refDisputedBy)
            for qualifier in evidence_qualifiers:
                temp_qualifier.append(qualifier)
            evidence_qualifiers = temp_qualifier
            prep["P3354"] = [wdi_core.WDItemID(value=statement["wd_drug"], prop_nr="P3354", references=[copy.deepcopy(evidence_references)], qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Negative therapeutic predictor
        if evidence_item["evidence_type"] == "Resistance or Non-Response" and evidence_item[
            "clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Supports":
            temp_qualifier = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P2175", is_qualifier=True)]
            for qualifier in evidence_qualifiers:
                temp_qualifier.append(qualifier)
            evidence_qualifiers = temp_qualifier
            prep["P3355"] = [wdi_core.WDItemID(value=statement["wd_drug"], prop_nr="P3355", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Evidence does not support Negative therapeutic predictor
        if evidence_item["evidence_type"] == "Predictive" and evidence_item[
            "clinical_significance"] == "Resistance or Non-Response" and evidence_item["evidence_direction"] == "Does Not Support":
            temp_qualifier = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P2175", is_qualifier=True)]
            temp_qualifier.append(refDisputedBy)
            for qualifier in evidence_qualifiers:
                temp_qualifier.append(qualifier)
            evidence_qualifiers = temp_qualifier
            prep["P3355"] = [wdi_core.WDItemID(value=statement["wd_drug"], prop_nr="P3355", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]


        # Positive diagnostic predictor
        if evidence_item["evidence_type"] == "Diagnostic" and evidence_item["clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Supports":
            prep["P3356"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3356", references=[copy.deepcopy(evidence_references)], qualifiers=evidence_qualifiers)]

        # Evidence does not support Positive diagnostic predictor
        if evidence_item["evidence_type"] == "Diagnostic" and evidence_item[
            "clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Does Not Support":
            evidence_qualifiers.append(refDisputedBy)
            prep["P3356"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3356", references=[copy.deepcopy(evidence_references)], qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Negative diagnostic predictor
        if evidence_item["evidence_type"] == "Diagnostic" and evidence_item[
            "clinical_significance"] == "Resistance or Non-Response" and evidence_item["evidence_direction"] == "Supports":
            prep["P3357"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3357", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Evidence does not support Negative diagnostic predictor
        if evidence_item["evidence_type"] == "Diagnostic" and evidence_item[
            "clinical_significance"] == "Resistance or Non-Response" and evidence_item["evidence_direction"] == "Does Not Support":
            evidence_qualifiers.append(refDisputedBy)
            prep["P3357"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3357", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Positive prognostic predictor
        if evidence_item["evidence_type"] == "Prognositc" and evidence_item[
            "clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Supports":
            prep["P3358"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3358", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Evidence does not support Positive prognostic predictor
        if evidence_item["evidence_type"] == "Prognositc" and evidence_item[
            "clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Does Not Support":
            evidence_qualifiers.append(refDisputedBy)
            prep["P3358"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3358", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Negative prognostic predictor
        if evidence_item["evidence_type"] == "Prognostic" and evidence_item[
            "clinical_significance"] == "Resistance or Non-Response" and evidence_item["evidence_direction"] == "Supports":
            prep["P3359"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3359", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

        # Evidence does not support Negative prognostic predictor
        if evidence_item["evidence_type"] == "Prognostic" and evidence_item[
            "clinical_significance"] == "Resistance or Non-Response" and evidence_item[
            "evidence_direction"] == "Does Not Support":
            evidence_qualifiers.append(refDisputedBy)
            prep["P3359"] = [wdi_core.WDItemID(value=statement["wd_disease"], prop_nr="P3359", references=[copy.deepcopy(evidence_references)],
                                  qualifiers=copy.deepcopy(evidence_qualifiers))]

data2add = []
for key in prep.keys():
    for statement in prep[key]:
        data2add.append(statement)
        print(statement.prop_nr, statement.value)

pprint(prep)
name = variant_data["name"]
wdPage = wdi_core.WDItemEngine(item_name=name, data=data2add, server="www.wikidata.org", domain="genes")
synonyms = []
if name not in ignore_synonym_list:
    synonyms.append(name)

wdPage.set_label(variant_data["entrez_name"] + " " + name, "en")
# wdPage.set_label(variant_data["entrez_name"]+" "+ name, "nl")


if wdPage.get_description(lang='en') == "":
    wdPage.set_description("genetic variant", "en")
if wdPage.get_description(lang='nl') == "":
    wdPage.set_description("gen variant", "nl")
if len(variant_data["variant_aliases"]) > 0:
    for alias in variant_data["variant_aliases"]:
        synonyms.append(alias)
if len(synonyms) > 0:
    wdPage.set_aliases(aliases=synonyms, lang='en', append=True)

pprint(wdPage.get_wd_json_representation())
# wdPage.write(logincreds)


file.close()