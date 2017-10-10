__author__ = 'andra'
import requests
import pprint

r = requests.get('https://civic.genome.wustl.edu/api/variants?count=10000')
variant_data = r.json()
pprint.pprint(variant_data)
for record in variant_data['records']:
    print(record['id'])