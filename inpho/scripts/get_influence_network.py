# Get influence network from OWL file
# PAUVERT Charlie
# 2020-12-14
from owlready2 import *

# Get the aggregated OWL file
onto = get_ontology(snakemake.input[0]).load()

# Search the philosophers with the influenced properties
influencers = onto.search(has_influenced = "*")
with open(snakemake.output[0],"w") as network:
    for phil in influencers[1:4]:
        for influencee in phil.has_influenced:
            print(phil.get_name(), "\t", influencee.get_name(),
                    file=network)

