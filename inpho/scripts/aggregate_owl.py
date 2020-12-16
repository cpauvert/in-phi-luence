# Aggregate OWL files into one
# PAUVERT Charlie
# 2020-12-14
from owlready2 import *

# Fetch the list of OWL files
owls = list(snakemake.input) # Necessary copy to avoid AttributeError from Snakemake
# Extract the very first OWL from the list to initiate the object
owls.sort()
owls.reverse()
first_owl = owls.pop()
owls.reverse()

# Load the first OWL
onto = get_ontology(first_owl).load()

# Update the initial OWL object with the others
for update in owls:
    try:
        owl_update = get_ontology(update).load()
        onto.imported_ontologies.append(owl_update)
    except:
        print("Skipped ", update)
        continue

# Write the aggregated file
onto.save(snakemake.output[0])
