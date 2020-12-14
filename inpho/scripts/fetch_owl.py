# Fetch OWL files from InPho site 
# PAUVERT Charlie
# 2020-12-14
from requests import get
from bs4 import BeautifulSoup
import re as re
import os as os

# Fetch the page with the list of OWL files
urls = ["https://www.inphoproject.org", "owl/"]
response = get("/".join(urls))

# Scrap the page to get the list of links
html_soup = BeautifulSoup(response.text, 'html.parser')
# Find the links
lks=html_soup.find_all("a")

# Filtering the links by 
# 1. Getting the destination of the link 
owls=[ x['href'] for x in lks ]
# 2. Filter with a comprehension list to get only OWL files
owls=[ x for x in owls if re.match("^/owl.+\.owl$", x) is not None ]

# Create the output directory
if not os.path.exists(snakemake.output[0]):
    os.mkdir(snakemake.output[0])

# Download the OWL file if not already present
for owl in owls[-5:]:
    # Get into a list 
    # 1. the output directory and
    # 2. the OWL filename (~basename like function)
    owl_path=os.path.join(snakemake.output[0], owl.split("/")[-1])
    if not os.path.exists(owl_path):
        html=get(urls[0] + owl)
        with open(owl_path, "wb") as dwl:
            dwl.write(html.content)
    else:
        continue

