# In-phi-luence: a network view of philosophers

A little project of data aggregation (R and Shiny and a pinch of Python) fetched from [Wikipedia (en)](https://en.wikipedia.org) and more recently from [The Internet Philosophy Ontology Project (InPho)](https://www.inphoproject.org/) . The original idea stems from my partner wondering whether it was possible to visualise connections between philosophers and their schools of thoughts.

## Why?

This app allows to visualise knowledge flow as a whole. It is also a way to promote the content of knowledge sources that are accessible once hovering on the nodes.
Scholars and students in philosophy could explore the network and even spot missing links that should be indicated in Wikipedia articles.

## How?

Two sources of philosophers influences were used:

* The Free Encyclopedia Wikipedia (en)
* The Internet Philosophy Ontology Project (InPho)

### Wikipedia (en)

Wikipedia articles were first fetched if they had a Philosopher Infobox and if they belong to the category of Philosophy of science.
Mentions of influences in the infobox were mined, collected and gathered into a network.
Additional philosophers, novelist and school of thoughts were also added after following the influences linked in the listed Wikipedia pages. The currently displayed network was generated on 2020-02-14.
The gathering of Wikipedia articles was made possible with the powerful search tool [PetScan](https://petscan.wmflabs.org).

R code was then used to fetch and clean data and is available at https://github.com/cpauvert/in-phi-luence/in-phi-luence.R

### InPho

This [amazing scholarly resource](https://www.inphoproject.org) compiles ontologies on philosophers and which are then made accessible through API or OWL files. Monthly archives of the InPho ontologies were fetched.
Automatic mining of the ontologies extracted relevant properties (such as `has_influenced` or `was_influenced`), and all results were gather into a network.
The currently displayed network was generated on 2020-12-16.
R and Python code (as well as a Snakemake workflow) were used to fetch and clean data and are available at https://github.com/cpauvert/in-phi-luence/inpho

# Future avenues and suggestions

Future avenues include an update of the sources used and a comparison between the two generated networks.
The latter could help precise missing influences in Wikipedia articles for instance.
Suggestions are welcomed in the issues section where some of them are already listed.
