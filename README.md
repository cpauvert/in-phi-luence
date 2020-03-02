# In-phi-luence: a network view of philosophers

A little project of data aggregation (R and Shiny) fetched from [Wikipedia (en)](https://en.wikipedia.org). The original idea stems from my partner wondering whether it was possible to visualise connections between philosophers and their schools of thoughts. I knew some influences were displayed on well-documented Wikipedia articles and so it began.

## How?

Wikipedia articles were first fetched if they had a Philosopher Infobox and if they belong to the category of Philosophy of science.
Mentions of influences in the infobox were mined, collected and gathered into a network.

The gathering of Wikipedia articles was made possible with the powerful search tool [PetScan](https://petscan.wmflabs.org).

## Why?


This app allows to visualise knowledge flow as a whole. It is also a way to promote the content of Wikipedia articles which are accessible once hovering on the nodes.
Scholars and students in philosophy could explore the network and even spot missing links that should be indicated in Wikipedia articles.

# Future avenues and suggestions

- [] Improve the visualisation (inward/outward links; list of influencers etc.)
- [] Improve the search in the list of philosophers to reduce the app load
- [] Expand the search to all philosophers
- [] Intersect the network with one built from articles in different languages
- [] Restrict the visualisation to philosophers only, not influences like novelists
- [] Add this list to proper issues in Github
