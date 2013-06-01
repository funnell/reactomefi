Reactome-FI
===========

R Reactome FI Cytoscape plugin wrapper


#### Get Started:

    network <- ReactomeFINetwork(genes, "2012")
    annotate(network, "Pathway")

    network <- cluster(network)
    annotateModules(network, "Pathway")
