Reactome-FI
===========

R Reactome FI Cytoscape plugin wrapper


#### Get Started:

    network <- ReactomeFINetwork("2012")
    network <- build(network, genes)
    annotate(network, "Pathway")

    network <- cluster(network)
    annotateModules(network, "Pathway")
