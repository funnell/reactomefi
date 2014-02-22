Reactome-FI
===========

R Reactome FI Cytoscape plugin wrapper


#### Get Started:

    network <- ReactomeFINetwork("2012", genes)
    annotate(network, "Pathway")
    plot(network)

    network <- ReactomeFINetwork("2012", genes, cluster=TRUE)
    annotateModules(network, "Pathway")
    plot(network)
