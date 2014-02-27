Reactome-FI
===========

R Reactome FI Cytoscape plugin wrapper.
This package currently implements a limited interface to ReactomeFI functionality including network construction, network clustering, and Pathway/GO term enrichment. It also implements some network plotting capabilities using either sna/ggplot2 or RCytoscape.

#### Dependencies

* plyr
* ggplot2
* sna
* RCurl
* XML

#### Get Started:

    network <- ReactomeFINetwork("2012", genes)
    annotate(network, "Pathway")
    plot(network)

    network <- ReactomeFINetwork("2012", genes, cluster=TRUE)
    annotateModules(network, "Pathway")
    plot(network)
