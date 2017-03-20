import operator
import sys
import numpy as N
import networkx as nx
import community
import io, unicodedata

# load graph - available at http://raw.github.com/poldrack/ohbm2011/full_author_graph_anon.graphml																	   
def getPageRanks(filename):
	Ga=nx.read_graphml(filename)
	pageranks = nx.pagerank_numpy(Ga)
	return pageranks

def getRandomPageRanks(filename):
	Ga=nx.read_graphml(sys.argv[1])

	# create a copy of the graph and extract giant component
	# get component size distribution
	cc=nx.connected_components(Ga)
	cc_dict={}
	for x in range(0,len(cc)):
		try:
			cc_dict[len(cc[x])].append(x)
		except KeyError:
			cc_dict[len(cc[x])]=[]
			cc_dict[len(cc[x])].append(x)

	isolates=nx.isolates(Ga)

	rg=nx.fast_gnp_random_graph(Ga.number_of_nodes(),2.0*Ga.number_of_edges()/(Ga.number_of_nodes()*(Ga.number_of_nodes()-1)))
	c_rg=nx.average_clustering(rg)
	rg_cc=nx.connected_component_subgraphs(rg)[0]
	rg_asp=nx.algorithms.shortest_paths.generic.average_shortest_path_length(rg_cc)

	p_rg=community.best_partition(rg_cc)
	m_rg=community.modularity(p_rg,rg_cc)

	pageranks = nx.pagerank_numpy(rg)
	return pageranks

if __name__ == "__main__":
	filename = sys.argv[1]
	pageranks = getPageRanks(filename)
	pageranks = sorted(pageranks.iteritems(), key=operator.itemgetter(1))
	pageranks.reverse()
	
	print "\t".join(["Title", "PageRank"])
	for x in pageranks[:100]:
		print "\t".join([str(x[0]), str(x[1])])
	