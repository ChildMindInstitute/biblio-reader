
""" Network analysis of authorship for ohbm2011
by Russ Poldrack, June 30, 2011

dependencies:
- numpy: 
- matplotlib: 
- networx: 
- community.py:

also requires graph data from http://raw.github.com/poldrack/ohbm2011/full_author_graph_anon.graphml

"""

## Copyright 2011, Russell Poldrack. All rights reserved.

## Redistribution and use in source and binary forms, with or without modification, are
## permitted provided that the following conditions are met:

##    1. Redistributions of source code must retain the above copyright notice, this list of
##       conditions and the following disclaimer.

##    2. Redistributions in binary form must reproduce the above copyright notice, this list
##       of conditions and the following disclaimer in the documentation and/or other materials
##       provided with the distribution.

## THIS SOFTWARE IS PROVIDED BY RUSSELL POLDRACK ``AS IS'' AND ANY EXPRESS OR IMPLIED
## WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
## FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RUSSELL POLDRACK OR
## CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
## ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
## NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import operator
import sys
import numpy as N
import matplotlib.pyplot as plt
import networkx as nx
import community


# load graph - available at http://raw.github.com/poldrack/ohbm2011/full_author_graph_anon.graphml

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

print 'Distribution of component sizes:'
for x in N.sort(cc_dict.keys())[::-1]:
    print '%d: %d components'%(x,len(cc_dict[x]))

isolates=nx.isolates(Ga)
print '%d isolate nodes'%len(isolates)
print ''

pageranks = nx.pagerank_numpy(Ga)
sorted_x = sorted(pageranks.iteritems(), key=operator.itemgetter(1))
for x in sorted_x[-100:]:
    print str(x[0]) + "\t" + str(x[1])

sys.exit(0)

print 'computing graph positions and saving graphs - can take a while'
pos=nx.spring_layout(Ga)
plt.clf()
nx.draw_networkx(Ga,pos,with_labels=False,node_color='k',node_size=1,width=0.25,edge_color='0.8')
plt.savefig('plots/full_network_citations.pdf',format='pdf',orientation='landscape',dpi=150,width=0.25)

# create a matched random graph and compute network statistics
print 'creating matched random graph and computing statistics - this could take a moment...'
rg=nx.fast_gnp_random_graph(Ga.number_of_nodes(),2.0*Ga.number_of_edges()/(Ga.number_of_nodes()*(Ga.number_of_nodes()-1)))
c_rg=nx.average_clustering(rg)
rg_cc=nx.connected_component_subgraphs(rg)[0]
rg_asp=nx.algorithms.shortest_paths.generic.average_shortest_path_length(rg_cc)

p_rg=community.best_partition(rg_cc)
m_rg=community.modularity(p_rg,rg_cc)

pageranks = nx.pagerank_numpy(rg)
sorted_x = sorted(pageranks.iteritems(), key=operator.itemgetter(1))
for x in sorted_x[-10:]:
    print x[0], x[1]

# compute network statistics
print 'Network statistics (values from matched random graph)'
c=nx.average_clustering(Ga[0])
print 'Average clustering: %f (%f)'%(c,c_rg)

asp=nx.algorithms.shortest_paths.generic.average_shortest_path_length(Ga)
print 'Average shortest path: %f (%f)'%(asp,rg_asp)

print 'plotting degree distribution'
degree_Ga=nx.degree_histogram(Ga)
degree_rg=nx.degree_histogram(rg)

log_degree=N.log10(degree_Ga[1:])
log_degree[N.where(log_degree<0)]=0
log_n=N.log10(range(1,len(degree_Ga)))
plt.clf()
plt.plot(log_n,log_degree,label='OHBM 2011')

log_degree_rg=N.log10(degree_rg[1:])
log_degree_rg[N.where(log_degree_rg<0)]=0
log_n_rg=N.log10(range(1,len(degree_rg)))
plt.plot(log_n_rg,log_degree_rg,label='random graph')
plt.xlabel('Degree (log10)')
plt.ylabel('Frequency (log10)')
plt.legend()
plt.title('Degree distribution')

plt.savefig('degree_distribution.pdf',format='pdf',orientation='landscape',dpi=150,width=0.25)

