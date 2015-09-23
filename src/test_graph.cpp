#include <iostream>
#include "Optimiser.h"
#include "MutableVertexPartition.h"
#include "SurpriseVertexPartition.h"
#include "ModularityVertexPartition.h"
#include "GraphHelper.h"

#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
    srand(time(0));
    igraph_t graph;

    FILE *f;
    f = std::fopen("../data/karate.edge","r");
    igraph_read_graph_edgelist(&graph,f,34,0);

    Graph *G = new Graph(&graph);
    Optimiser* opt = new Optimiser();
    opt->consider_comms = Optimiser::ALL_NEIGH_COMMS;
    opt->random_order = true;

    MutableVertexPartition* partition = new SurpriseVertexPartition(G);
    opt->optimize_partition(partition);
    delete opt;
    delete partition;
    delete G;
    igraph_destroy(&graph);

    return 0;
}
