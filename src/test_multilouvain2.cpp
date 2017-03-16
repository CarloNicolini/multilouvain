#include <iostream>
#include <fstream>
#include "ModularityVertexPartition.h"
#include "igraph_utils.h"
#include "GraphHelper.h"
#include "Optimiser.h"

using namespace std;


int main(int argc, char *argv[])
{
    srand(time(0));
    // Create the graph from the adjacency matrix
    // Fill the adjacency matrix of the graph
    Eigen::MatrixXd A = read_adj_matrix<Eigen::MatrixXd>(std::string(argv[1]));
    cout << A << endl;
    
    // Create the Graph helper object specifying edge weights too
    Graph *G = init(A.data(),A.rows(),A.cols());

    // Create the partition function
    MutableVertexPartition *partition = new ModularityVertexPartition(G);

    Optimiser *opt = new Optimiser;
    double qual = opt->optimize_partition(partition);

    cout << "Q = " <<  qual << endl;
    for (size_t i=0; i<partition->membership().size(); ++i)
        cout << partition->membership(i) << " ";
    cout << endl;

    delete opt;
    delete partition;
    G->dispose();
    delete G;
    return 0;
}
