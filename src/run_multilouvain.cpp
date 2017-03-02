#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include "ModularityVertexPartition.h"
#include "GraphHelper.h"
#include "Optimiser.h"
#include "igraph_utils.h"

using namespace std;

void exit_with_help()
{
    std::printf(
                "Usage: run_multilouvain graph_file [options]\n"
                "Version 0.1 19 January 2017\n"
                "graph_file the file containing the graph. Accepted formats are pajek, graph_ml, adjacency matrix or"
                "\nedges list (the ncol format), additionally with a third column with edge weights"
                "options:\n"
                "-q [quality]\n"
                "   0 Asymptotical Surprise\n"


                "-c [consider_comms]:"
                "   0 Agglomerative Optimizer\n"
                "   1 Random\n"
                "   2 Simulated Annealing\n"
                "-V [report_level] ERROR=0, WARNING=1, INFO=2, DEBUG=3, DEBUG1=4, DEBUG2=5, DEBUG3=6, DEBUG4=7\n"
                "-S [seed] specify the random seed, default time(0)\n"
                "-b [bool] wheter to start with initial random cluster or every node in its community\n"
                "-r [repetitions], number of repetitions of PACO, default=1\n"
                "-p [print solution]\n"
                "\n"
                );
    exit(1);
}

enum error_type
{
    NO_ERROR = 0,
    ERROR_TOO_MANY_OUTPUT_ARGS = 1,
    ERROR_NOT_ENOUGH_ARGS = 2,
    ERROR_ARG_VALUE = 3,
    ERROR_ARG_TYPE = 4,
    ERROR_MATRIX = 5,
    ERROR_ARG_EMPTY=6,
    ERROR_ARG_UNKNOWN=7
};

//static const char *error_strings[] =
//{
//    "",
//    "Too many output arguments.",
//    "Not enough input arguments.",
//    "Non valid argument value.",
//    "Non valid argument type.",
//    "Non valid input adjacency matrix. PACO accepts symmetric real dense-type (n x n) matrices or sparse edges-list representation \
//    [num_edges x 3] array of edges list with edge endpoints and weight.",
//    "Expected some argument value but empty found.",
//    "Unkwown argument."
//};


void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}


int main(int argc, char *argv[])
{

    vector<double> edges_weights(0);
    std::vector<double> edges_list(0);
    // The container structure igraph_t

    std::vector<double> data = read_adj_matrix(std::string(argv[1]));

    int N = sqrt(data.size());
    Eigen::MatrixXd EW = Eigen::Map<Eigen::MatrixXd>(data.data(),N,N);
    if (EW.trace() > 0)
    {
        throw std::logic_error("Adjacency matrix has self loops, only simple graphs allowed.");
    }
    for (int i=0; i<N; ++i)
    {
        for (int j=i+1; j<N; j++)
        {
            double w = std::max(EW.coeffRef(i,j),EW.coeffRef(j,i));
            if (w>0)
            {
                edges_list.push_back(i);
                edges_list.push_back(j);
                edges_weights.push_back(w);
            }
            if (w<0)
            {
                throw std::logic_error("Negative edge weight found. Only positive weights supported.");
            }
        }
    }

    if (edges_list.empty())
        throw std::logic_error("Empty graph provided.");

    // Create the Graph object from the igraph data structure
    // Fill the edges into the igraph IG
    igraph_vector_t igedges_list;
    igraph_vector_view(&igedges_list, edges_list.data(), edges_list.size());
    igraph_t IG;
    IGRAPH_TRY(igraph_create(&IG, &igedges_list, 0, 0));
    // cout << "=========" << endl;
    // cerr << "G=(V,E)=" << igraph_vcount(&IG) << " " << igraph_ecount(&IG) << endl;
    //Graph *G  = new Graph(&IG,edges_weights);
    //cerr << G->ecount() << endl;
    //delete G;

    igraph_vector_destroy(&igedges_list);
    return 0;
}

