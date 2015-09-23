#include "CPMVertexPartition.h"

CPMVertexPartition::CPMVertexPartition(Graph* graph,
                                       vector<size_t> membership, double resolution_parameter) :
    LinearResolutionParameterVertexPartition(graph,
            membership, resolution_parameter)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph,
                                       vector<size_t> membership) :
    LinearResolutionParameterVertexPartition(graph,
            membership)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph,
                                       double resolution_parameter) :
    LinearResolutionParameterVertexPartition(graph, resolution_parameter)
{ }

CPMVertexPartition::CPMVertexPartition(Graph* graph) :
    LinearResolutionParameterVertexPartition(graph)
{ }

CPMVertexPartition::~CPMVertexPartition()
{ }

CPMVertexPartition* CPMVertexPartition::create(Graph* graph)
{
    return new CPMVertexPartition(graph, this->resolution_parameter);
}

/********************************************************************************
  RBER implementation of a vertex partition
  (which includes a resolution parameter).
 ********************************************************************************/
double CPMVertexPartition::diff_move(size_t v, size_t new_comm)
{
#ifdef DEBUG
    cerr << "double CPMVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "Using resolution parameter: " << this->resolution_parameter << "." << endl;
#endif
    size_t old_comm = this->membership(v);
    double diff = 0.0;
    if (new_comm != old_comm)
    {
        double w_to_old = this->weight_to_comm(v, old_comm);
#ifdef DEBUG
        cerr << "\t" << "w_to_old: " << w_to_old << endl;
#endif
        double w_to_new = this->weight_to_comm(v, new_comm);
#ifdef DEBUG
        cerr << "\t" << "w_to_new: " << w_to_new << endl;
#endif
        double w_from_old = this->weight_from_comm(v, old_comm);
#ifdef DEBUG
        cerr << "\t" << "w_from_old: " << w_from_old << endl;
#endif
        double w_from_new = this->weight_from_comm(v, new_comm);
#ifdef DEBUG
        cerr << "\t" << "w_from_new: " << w_from_new << endl;
#endif
        size_t nsize = this->graph->node_size(v);
#ifdef DEBUG
        cerr << "\t" << "nsize: " << nsize << endl;
#endif
        size_t csize_old = this->csize(old_comm);
#ifdef DEBUG
        cerr << "\t" << "csize_old: " << csize_old << endl;
#endif
        size_t csize_new = this->csize(new_comm);
#ifdef DEBUG
        cerr << "\t" << "csize_new: " << csize_new << endl;
#endif
        double self_weight = this->graph->node_self_weight(v);
#ifdef DEBUG
        cerr << "\t" << "self_weight: " << self_weight << endl;
        cerr << "\t" << "density: " << this->graph->density() << endl;
#endif
        double possible_edge_difference_old = 0.0;
        if (this->graph->correct_self_loops())
            possible_edge_difference_old = nsize*(2.0*csize_old - nsize);
        else
            possible_edge_difference_old = nsize*(2.0*csize_old - nsize - 1.0);
#ifdef DEBUG
        cerr << "\t" << "possible_edge_difference_old: " << possible_edge_difference_old << endl;
#endif
        double diff_old = w_to_old + w_from_old -
                          self_weight - this->resolution_parameter*possible_edge_difference_old;
#ifdef DEBUG
        cerr << "\t" << "diff_old: " << diff_old << endl;
#endif
        double possible_edge_difference_new = 0.0;
        if (this->graph->correct_self_loops())
            possible_edge_difference_new = nsize*(2.0*csize_new + nsize);
        else
            possible_edge_difference_new = nsize*(2.0*csize_new + nsize - 1.0);
#ifdef DEBUG
        cerr << "\t" << "possible_edge_difference_new: " << possible_edge_difference_new << endl;
#endif
        double diff_new = w_to_new + w_from_new + self_weight -
                          this->resolution_parameter*possible_edge_difference_new;
#ifdef DEBUG
        cerr << "\t" << "diff_new: " << diff_new << endl;
#endif
        diff = diff_new - diff_old;
#ifdef DEBUG
        cerr << "\t" << "diff: " << diff << endl;;
#endif
    }
#ifdef DEBUG
    cerr << "exit CPMVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
#endif
    return diff;
}

double CPMVertexPartition::quality()
{
#ifdef DEBUG
    cerr << "double CPMVertexPartition::quality()" << endl;
#endif
    double mod = 0.0;
    for (size_t c = 0; c < this->nb_communities(); c++)
    {
        size_t csize = this->csize(c);
        double w = this->total_weight_in_comm(c);
        size_t comm_possible_edges;
        if (this->graph->correct_self_loops())
            comm_possible_edges = csize*csize;
        else
            comm_possible_edges = csize*(csize - 1);
        if (!this->graph->is_directed())
            comm_possible_edges /= 2;
#ifdef DEBUG
        cerr << "\t" << "Comm: " << c << ", w_c=" << w << ", n_c=" << csize << ", comm_possible_edges=" << comm_possible_edges << ", p=" << this->graph->density() << "." << endl;
#endif
        mod += w - this->resolution_parameter*comm_possible_edges;
    }
#ifdef DEBUG
    cerr << "exit double CPMVertexPartition::quality()" << endl;
    cerr << "return " << mod << endl << endl;
#endif
    return (2.0 - this->graph->is_directed())*mod;
}
