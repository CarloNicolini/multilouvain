#include "KLModularityVertexPartition.h"
#include <iostream>
using std::cerr;
using std::endl;


KLModularityVertexPartition::KLModularityVertexPartition(Graph* graph,
        vector<size_t> membership) :
    MutableVertexPartition(graph,
                           membership)
{ }

KLModularityVertexPartition::KLModularityVertexPartition(Graph* graph) :
    MutableVertexPartition(graph)
{ }

KLModularityVertexPartition* KLModularityVertexPartition::create(Graph* graph)
{
    return new KLModularityVertexPartition(graph);
}

KLModularityVertexPartition::~KLModularityVertexPartition()
{ }

double KLModularityVertexPartition::diff_move(size_t v, size_t new_comm)
{
    size_t old_comm = this->_membership[v];
    double diff = 0.0;
    if (new_comm != old_comm)
    {
        double normalise = (2.0 - this->graph->is_directed());
        double total_weight = this->graph->total_weight()*(normalise);

        double w_to_old = this->weight_to_comm(v, old_comm);
        double w_from_old = this->weight_from_comm(v, old_comm);
        double w_to_new = this->weight_to_comm(v, new_comm);
        double w_from_new = this->weight_from_comm(v, new_comm);
        double k_out = this->graph->strength(v, IGRAPH_OUT);
        double k_in = this->graph->strength(v, IGRAPH_IN);
        double self_weight = this->graph->node_self_weight(v);

        double K_out_old = this->total_weight_from_comm(old_comm);
        double K_in_old = this->total_weight_to_comm(old_comm);

        double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
        double K_in_new = this->total_weight_to_comm(new_comm) + k_in;


/*
        double diff_old = (w_to_old - k_out*K_in_old/total_weight) +
                          (w_from_old - k_in*K_out_old/total_weight);

        double diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) +
                          (w_from_new + self_weight - k_in*K_out_new/total_weight);
*/
        double mnew

        diff = diff_new - diff_old;
    }

    return diff;
}

double KLModularityVertexPartition::quality()
{
    double mc = this->total_weight_in_all_comms();
    double m = this->graph->total_weight();

    // COMPUTE MODULARITY (calcolo corretto ho verificato)
    double mod = 0.0;
    for (size_t c = 0; c < this->nb_communities(); c++)
    {
        double w = this->total_weight_in_comm(c);
        double w_out = this->total_weight_from_comm(c);
        double w_in = this->total_weight_to_comm(c);
        double configurationProb = w_out*w_in/pow((this->graph->is_directed() ? 1.0 : 2.0)*this->graph->total_weight(),2);
        cerr << "Comm=" << c << " w_in=" << w << " Kc^2/4m^2=" << configurationProb << endl;
        mod += KL( w/m, configurationProb );
    }
    cerr << "KLmod=" << mod << endl;
    return mod;
}
