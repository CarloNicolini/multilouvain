#include "DCSurpriseVertexPartition.h"

//#ifdef DEBUG
#include <iostream>
using std::cerr;
using std::endl;
//#endif

DCSurpriseVertexPartition::DCSurpriseVertexPartition(Graph* graph,
                                                     vector<size_t> membership) :
    MutableVertexPartition(graph,
                           membership)
{ }

DCSurpriseVertexPartition::DCSurpriseVertexPartition(Graph* graph) :
    MutableVertexPartition(graph)
{ }

DCSurpriseVertexPartition* DCSurpriseVertexPartition::create(Graph* graph)
{
    return new DCSurpriseVertexPartition(graph);
}

DCSurpriseVertexPartition::~DCSurpriseVertexPartition()
{ }

double DCSurpriseVertexPartition::diff_move(size_t v, size_t new_comm)
{
    size_t old_comm = this->membership(v);
    size_t nsize = this->graph->node_size(v);

    double diff = 0.0;
    if (new_comm != old_comm)
    {
        double normalise = (2.0 - this->graph->is_directed());
        double m = this->graph->total_weight();
        size_t n = this->graph->total_size();
        size_t n2 = 0;

        if (this->graph->correct_self_loops())
            n2 = n*n/normalise;
        else
            n2 = n*(n-1)/normalise;
        // Before move
        double mc = this->total_weight_in_all_comms();
        size_t nc2 = this->total_possible_edges_in_all_comms();
        // To old comm
        size_t n_old = this->csize(old_comm);
        double sw = this->graph->node_self_weight(v);
        double wtc = this->weight_to_comm(v, old_comm) - sw;
        double wfc = this->weight_from_comm(v, old_comm) - sw;
        double m_old = wtc/normalise + wfc/normalise + sw;

        // To new comm
        size_t n_new = this->csize(new_comm);
        wtc = this->weight_to_comm(v, new_comm);
        wfc = this->weight_from_comm(v, new_comm);
        sw = this->graph->node_self_weight(v);
        double m_new = wtc/normalise + wfc/normalise + sw;


        //***************************** CONFIGURATION MODEL PART ************
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
        double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());
        double diff_old = (w_to_old - k_out*K_in_old/total_weight) + \
          (w_from_old - k_in*K_out_old/total_weight);
        double diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) + \
          (w_from_new + self_weight - k_in*K_out_new/total_weight);
        diff = diff_new - diff_old;

        //********************** END CONFIGURATION MODEL PART **************
        double q = mc/m;
        double s = k_out*K_in_old/total_weight + k_in*K_out_old/total_weight;
        double q_new = (mc - m_old + m_new)/m;
        double s_new = k_out*K_in_new/total_weight + k_in*K_out_new/total_weight;
        //double s_new = (double)(nc2 + 2*nsize*(n_new - n_old + nsize)/normalise)/(double)n2;
        //cerr << "\t" << "nc2 + 2*nsize*(n_new - n_old + nsize)/normalise=" << nc2 + 2*nsize*(n_new - n_old + nsize)/normalise << endl;
        diff = m*(KL(q_new, s_new/(4*m)) - KL(q, s/(4*m)));
    }

    return diff;
}

double DCSurpriseVertexPartition::quality()
{
#ifdef DEBUG
    cerr << "double DCSurpriseVertexPartition::quality()" << endl;
#endif

    double mc = this->total_weight_in_all_comms();
    double m = this->graph->total_weight();

    double sumKc2 = 0.0;
    for (size_t c = 0; c < this->nb_communities(); c++)
    {
        double w_out = this->total_weight_from_comm(c);
        double w_in = this->total_weight_to_comm(c);
//#ifdef DEBUG
        //size_t csize = this->csize(c);
        //cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
//#endif
        //sumKc2 += w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
        sumKc2 += w_out*w_in/(4*m*m);
    }

    return m*KL(mc/m,sumKc2/m);

#ifdef DEBUG
    cerr << "\t" << "mc=" << mc << ", m=" << m << ", nc2=" << nc2 << ", n2=" << n2 << "." << endl;
#endif
    double q = mc/m;
    double s = sumKc2;
//#ifdef DEBUG
    cerr << "\t" << "q:\t" << q << ", s:\t"  << s << "." << endl;
//#endif
    double S = m*KL(q,s);
#ifdef DEBUG
    cerr << "exit DCSurpriseVertexPartition::quality()" << endl;
    cerr << "return " << S << endl << endl;
#endif
    return S;
}
