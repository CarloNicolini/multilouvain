#include "DCSurpriseVertexPartition.h"

#ifdef DEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

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
#ifdef DEBUG
    cerr << "virtual double DCSurpriseVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
#endif
    size_t old_comm = this->membership(v);
    size_t nsize = this->graph->node_size(v);
#ifdef DEBUG
    cerr << "\t" << "nsize: " << nsize << endl;
#endif
    double diff = 0.0;
    if (new_comm != old_comm)
    {
        double normalise = (2.0 - this->graph->is_directed());
        double m = this->graph->total_weight();

#ifdef DEBUG
        cerr << "\t" << "Community: " << old_comm << " => " << new_comm << "." << endl;
        cerr << "\t" << "m: " << m << ", n2: " << n2 << "." << endl;
#endif
        // Before move
        double mc = this->total_weight_in_all_comms();
        double sumKc2 = 0.0;
        for (size_t c = 0; c < this->nb_communities(); c++)
        {
            double w_out = this->total_weight_from_comm(c);
            double w_in = this->total_weight_to_comm(c);
#ifdef DEBUG
            size_t csize = this->csize(c);
            cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
#endif
            sumKc2 += w_out*w_in;
        }

#ifdef DEBUG
        cerr << "\t" << "mc: " << mc << ", nc2: " << nc2 << "." << endl;
#endif
        double sw = this->graph->node_self_weight(v);
        double wtc = this->weight_to_comm(v, old_comm) - sw;
        double wfc = this->weight_from_comm(v, old_comm) - sw;
#ifdef DEBUG
        cerr << "\t"  << "wtc: " << wtc << ", wfc: " << wfc << ", sw: " << sw << "." << endl;
#endif
        double m_old = wtc/normalise + wfc/normalise + sw;
#ifdef DEBUG
        cerr << "\t" << "m_old: " << m_old << ", n_old: " << n_old << "." << endl;
#endif

        // To new comm
        size_t n_new = this->csize(new_comm);
        wtc = this->weight_to_comm(v, new_comm);
        wfc = this->weight_from_comm(v, new_comm);
        sw = this->graph->node_self_weight(v);
#ifdef DEBUG
        cerr << "\t"  << "wtc: " << wtc << ", wfc: " << wfc << ", sw: " << sw << "." << endl;
#endif
        double m_new = wtc/normalise + wfc/normalise + sw;
#ifdef DEBUG
        cerr << "\t" << "m_new: " << m_new << ", n_new: " << n_new << "." << endl;
#endif

        double q = mc/m;
        double s = (double)nc2/(double)n2;
        double q_new = (mc - m_old + m_new)/m;
#ifdef DEBUG
        cerr << "\t" << "mc - m_old + m_new=" << (mc - m_old + m_new) << endl;
#endif
        double s_new = (double)(nc2 + 2*nsize*(n_new - n_old + nsize)/normalise)/(double)n2;
#ifdef DEBUG
        cerr << "\t" << "nc2 + 2*nsize*(n_new - n_old + nsize)/normalise=" << nc2 + 2*nsize*(n_new - n_old + nsize)/normalise << endl;
#endif
#ifdef DEBUG
        cerr << "\t" << "q:\t" << q << ", s:\t"  << s << "." << endl;
        cerr << "\t" << "q_new:\t" << q_new << ", s_new:\t"  << s_new << "." << endl;
#endif
        diff = m*(KL(q_new, s_new) - KL(q, s));

#ifdef DEBUG
        cerr << "\t" << "diff: " << diff << "." << endl;
#endif
    }
#ifdef DEBUG
    cerr << "exit double DCSurpriseVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
#endif
    return diff;
}

double DCSurpriseVertexPartition::quality()
{
#ifdef DEBUG
    cerr << "double DCSurpriseVertexPartition::quality()" << endl;
#endif

    double sumKc2 = 0.0;
    for (size_t c = 0; c < this->nb_communities(); c++)
    {
        double w_out = this->total_weight_from_comm(c);
        double w_in = this->total_weight_to_comm(c);
#ifdef DEBUG
        size_t csize = this->csize(c);
        cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
#endif
        sumKc2 += w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
    }

    double mc = this->total_weight_in_all_comms();
    double m = this->graph->total_weight();

#ifdef DEBUG
    cerr << "\t" << "mc=" << mc << ", m=" << m << ", nc2=" << nc2 << ", n2=" << n2 << "." << endl;
#endif
    double q = mc/m;
    double s = sumKc2/(4*m*m);
#ifdef DEBUG
    cerr << "\t" << "q:\t" << q << ", s:\t"  << s << "." << endl;
#endif
    double S = m*KL(q,s);
#ifdef DEBUG
    cerr << "exit DCSurpriseVertexPartition::quality()" << endl;
    cerr << "return " << S << endl << endl;
#endif
    return S;
}
