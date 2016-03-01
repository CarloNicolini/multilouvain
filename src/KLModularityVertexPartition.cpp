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
    /*
    cerr << "double KLModularityVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    size_t old_comm = this->_membership[v];
    double diff = 0.0;
    double m = this->graph->total_weight();

    if (new_comm != old_comm)
    {
        double w_to_old = this->weight_to_comm(v, old_comm);
        cerr << "\t" << "w_to_old: " << w_to_old << endl;

        double w_to_new = this->weight_to_comm(v, new_comm);
        cerr << "\t" << "w_to_new: " << w_to_new << endl;

        double k = this->graph->strength(v, IGRAPH_OUT);
        cerr << "\t" << "k: " << k << endl;

        double self_weight = this->graph->node_self_weight(v);
        cerr << "\t" << "self_weight: " << self_weight << endl;

        double K_old = this->total_weight_from_comm(old_comm);
        cerr << "\t" << "K_old: " << K_old << endl;

        double K_new = this->total_weight_from_comm(new_comm) + k;
        cerr << "\t" << "K_new: " << K_new << endl;
        cerr << "\t" << "total_weight: " << m << endl;

        double diff_old = (w_to_old - k*K_old/(2*m)) + (w_to_old - k*K_old/(2*m));

        double diff_new = (w_to_new + self_weight - k*K_new/(2*m)) + (w_to_new + self_weight - k*K_new/(2*m));

        double Qold = 2*(this->total_weight_in_comm(old_comm)/m - pow(total_weight_from_comm(old_comm)/(2*m),2));
        double Qnew = 2*((this->total_weight_in_comm(new_comm)+)/m - pow((total_weight_from_comm(new_comm)-)/(2*m),2));

        cerr << Qnew << " " << Qold << endl;
        //cerr << "Q0=" << this->total_weight_in_comm(old_comm) << " " << total_weight_from_comm(old_comm) << endl;
        //cerr << "Q1=" << this->total_weight_in_comm(new_comm) << " " << total_weight_from_comm(new_comm) << endl;

        diff = (diff_new - diff_old);
    }

    cerr << "exit double KLModularityVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
    return diff;
    */
#ifdef DEBUG
    cerr << "virtual double SignificanceVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
#endif
    size_t old_comm = this->membership(v);
    double m = graph->total_weight();
    double diff = 0.0;
    if (graph->is_directed())
        throw std::runtime_error("Not working with directed graphs");
    double normalise = (2.0 - this->graph->is_directed());
    double k = this->graph->strength(v, IGRAPH_OUT);

    if (new_comm != old_comm)
    {
        //Old comm
        double m_old = this->total_weight_in_comm(old_comm); // number of edges in old community
        double q_old = m_old/m; // fraction of edges in old community
        double K_old = this->total_weight_from_comm(old_comm); // number of edge stubs in old community
        // da mettere total_weight_from * total_weight_to se si vuole direzionato
        double p_old = (K_old*K_old)/(4*m*m);

        // Old comm after move
        double sw = this->graph->node_self_weight(v);
        // Be careful to exclude the self weight here, because this is include in the weight_to_comm function.
        double wtc = this->weight_to_comm(v, old_comm) - sw;
        double wfc = this->weight_from_comm(v, old_comm) - sw;

        double m_oldx = m_old - wtc/normalise - wfc/normalise - sw;
        double q_oldx = m_oldx/m;
        double K_oldx = K_old - k;
        double p_oldx = (K_oldx*K_oldx)/(4*m*m);

        // New comm
        double m_new = this->total_weight_in_comm(new_comm);
        double q_new = m_new/m;
        double K_new = this->total_weight_from_comm(new_comm);
        double p_new = (K_new*K_new)/(4*m*m);

        // New comm after move
        wtc = this->weight_to_comm(v, new_comm);
        wfc = this->weight_from_comm(v, new_comm);
        sw = this->graph->node_self_weight(v);

        double m_newx = m_new + wtc/normalise + wfc/normalise + sw;
        double q_newx = m_newx/m;
        double K_newx = K_new + k;
        double p_newx = (K_newx*K_newx)/(4*m*m);

        double d1,d2,d3,d4;
        d1 = KL(q_old, p_old);
        d2 = KL(q_oldx, p_oldx);
        d3 = KL(q_new, p_new);
        d4 = KL(q_newx, p_newx);

        // Calculate actual diff
        diff = - KL(q_old, p_old)
               + KL(q_oldx, p_oldx)
               - KL(q_new, p_new)
               + KL(q_newx, p_newx);
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
        mod +=  KL(w/m , configurationProb);
    }
    cerr << "mod=" << mod << endl;
    return mod;
}
