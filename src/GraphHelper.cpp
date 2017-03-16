#include <Eigen/Core>
#include "GraphHelper.h"

#ifdef MATLAB_SUPPORT
#include "mexInterrupt.h"
#endif

#ifdef DEBUG
using std::cerr;
using std::endl;
#endif

vector<size_t> range(size_t n)
{
    vector<size_t> range_vec(n);
    for (size_t i = 0; i < n; i++)
        range_vec[i] = i;
    return range_vec;
}

/****************************************************************************
  The binary Kullback-Leibler divergence.
****************************************************************************/
double KL(double q, double p)
{
    double KL = 0.0;
    if (q > 0.0 && p > 0.0)
        KL += q * log(q / p);
    if (q < 1.0 && p < 1.0)
        KL += (1.0 - q) * log((1.0 - q) / (1.0 - p));
    return KL;
}

Graph::Graph(igraph_t* graph, const vector<double> &edge_weights, const vector<size_t> &node_sizes, const vector<double> &node_self_weights, int correct_self_loops)
{
    this->_graph = graph;
    this->_remove_graph = false;

    if (edge_weights.size() != this->ecount())
        throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
    this->_edge_weights = edge_weights;
    this->_is_weighted = true;

    if (node_sizes.size() != this->vcount())
        throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
    this->_node_sizes = node_sizes;

    if (node_self_weights.size() != this->vcount())
        throw Exception("Node self weights vector inconsistent length with the vertex count of the graph.");
    this->_node_self_weights = node_self_weights;

    this->_correct_self_loops = correct_self_loops;
    this->init_admin();
}

Graph::Graph(igraph_t* graph, const vector<double> &edge_weights, const vector<size_t> &node_sizes, const vector<double> &node_self_weights)
{
    this->_graph = graph;
    this->_remove_graph = false;

    if (edge_weights.size() != this->ecount())
        throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
    this->_edge_weights = edge_weights;
    this->_is_weighted = true;

    if (node_sizes.size() != this->vcount())
        throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
    this->_node_sizes = node_sizes;

    this->_node_self_weights = node_self_weights;
    this->init_admin();
}

Graph::Graph(igraph_t* graph, const vector<double> &edge_weights, const vector<size_t> &node_sizes, int correct_self_loops)
{
    this->_graph = graph;
    this->_remove_graph = false;

    if (edge_weights.size() != this->ecount())
        throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
    this->_edge_weights = edge_weights;
    this->_is_weighted = true;

    if (node_sizes.size() != this->vcount())
        throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
    this->_node_sizes = node_sizes;

    this->_correct_self_loops = correct_self_loops;
    this->init_admin();
    this->set_self_weights();
}

Graph::Graph(igraph_t* graph, const vector<double> &edge_weights, const vector<size_t> &node_sizes)
{
    this->_graph = graph;
    this->_remove_graph = false;
    if (edge_weights.size() != this->ecount())
        throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
    this->_edge_weights = edge_weights;
    this->_is_weighted = true;

    if (node_sizes.size() != this->vcount())
        throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
    this->_node_sizes = node_sizes;

    this->init_admin();
    this->set_self_weights();
}

Graph::Graph(igraph_t* graph, const vector<double> &edge_weights, int correct_self_loops)
{
    this->_graph = graph;
    this->_remove_graph = false;
    this->_correct_self_loops = correct_self_loops;
    if (edge_weights.size() != this->ecount())
        throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
    this->_edge_weights = edge_weights;
    this->_is_weighted = true;
    this->set_default_node_size();
    this->init_admin();
    this->set_self_weights();
}

Graph::Graph(igraph_t* graph, const vector<double> &edge_weights)
{
    this->_graph = graph;
    this->_remove_graph = false;
    if (edge_weights.size() != this->ecount())
        throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
    this->_edge_weights = edge_weights;
    this->_is_weighted = true;
    this->set_default_node_size();
    this->init_admin();
    this->set_self_weights();
    /*
        igraph_adjlist_t adjlist;
        igraph_neimode_t mode = IGRAPH_TOTAL;
        igraph_adjlist_init(_graph,&adjlist,mode);
        igraph_adjlist_print(&adjlist);
        igraph_adjlist_destroy(&adjlist);
    */
}

Graph::Graph(igraph_t* graph, int correct_self_loops)
{
    this->_graph = graph;
    this->_remove_graph = false;
    this->_correct_self_loops = correct_self_loops;
    this->set_defaults();
    this->_is_weighted = false;
    this->init_admin();
    this->set_self_weights();
}

Graph::Graph(igraph_t* graph)
{
    this->_graph = graph;
    this->_remove_graph = false;
    this->set_defaults();
    this->_is_weighted = false;
    this->init_admin();
    this->set_self_weights();
}

Graph::Graph()
{
    this->_graph = new igraph_t();
    this->_remove_graph = true;
    this->set_defaults();
    this->_is_weighted = false;
    this->init_admin();
    this->set_self_weights();
}

Graph::~Graph()
{
    if (this->_remove_graph)
    {
        igraph_destroy(this->_graph);
        delete this->_graph;
    }
}

void Graph::set_defaults()
{
    this->set_default_edge_weight();
    this->set_default_node_size();
}

void Graph::set_default_edge_weight()
{
    size_t m = this->ecount();

    // Set default edge weight of 1.0
    this->_edge_weights.clear();
    this->_edge_weights.resize(m);
    fill(this->_edge_weights.begin(), this->_edge_weights.end(), 1.0);
    this->_is_weighted = false;
}

void Graph::set_default_node_size()
{
    size_t n = this->vcount();

    // Set default node size of 1
    this->_node_sizes.clear();
    this->_node_sizes.resize(n);
    fill(this->_node_sizes.begin(), this->_node_sizes.end(), 1);
}

void Graph::set_self_weights()
{
    size_t n = this->vcount();

    // Set default self_weights of the total weight of any possible self-loops
    this->_node_self_weights.clear();
    this->_node_self_weights.resize(n);
    for (size_t v = 0; v < n; v++)
    {
#ifdef DEBUG
        cerr << "\t" << "Size node " << v << ": " << this->node_size(v) << endl;
#endif
        double self_weight = 0.0;
        // There should be only one self loop
        igraph_integer_t eid;
        // Get edge id for self loop
        igraph_get_eid(this->_graph, &eid, v, v, this->is_directed(), false);
        if (eid >= 0)
            self_weight = this->edge_weight(eid);

        this->_node_self_weights[v] = self_weight;
#ifdef DEBUG
        cerr << "\t" << "Self weight node " << v << ": " << self_weight << endl;
#endif
    }
}

void Graph::init_admin()
{

    size_t m = this->ecount();

    // Determine total weight in the graph.
    this->_total_weight = 0.0;
    for (size_t e = 0; e < m; e++)
        this->_total_weight += this->edge_weight(e);

    // Make sure to multiply by 2 for undirected graphs
    //if (!this->is_directed())
    //  this->_total_weight *= 2.0;

    size_t n = this->vcount();

    this->_total_size = 0;
    for (size_t v = 0; v < n; v++)
        this->_total_size += this->node_size(v);

    igraph_vector_t weights;
    igraph_vector_t res;

    // Strength IN
    igraph_vector_init(&res, n);
    // Copy weights to an igraph_vector_t
    igraph_vector_init_copy(&weights, &this->_edge_weights[0], this->ecount());
    // Calculcate strength
    igraph_strength(this->_graph, &res, igraph_vss_all(), IGRAPH_IN, true, &weights);
    igraph_vector_destroy(&weights);

    // Assign to strength vector
    this->_strength_in.clear();
    this->_strength_in.resize(n);
    for (size_t v = 0; v < n; v++)
        this->_strength_in[v] = VECTOR(res)[v];
    igraph_vector_destroy(&res);

    // Strength OUT
    igraph_vector_init(&res, n);
    // Copy weights to an igraph_vector_t
    igraph_vector_init_copy(&weights, &this->_edge_weights[0], this->ecount());
    // Calculcate strength
    igraph_strength(this->_graph, &res, igraph_vss_all(), IGRAPH_OUT, true, &weights);
    igraph_vector_destroy(&weights);

    // Assign to strength vector
    this->_strength_out.clear();
    this->_strength_out.resize(n);
    for (size_t v = 0; v < n; v++)
        this->_strength_out[v] = VECTOR(res)[v];
    igraph_vector_destroy(&res);

    // Degree IN
    igraph_vector_init(&res, n);
    igraph_degree(this->_graph, &res, igraph_vss_all(), IGRAPH_IN, true);
    this->_degree_in.clear();
    this->_degree_in.resize(n);
    for (size_t v = 0; v < n; v++)
        this->_degree_in[v] = VECTOR(res)[v];
    igraph_vector_destroy(&res);

    // Degree OUT
    igraph_vector_init(&res, n);
    igraph_degree(this->_graph, &res, igraph_vss_all(), IGRAPH_OUT, true);
    this->_degree_out.clear();
    this->_degree_out.resize(n);
    for (size_t v = 0; v < n; v++)
        this->_degree_out[v] = VECTOR(res)[v];
    igraph_vector_destroy(&res);

    // Degree ALL
    igraph_vector_init(&res, n);
    igraph_degree(this->_graph, &res, igraph_vss_all(), IGRAPH_ALL, true);
    this->_degree_all.clear();
    this->_degree_all.resize(n);
    for (size_t v = 0; v < n; v++)
        this->_degree_all[v] = VECTOR(res)[v];
    igraph_vector_destroy(&res);

    // Calculate density;
    double w = this->total_weight();
    size_t n_size = this->total_size();

    // For now we default to not correcting self loops.
    this->_correct_self_loops = false;

    double normalise = 0.0;
    if (this->_correct_self_loops)
        normalise = n_size * n_size;
    else
        normalise = n_size * (n_size - 1);

    if (this->is_directed())
        this->_density = w / normalise;
    else
        this->_density = 2 * w / normalise;
}

double Graph::weight_tofrom_community(size_t v, size_t comm, vector<size_t>* membership, igraph_neimode_t mode)
{
    // Weight between vertex and community
#ifdef DEBUG
    cerr << "double Graph::weight_tofrom_vertex_set(" << v << ", " << comm << ", " << mode << ")." << endl;
#endif
    double total_w = 0.0;
    size_t degree = this->degree(v, mode);
    igraph_vector_t incident_edges, neighbours;
    igraph_vector_init(&incident_edges, degree);
    igraph_vector_init(&neighbours, degree);
    igraph_incident(this->_graph, &incident_edges, v, mode);
    igraph_neighbors(this->_graph, &neighbours, v, mode);
    for (size_t i = 0; i < degree; i++)
    {
        size_t u = VECTOR(neighbours)[i];

        // If it is an edge to the requested community
#ifdef DEBUG
        size_t u_comm = (*membership)[u];
#endif
        if ((*membership)[u] == comm)
        {
#ifdef DEBUG
            //cerr << "\t" << "Sum edge (" << v << "-" << u << "), Comm (" << comm << "-" << u_comm << ") weight: " << w << "." << endl;
#endif
            size_t e = VECTOR(incident_edges)[i];
            // Get the weight of the edge
            double w = this->_edge_weights[e];
            // Self loops appear twice here if the graph is undirected, so divide by 2.0 in that case.
            if (u == v && !this->is_directed())
                w /= 2.0;

            total_w += w;
        }
#ifdef DEBUG
        else
        {
            cerr << "\t" << "Ignore edge (" << v << "-" << u << "), Comm (" << comm << ") weight: " << this->_edge_weights[VECTOR(incident_edges)[i]] << "." << endl;
        }
#endif
    }
    igraph_vector_destroy(&incident_edges);
    igraph_vector_destroy(&neighbours);
#ifdef DEBUG
    cerr << "exit Graph::weight_tofrom_vertex_set(" << v << ", " << comm << ", " << mode << ")." << endl;
#endif
    return total_w;
}

vector< pair<size_t, size_t> >*
Graph::get_neighbour_edges(size_t v, igraph_neimode_t mode)
{
    size_t degree = this->degree(v, mode);
    vector< pair<size_t, size_t> >* neigh_edges
        = new vector< pair<size_t, size_t> >(degree);

    igraph_vector_t incident_edges, neighbours;
    igraph_vector_init(&incident_edges, degree);
    igraph_vector_init(&neighbours, degree);
    igraph_incident(this->_graph, &incident_edges, v, mode);
    igraph_neighbors(this->_graph, &neighbours, v, mode);
    for (size_t i = 0; i < degree; i++)
    {
        size_t e = VECTOR(incident_edges)[i];
        size_t u = VECTOR(neighbours)[i];
        (*neigh_edges)[i] = make_pair(u, e);
    }
    igraph_vector_destroy(&incident_edges);
    igraph_vector_destroy(&neighbours);
    return neigh_edges;
}

vector< size_t >*
Graph::get_neighbours(size_t v, igraph_neimode_t mode)
{
    size_t degree = this->degree(v, mode);

    igraph_vector_t neighbours;
    igraph_vector_init(&neighbours, degree);
    igraph_neighbors(this->_graph, &neighbours, v, mode);
    vector< size_t >* neighs = new vector< size_t >(
        igraph_vector_e_ptr(&neighbours, 0),
        igraph_vector_e_ptr(&neighbours, degree));
    igraph_vector_destroy(&neighbours);
    return neighs;
}

/********************************************************************************
 * This should return a random neighbour in O(1)
 ********************************************************************************/
size_t Graph::get_random_neighbour(size_t v, igraph_neimode_t mode)
{
    size_t node = v;
    size_t rand_neigh = -1;

    if (this->degree(v, mode) <= 0)
        throw Exception("Cannot select a random neighbour for an isolated node.");

    if (igraph_is_directed(this->_graph) && mode != IGRAPH_ALL)
    {
        if (mode == IGRAPH_OUT)
        {
            // Get indices of where neighbours are
            size_t cum_degree_this_node = (size_t) VECTOR(this->_graph->os)[node];
            size_t cum_degree_next_node = (size_t) VECTOR(this->_graph->os)[node + 1];
            // Get a random index from them
            size_t rand_neigh_idx = igraph_rng_get_integer(igraph_rng_default(), cum_degree_this_node, cum_degree_next_node - 1);
            // Return the neighbour at that index
#ifdef DEBUG
            cerr << "Degree: " << this->degree(node, mode) << " diff in cumulative: " << cum_degree_next_node - cum_degree_this_node << endl;
#endif
            rand_neigh = VECTOR(this->_graph->to)[ (size_t)VECTOR(this->_graph->oi)[rand_neigh_idx] ];
        }
        else if (mode == IGRAPH_IN)
        {
            // Get indices of where neighbours are
            size_t cum_degree_this_node = (size_t) VECTOR(this->_graph->is)[node];
            size_t cum_degree_next_node = (size_t) VECTOR(this->_graph->is)[node + 1];
            // Get a random index from them
            size_t rand_neigh_idx = igraph_rng_get_integer(igraph_rng_default(), cum_degree_this_node, cum_degree_next_node - 1);
#ifdef DEBUG
            cerr << "Degree: " << this->degree(node, mode) << " diff in cumulative: " << cum_degree_next_node - cum_degree_this_node << endl;
#endif
            // Return the neighbour at that index
            rand_neigh = VECTOR(this->_graph->from)[ (size_t)VECTOR(this->_graph->ii)[rand_neigh_idx] ];
        }
    }
    else
    {
        // both in- and out- neighbors in a directed graph.
        size_t cum_outdegree_this_node = (size_t)VECTOR(this->_graph->os)[node];
        size_t cum_indegree_this_node  = (size_t)VECTOR(this->_graph->is)[node];

        size_t cum_outdegree_next_node = (size_t)VECTOR(this->_graph->os)[node + 1];
        size_t cum_indegree_next_node  = (size_t)VECTOR(this->_graph->is)[node + 1];

        size_t total_outdegree = cum_outdegree_next_node - cum_outdegree_this_node;
        size_t total_indegree = cum_indegree_next_node - cum_indegree_this_node;

        size_t rand_idx = igraph_rng_get_integer(igraph_rng_default(), 0, total_outdegree + total_indegree - 1);

#ifdef DEBUG
        cerr << "Degree: " << this->degree(node, mode) << " diff in cumulative: " << total_outdegree + total_indegree << endl;
#endif
        // From among in or out neighbours?
        if (rand_idx < total_outdegree)
        {
            // From among outgoing neighbours
            size_t rand_neigh_idx = cum_outdegree_this_node + rand_idx;
            rand_neigh = VECTOR(this->_graph->to)[ (size_t)VECTOR(this->_graph->oi)[rand_neigh_idx] ];
        }
        else
        {
            // From among incoming neighbours
            size_t rand_neigh_idx = cum_indegree_this_node + rand_idx - total_outdegree;
            rand_neigh = VECTOR(this->_graph->from)[ (size_t)VECTOR(this->_graph->ii)[rand_neigh_idx] ];
        }
    }

    return rand_neigh;
}

#include <iostream>
using namespace std;
/****************************************************************************
  Creates a graph with communities as node and links as weights between
  communities.

  The weight of the edges in the new graph is simply the sum of the weight
  of the edges between the communities. The self weight of a node (i.e. the
  weight of its self loop) is the internal weight of a community. The size
  of a node in the new graph is simply the size of the community in the old
  graph.
*****************************************************************************/
Graph* Graph::collapse_graph(MutableVertexPartition* partition)
{
#ifdef DEBUG
    cerr << "Graph* Graph::collapse_graph(vector<size_t> membership)" << endl;
#endif
    //size_t n = this->vcount();
    size_t m = this->ecount();

    vector< map<size_t, double> > collapsed_edge_weights(partition->nb_communities());

    igraph_integer_t v, u;
    for (size_t e = 0; e < m; e++)
    {
#ifdef MATLAB_SUPPORT
            ctrlcCheckPoint(__FILE__, __LINE__); // Interrupt here
#endif
        double w = this->edge_weight(e);
        igraph_edge(this->_graph, e, &v, &u);
        size_t v_comm = partition->membership((size_t)v);
        size_t u_comm = partition->membership((size_t)u);
        if (collapsed_edge_weights[v_comm].count(u_comm) > 0)
            collapsed_edge_weights[v_comm][u_comm] += w;
        else
            collapsed_edge_weights[v_comm][u_comm] = w;
    }

    // Now create vector for edges, first determined the number of edges
    size_t m_collapsed = 0;
    size_t n_collapsed = partition->nb_communities();

    for (vector< map<size_t, double> >::iterator itr = collapsed_edge_weights.begin();
            itr != collapsed_edge_weights.end(); itr++)
    {
        m_collapsed += itr->size();
    }

    igraph_vector_t edges;
    vector<double> collapsed_weights(m_collapsed, 0.0);
    double total_collapsed_weight = 0.0;

    igraph_vector_init(&edges, 2 * m_collapsed); // Vector or edges with edges (edge[0], edge[1]), (edge[2], edge[3]), etc...

    size_t e_idx = 0;
    for (size_t v = 0; v < n_collapsed; v++)
    {
#ifdef MATLAB_SUPPORT
            ctrlcCheckPoint(__FILE__, __LINE__); // Interrupt here
#endif
        for (map<size_t, double>::iterator itr = collapsed_edge_weights[v].begin();
                itr != collapsed_edge_weights[v].end(); itr++)
        {
            size_t u = itr->first;
            double w = itr->second;
            VECTOR(edges)[2 * e_idx] = v;
            VECTOR(edges)[2 * e_idx + 1] = u;
            collapsed_weights[e_idx] = w;
            total_collapsed_weight += w;
            if (e_idx >= m_collapsed)
                throw Exception("Maximum number of possible edges exceeded.");
            // next edge
            e_idx += 1;
        }
    }

    double const eps = 1e-6;
    if (fabs(total_collapsed_weight - this->total_weight()) > eps)
        throw Exception("Total collapsed weight is not equal to original weight.");

    // Create graph based on edges
    igraph_t* graph = new igraph_t();
    igraph_create(graph, &edges, n_collapsed, this->is_directed());
    igraph_vector_destroy(&edges);

    if ((size_t) igraph_vcount(graph) != partition->nb_communities())
        throw Exception("Something went wrong with collapsing the graph.");

    // Calculate new node sizes
    vector<size_t> csizes(n_collapsed, 0);
    for (size_t c = 0; c < partition->nb_communities(); c++)
        csizes[c] = partition->csize(c);

    Graph* G = new Graph(graph, collapsed_weights, csizes);
    G->_remove_graph = true;
#ifdef DEBUG
    cerr << "exit Graph::collapse_graph(vector<size_t> membership)" << endl << endl;
#endif
    return G;
}

// ##############################################
Graph * init(double *W, int N, int M)
{
    Eigen::MatrixXd EW = Eigen::Map<Eigen::MatrixXd>(W, N, M); //check N,M rows or cols
    bool feedingSparseMatrix = false;
    if (N > 3 && M == 3)
    {
        feedingSparseMatrix = true;
    }
    vector<double> edges_weights;
    vector<double> edges_list;
    if (feedingSparseMatrix) // the input matrix is a sparse matrix with 2 or 3 columns. If 2 columns the edges list is given, if 3 columns the third column is the edge weight
    {
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> B1, B2;
        B1 = (EW.col(0).array() >= EW.col(1).array()).cast<int>();
        B2 = (EW.col(0).array() < EW.col(1).array()).cast<int>();
        bool isUpperTriangular = false;
        bool isLowerTriangular = false;
        bool isSymmetric = false;
        int sum1 = B1.sum();
        int sum2 = B2.sum();
        // Condizione semplice da verificare facendo [i j w]=find(A), oppure [i j w]=find(triu(A)) oppure [i j w]=find(tril(A))
        if (sum1 == sum2)
            isSymmetric = true;
        if (sum1 == 0 && sum2 == M)
            isUpperTriangular = true;
        if (sum1 == M && sum2 == 0)
            isLowerTriangular = true;

        if (!isSymmetric && !isUpperTriangular && !isLowerTriangular)
        {
            throw Exception("Matrix is not symmetric, nor triangular lower or upper triangular. Check diagonal and non symmetric values.");
        }
        for (int l = 0; l < M; ++l)
        {
            double row_node = EW(l, 0); //index of row from MATLAB find command
            double column_node = EW(l, 1); //index of column from MATLAB find command
            double w = EW(l, 2);
            if ( isUpperTriangular || isLowerTriangular) // keeps only symmetric and also avoid self-loops (implicitly inserting upper triangular)
            {
                edges_list.push_back(column_node - 1);
                edges_list.push_back(row_node - 1);
                edges_weights.push_back(w);
            }
            else if (isSymmetric)
            {
                if (row_node < column_node)
                {
                    edges_list.push_back(column_node - 1);
                    edges_list.push_back(row_node - 1);
                    edges_weights.push_back(w);
                }
            }
        }
    }
    else // the input matrix is square adjacency matrix
    {
        if (EW.trace() > 0)
        {
            throw Exception("Adjacency matrix has self loops, only simple graphs allowed.");
        }
        for (int i = 0; i < N; ++i)
        {
#ifdef MATLAB_SUPPORT
            ctrlcCheckPoint(__FILE__, __LINE__); // Interrupt here
#endif
            for (int j = i + 1; j < N; j++)
            {
                double w = std::max(EW.coeffRef(i, j), EW.coeffRef(j, i));
                if (w > 0)
                {
                    edges_list.push_back(i);
                    edges_list.push_back(j);
                    edges_weights.push_back(w);
                }
                if (w < 0)
                {
                    throw Exception("Negative edge weight found. Only positive weights supported.");
                }
            }
        }
    }

    if (edges_list.empty())
        throw Exception("Empty graph provided.");

    // Create the Graph object from the igraph data structure
    // Fill the edges into the igraph IG
    igraph_t *IG = new igraph_t;
    igraph_vector_t igedges_list;
    igraph_vector_view(&igedges_list, edges_list.data(), edges_list.size());
    igraph_create(IG, &igedges_list, std::max(N,M), 0);

    Graph *G = new Graph(IG, edges_weights);
    
    //igraph_destroy(&IG);

    return G;
}
