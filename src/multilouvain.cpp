/* This file is part of FAGSO, a program to find network partitions
*
*  Copyright (C) 2014-2015 Carlo Nicolini <carlo.nicolini@iit.it>
*
*  FAGSO is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 3 of the License, or (at your option) any later version.
*
*  Alternatively, you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  FAGSO is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License and a copy of the GNU General Public License along with
*  FAGSO. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <string.h>
#include <cmath>
#include <sstream>
#include <igraph.h>
#include <sys/time.h>

#include "GraphHelper.h"
#include "Optimiser.h"
#include "MutableVertexPartition.h"
#include "SurpriseVertexPartition.h"
#include "SignificanceVertexPartition.h"
#include "RBConfigurationVertexPartition.h"
#include "RBERVertexPartition.h"
#include "CPMVertexPartition.h"
#include "ModularityVertexPartition.h"
#include "KLModularityVertexPartition.h"
#include "igraph_utils.h"

#ifdef __linux__
    #include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;
#define IS_ADJACENCY_MATRIX(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && (mxIsDouble(P) || mxIsUint8(P) || mxIsUint8(P) || mxIsLogical(P)))

enum LouvainMethod
{
    MethodSurprise = 0,
    MethodSignificance = 1,
    MethodRBER = 2,
    MethodRBConfiguration = 3,
    MethodCPM = 4 ,
    MethodModularity = 5,
    MethodKLModularity = 6
};


void printUsage()
{
    mexPrintf("LOUVAIN Louvain Algorithm for community detection.\n");
    mexPrintf("[membership, qual] =louvain(W);\n");
    mexPrintf("Input:\n");
    mexPrintf("	W: an undirected weighted network with positive edge weights. Negative edge weights are not seen as edges and therefore discarded. Remember to use real matrices, logical matrices throw error.\n");
    mexPrintf("Output:\n");
    mexPrintf("	membership: the membership vector that represents the community which every vertex belongs to after quality optimization.\n");
    mexPrintf("	qual: the current quality of the partition.\n");
    mexPrintf("\n");
    mexPrintf("Options:\n");
    mexPrintf("louvain accepts additional arguments to control the optimization process\n");
    mexPrintf("[m, qual] = louvain(W,'method',val);\n");
    mexPrintf("	val is one of the following integers: {1,2,3,4,5}:\n");
    mexPrintf("		0: Surprise\n");
    mexPrintf("		1: Significance\n");
    mexPrintf("		2: Reichardt-Bornholdt with Erdos-Renyi null model\n");
    mexPrintf("		3: Reichardt-Bornholdt with configuration model null model\n");
    mexPrintf("		4: Constant Potts Model (gamma=0.5), to specify as further argument.\n");
    mexPrintf("		5: Newman's modularity\n");
    mexPrintf("[m, qual] = louvain(W,'consider_comms',val);\n");
    mexPrintf("	consider_comms is one of the following integers: {1,2,3,4}:\n");
    mexPrintf("		1: ALL_COMMS. Looks for improvements in all communities. Slower but better results.\n");
    mexPrintf("		2: ALL_NEIGH_COMMS Looks for improvements just in neighboring communities. Faster and still good results.\n");
    mexPrintf("		2: RAND_COMMS Choose a random communitiy to look for improvements.\n");
    mexPrintf("		3: RAND_NEIGH_COMMS Choose a random neighbor community. Fastest option but poor results.\n");
    mexPrintf("[m, qual] = louvain(W,'cpm_gamma',val);\n");
    mexPrintf("		val: is the resolution parameter for CPM method.\n");
    mexPrintf("[m, qual] = louvain(W,'max_itr',val);\n");
    mexPrintf("		val: is the maximum number of iterations of Louvain algorithm. Default 100000.\n");
    mexPrintf("[m, qual] = louvain(W,'random_order',val);\n");
    mexPrintf("		val: whether to consider nodes in random order or not. Default true.\n");
    mexPrintf("[m, qual] = louvain(W,'seed',val);\n");
    mexPrintf("		val: to provide a specific random seed to the algorithm, in order to have reproducible results.\n");
    mexPrintf("\n\n");
    mexPrintf("Example:\n");
    mexPrintf(">> A=rand(100,100); A=(A+A')/2; A=A.*(A>0.5);\n");
    mexPrintf(">> [memb, qual] = louvain(A,'method',2,'consider_comms',2);\n");
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

static const char *error_strings[] =
{
    "",
    "Too many output arguments.",
    "Not enough input arguments.",
    "Non valid argument value.",
    "Non valid argument type.",
    "Non valid input adjacency matrix. Must be symmetric real square matrix.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};

struct LouvainParams
{
    LouvainMethod method;
    int consider_comms;  // Indicates how communities will be considered for improvement. Should be one of the parameters below
    double eps;          // If the improvement falls below this threshold, stop iterating.
    double delta;        // If the number of nodes that moves falls below this threshold, stop iterating.
    size_t max_itr;      // Maximum number of iterations to perform.
    int random_order;    // If True the nodes will be traversed in a random order when optimising a quality function.
    double cpmgamma; // resolution parameter in  CPM
    unsigned int rand_seed; // random seed for the louvain algorithm
};

error_type parse_args(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[], LouvainParams *pars, int *argposerr )
{
    if (nInputArgs < 1)
    {
        *argposerr = 0;
        return ERROR_NOT_ENOUGH_ARGS;
    }

    if (nOutputArgs>2)
    {
        *argposerr = 0;
        return ERROR_TOO_MANY_OUTPUT_ARGS;
    }

    const mxArray *W = inputArgs[0];

    int M = mxGetM(W);
    int N = mxGetN(W);

    if (M!=N  || mxIsComplex(W) || mxIsEmpty(W) || mxIsCell(W) || !mxIsNumeric(W))
    {
        *argposerr = 0;
        return ERROR_MATRIX;
    }

    // Iterate on function arguments
    int argcount=1;
    while (argcount<nInputArgs)
    {
        // Be sure that something exists after c-th argument
        if (argcount+1 >= nInputArgs)
        {
            *argposerr = argcount;
            return ERROR_ARG_EMPTY;
        }
        // Couple argument type - argument value
        const mxArray *partype = inputArgs[argcount];
        const mxArray *parval = inputArgs[argcount+1];
        char* cpartype;
        // To be a valid parameter specification it must be a couple ['char',real]
        if (mxIsChar(partype) && !mxIsChar(parval))
        {
            cpartype = mxArrayToString(partype);
#ifdef _DEBUG
            mexPrintf("ARGUMENT: %s VALUE=%g\n", cpartype,*mxGetPr(parval));
#endif
            // Parse string value inputArgs[c]
            if ( strcasecmp(cpartype,"Method")==0 )
            {
                pars->method = static_cast<LouvainMethod>(*mxGetPr(parval));
                if (pars->method<0 || pars->method>MethodKLModularity)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("consider_comms"))==0 )
            {
                pars->consider_comms = static_cast<int>(*mxGetPr(parval));
                if (pars->consider_comms <1 || pars->consider_comms > 4)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("cpm_gamma"))==0 )
            {
                pars->cpmgamma = (*mxGetPr(parval));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("delta"))==0 )
            {
                pars->delta = *mxGetPr(parval);
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("max_itr"))==0 )
            {
                pars->max_itr = static_cast<size_t>(std::floor(*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("random_order"))==0 )
            {
                pars->random_order = static_cast<bool>((*mxGetPr(parval))>0);
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("seed"))==0 )
            {
                pars->rand_seed = static_cast<int>(std::floor(*mxGetPr(parval)));
                argcount+=2;
            }
            else
            {
                *argposerr = argcount;
                return ERROR_ARG_UNKNOWN;
            }
        }
        else //else return the position of the argument and type of error
        {
            *argposerr = argcount;
            return ERROR_ARG_TYPE;
        }
        mxFree(cpartype); // free the converted argument
    }
    return NO_ERROR;
}

void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    LouvainParams pars;
    // Set default values for parameters
    pars.method = MethodSurprise;
    pars.consider_comms = Optimiser::ALL_COMMS;
    pars.eps = 1E-5;
    pars.cpmgamma = 0.5;
    pars.delta = 1E-2;
    pars.max_itr = 1E4;
    pars.random_order = true;
    pars.rand_seed = 0; // default value for the random seed, if 0 than microseconds time is used.

    // Check the arguments of the function
    int error_arg_pos=-1;
    error_type err = parse_args(nOutputArgs, outputArgs, nInputArgs, inputArgs, &pars, &error_arg_pos);

    if (err!=NO_ERROR)
    {
        std::stringstream ss;
        ss << "Error at " << error_arg_pos << "-th argument: " << error_strings[err] ;
        if (err == ERROR_NOT_ENOUGH_ARGS)
            printUsage();
        mexErrMsgTxt(ss.str().c_str());
    }

#ifdef _DEBUG
    printf("[INFO] Method=%d\n[INFO] Consider_comms=%d\n[INFO] CPMgamma=%f\n[INFO] Delta=%f\n[INFO] Max_itr=%zu\n[INFO] Random_order=%d rand_seed=%d\n",pars.method, pars.consider_comms, pars.cpmgamma, pars.delta, pars.max_itr, pars.random_order, pars.rand_seed);
#endif
    // Get number of vertices in the network
    int N = mxGetN(inputArgs[0]);

    // Set the random seed on the current time in microseconds, if not specified
    if (pars.rand_seed==0)
    {
#ifdef WIN32
        QueryPerformanceCounter(&endCount);
        std::srand(startCount.QuadPart);
#endif
#ifdef __linux__
        struct timeval start;
        gettimeofday(&start, NULL);
        std::srand(start.tv_usec);
#endif
#ifdef __apple__
        struct timeval start;
        gettimeofday(&start, NULL);
        std::srand(start.tv_usec);
#endif
    }
    else
    {
        srand(pars.rand_seed); // initialize random seed from parameters
    }

    // Fill the adjacency matrix of the graph
    igraph_matrix_t adj;
    igraph_matrix_view(&adj,mxGetPr(inputArgs[0]),N,N);

    vector<double> edge_weights;

    for (int i=0; i<N; ++i)
    {
        for (int j=i+1; j<N; j++)
        {
            double w = (*(mxGetPr(inputArgs[0])+i*N+j));
            if (w>0)
                edge_weights.push_back(w);
        }
    }

    // Create the graph from the adjacency matrix
    igraph_t graph;
    igraph_weighted_adjacency(&graph,&adj,IGRAPH_ADJ_UNDIRECTED,NULL,true);

    // Create the Graph helper object specifying edge weights too
    Graph *G;
    try
    {
        G  = new Graph(&graph,edge_weights);
    }
    catch (Exception &e)
    {
        //mexErrMsgTxt("Input network has diagonal entries. Set them to zero.");
        mexErrMsgTxt(e.what());
    }
    // Create the partition function
    MutableVertexPartition *partition;
    // Create the optimizer instance
    Optimiser *opt = new Optimiser;
    opt->consider_comms = Optimiser::ALL_COMMS;

    // Allocate partition method and optimization
    switch ( pars.method )
    {
    case MethodSurprise:
    {
        partition = new SurpriseVertexPartition(G);
        break;
    }
    case MethodSignificance:
    {
        partition = new SignificanceVertexPartition(G);
        break;
    }
    case MethodRBER:
    {
        partition = new RBERVertexPartition(G);
        break;
    }
    case MethodRBConfiguration:
    {
        partition = new RBConfigurationVertexPartition(G);
        break;
    }
    case MethodCPM:
    {
        partition = new CPMVertexPartition(G,pars.cpmgamma);
        break;
    }
    case MethodKLModularity:
    {
        partition = new KLModularityVertexPartition(G);
        break;
    }
    }
    // Set optimization things
    opt->consider_comms = pars.consider_comms;
    opt->random_order = pars.random_order;
    opt->delta = pars.delta;
    opt->max_itr = pars.max_itr;
    opt->eps = pars.eps;

    // Finally optimize the partition
    double qual = opt->optimize_partition(partition);
#ifdef _DEBUG
    mexPrintf("G=(V,E)=%d %d\n",G->vcount(),G->ecount());
    mexPrintf("|C|=%d Qual=%g\n",partition->nb_communities(),partition->quality());
#endif
    // Prepare output
    outputArgs[0] = mxCreateDoubleMatrix(1,(mwSize)N, mxREAL);
    double *memb = mxGetPr(outputArgs[0]);
    // Copy the membership vector to outputArgs[0] which has been already preallocated
    for (size_t i = 0; i<N ; ++i)
        memb[i] = static_cast<double>(partition->membership(i)+1);

    // Copy the value of partition quality
    outputArgs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *q = mxGetPr(outputArgs[1]);
    q[0] = qual;

#ifdef _DEBUG
    size_t n = partition->graph->total_size();
    size_t nc2 = partition->total_possible_edges_in_all_comms();
    double mc = partition->total_weight_in_all_comms();
    double m = partition->graph->total_weight();
    double qual = partition->quality();
#endif

    // Cleanup the memory (follow this order)
    delete opt;
    delete partition;
    delete G;
    igraph_destroy(&graph); //always clean the igraph object, because standard destructor of GraphHelper doesn't do it.

    // Finish the function
    return;
}
