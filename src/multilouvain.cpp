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
#include "igraph_utils.h"

#include <Eigen/Core>

#ifdef __linux__
#include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;
#define IS_ADJACENCY_MATRIX(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && (mxIsDouble(P) || mxIsUint8(P) || mxIsUint8(P) || mxIsLogical(P)))


/**
 * @brief printUsage
 */
void printUsage()
{
    mexPrintf("LOUVAIN Louvain Algorithm for community detection.\n");
    mexPrintf("Version 0.2.2 17 March 2017\n");
    mexPrintf("[membership, qual] = multilouvain(W);\n");
    mexPrintf("Input:\n");
    mexPrintf("	W: an undirected weighted network with positive edge weights. Negative edge weights are not supported and an error is thrown. Remember to use real matrices as logical matrices throw error.\n");
    mexPrintf(" If the input matrix is not symmetric the maximum element between edge (i,j) and edge (j,i) is taken as policy.");
    mexPrintf("Output:\n");
    mexPrintf("	membership: the membership vector that represents the community which every vertex belongs to after quality optimization.\n");
    mexPrintf("	qual: the current quality of the partition.\n");
    mexPrintf("\n");
    mexPrintf("Options:\n");
    mexPrintf("multilouvain accepts additional arguments to control the optimization process\n");
    mexPrintf("[m, qual] = multilouvain(W,'quality',val);\n");
    mexPrintf("	val is one of the following integers: {0,1,2,3,4,5,6}:\n");
    mexPrintf("		0: Asymptotical Surprise\n");
    mexPrintf("		1: Significance\n");
    mexPrintf("		2: Reichardt-Bornholdt with Erdos-Renyi null model\n");
    mexPrintf("		3: Reichardt-Bornholdt with configuration model null model\n");
    mexPrintf("		4: Constant Potts Model (gamma=0.5), to specify as further argument.\n");
    mexPrintf("		5: Newman's modularity\n");
    mexPrintf("		6: Degree Corrected Surprise\n");
    mexPrintf("[m, qual] = multilouvain(W,'consider_comms',val);\n");
    mexPrintf("	consider_comms is one of the following integers: {1,2,3,4}:\n");
    mexPrintf("		1: ALL_COMMS. Looks for improvements in all communities. Slower but better results.\n");
    mexPrintf("		2: ALL_NEIGH_COMMS Looks for improvements just in neighboring communities. Faster and still good results.\n");
    mexPrintf("		2: RAND_COMMS Choose a random communitiy to look for improvements.\n");
    mexPrintf("		4: RAND_NEIGH_COMMS Choose a random neighbor community. Fastest option but poor results.\n");
    mexPrintf("[m, qual] = multilouvain(W,'cpm_gamma',val);\n");
    mexPrintf("		val: is the resolution parameter for CPM method.\n");
    mexPrintf("[m, qual] = multilouvain(W,'max_itr',val);\n");
    mexPrintf("		val: is the maximum number of iterations of Louvain algorithm. Default 100000.\n");
    mexPrintf("[m, qual] = multilouvain(W,'random_order',val);\n");
    mexPrintf("		val: whether to consider nodes in random order or not. Default true.\n");
    mexPrintf("[m, qual] = multilouvain(W,'seed',val);\n");
    mexPrintf("		val: to provide a specific random seed to the algorithm, in order to have reproducible results.\n");
}

/**
 * @brief The error_type enum
 */
enum error_type
{
    NO_ERROR = 0,
    ERROR_TOO_MANY_OUTPUT_ARGS = 1,
    ERROR_NOT_ENOUGH_ARGS = 2,
    ERROR_ARG_VALUE = 3,
    ERROR_ARG_TYPE = 4,
    ERROR_MATRIX = 5,
    ERROR_ARG_EMPTY = 6,
    ERROR_ARG_UNKNOWN = 7
};

static const char *error_strings[] =
{
    "",
    "Too many output arguments.",
    "Not enough input arguments.",
    "Non valid argument value.",
    "Non valid argument type.",
    "Non valid input adjacency matrix. MULTILOUVAIN accepts symmetric real dense-type (n x n) matrices or sparse edges-list representation \
    [num_edges x 3] array of edges list with edge endpoints and weight.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};



/**
 * @brief parse_args
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 * @param pars
 * @param argposerr
 * @return
 */
error_type parse_args(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[], LouvainParams *pars, int *argposerr )
{
    if (nInputArgs < 1)
    {
        *argposerr = 0;
        return ERROR_NOT_ENOUGH_ARGS;
    }

    if (nOutputArgs > 2)
    {
        *argposerr = 0;
        return ERROR_TOO_MANY_OUTPUT_ARGS;
    }

    const mxArray *W = inputArgs[0];
    int M = mxGetM(W);
    int N = mxGetN(W);
    // In this case we are feeding instead of the adjacency matrix, the result of [i j w]=find(A);
    bool feedingSparseMatrix = false;
    if (N == 3 && M > 3)
    {
        feedingSparseMatrix = true;
    }

    bool v1 = M != N;
    if (feedingSparseMatrix)
        v1 = false;
    bool v2 = mxIsComplex(W);
    bool v3 = mxIsEmpty(W);
    bool v4 = mxIsCell(W);
    bool v5 = !mxIsNumeric(W);
    //bool v6 = !mxIsSparse(W);

    if ( v1 || v2 || v3 || v4 || v5 )
    {
        *argposerr = 0;
        return ERROR_MATRIX;
    }

    // Iterate on function arguments
    int argcount = 1;
    while (argcount < nInputArgs)
    {
        // Be sure that something exists after c-th argument
        if (argcount + 1 >= nInputArgs)
        {
            *argposerr = argcount;
            return ERROR_ARG_EMPTY;
        }
        // Couple argument type - argument value
        const mxArray *partype = inputArgs[argcount];
        const mxArray *parval = inputArgs[argcount + 1];
        char* cpartype;
        // To be a valid parameter specification it must be a couple ['char',real]
        if (mxIsChar(partype) && !mxIsChar(parval))
        {
            cpartype = mxArrayToString(partype);
#ifdef _DEBUG
            mexPrintf("ARGUMENT: %s VALUE=%g\n", cpartype, *mxGetPr(parval));
#endif
            // Parse string value inputArgs[c]
            if ( strcasecmp(cpartype, "quality") == 0 )
            {
                pars->quality = static_cast<QualityFunction>(*mxGetPr(parval));
                if (pars->quality < 0 || pars->quality > QualityDCSurprise )
                {
                    *argposerr = argcount + 1;
                    return ERROR_ARG_VALUE;
                }
                argcount += 2;
            }
            else if ( strcasecmp(cpartype, static_cast<const char*>("consider_comms")) == 0 )
            {
                pars->consider_comms = static_cast<int>(*mxGetPr(parval));
                if (pars->consider_comms < 1 || pars->consider_comms > 4)
                {
                    *argposerr = argcount + 1;
                    return ERROR_ARG_VALUE;
                }
                argcount += 2;
            }
            else if ( strcasecmp(cpartype, static_cast<const char*>("cpm_gamma")) == 0 )
            {
                pars->cpmgamma = (*mxGetPr(parval));
                argcount += 2;
            }
            else if ( strcasecmp(cpartype, static_cast<const char*>("delta")) == 0 )
            {
                pars->delta = *mxGetPr(parval);
                argcount += 2;
            }
            else if ( strcasecmp(cpartype, static_cast<const char*>("max_itr")) == 0 )
            {
                pars->max_itr = static_cast<size_t>(std::floor(*mxGetPr(parval)));
                argcount += 2;
            }
            else if ( strcasecmp(cpartype, static_cast<const char*>("random_order")) == 0 )
            {
                pars->random_order = static_cast<bool>((*mxGetPr(parval)) > 0);
                argcount += 2;
            }
            else if ( strcasecmp(cpartype, static_cast<const char*>("seed")) == 0 )
            {
                pars->rand_seed = static_cast<int>(std::floor(*mxGetPr(parval)));
                argcount += 2;
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

/**
 * @brief mexFunction
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 */
void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    LouvainParams pars;
    // Set default values for parameters
    pars.quality = QualitySurprise;
    pars.consider_comms = Optimiser::ALL_COMMS;
    pars.eps = 1E-5;
    pars.cpmgamma = 0.5;
    pars.delta = 1E-2;
    pars.max_itr = 1E4;
    pars.random_order = true;
    pars.rand_seed = -1; // default value for the random seed, if -1 than microseconds time is used.

    // Check the arguments of the function
    int error_arg_pos = -1;
    error_type err = parse_args(nOutputArgs, outputArgs, nInputArgs, inputArgs, &pars, &error_arg_pos);

    if (err != NO_ERROR)
    {
        std::stringstream ss;
        ss << "Error at " << error_arg_pos << "-th argument: " << error_strings[err] ;
        if (err == ERROR_NOT_ENOUGH_ARGS)
            printUsage();
        mexErrMsgTxt(ss.str().c_str());
    }

#ifdef _DEBUG
    printf("[INFO] Method=%d\n[INFO] Consider_comms=%d\n[INFO] CPMgamma=%f\n[INFO] Delta=%f\n[INFO] Max_itr=%zu\n[INFO] Random_order=%d rand_seed=%d\n", pars.method, pars.consider_comms, pars.cpmgamma, pars.delta, pars.max_itr, pars.random_order, pars.rand_seed);
#endif

    // Set the random seed on the current time in microseconds, if not specified
    if (pars.rand_seed == -1)
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

    // Try to understand the type of the matrix if sparse or not
    // Get number of vertices in the network
    int ncols = mxGetN(inputArgs[0]); // number of columns
    int nrows = mxGetM(inputArgs[0]); // number of rows
    double *W = mxGetPr(inputArgs[0]);
    
    try
    {
        Graph *G = init(W, nrows, ncols);
        cerr << G->vcount() << " " << G->ecount() << endl;
        // Now that the Graph has been built it's time to run the solver.
        // Create the partition function
        MutableVertexPartition *partition;
        // Create the optimizer instance
        Optimiser *opt = new Optimiser;

        // Allocate partition method and optimization
        switch ( pars.quality )
        {
        case QualitySurprise:
        {
            partition = new SurpriseVertexPartition(G);
            break;
        }
        case QualitySignificance:
        {
            partition = new SignificanceVertexPartition(G);
            break;
        }
        case QualityRBER:
        {
            partition = new RBERVertexPartition(G);
            break;
        }
        case QualityRBConfiguration:
        {
            partition = new RBConfigurationVertexPartition(G);
            break;
        }
        case QualityCPM:
        {
            partition = new CPMVertexPartition(G, pars.cpmgamma);
            break;
        }
        case QualityModularity:
        {
            partition = new ModularityVertexPartition(G);
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
        // Prepare output
        outputArgs[0] = mxCreateDoubleMatrix(1, (mwSize)G->vcount(), mxREAL);
        double *memb = mxGetPr(outputArgs[0]);
        // Copy the membership vector to outputArgs[0] which has been already preallocated
        for (size_t i = 0; i < G->vcount() ; ++i)
            memb[i] = static_cast<double>(partition->membership(i) + 1);

        // Copy the value of partition quality
        outputArgs[1] = mxCreateDoubleScalar(qual);

        // Cleanup the memory (follow this order)
        delete opt;
        delete partition;
        G->dispose(); // must manually clear the underlying igraph_t
        delete G; // then delete the pointer
    }
    catch (std::exception &e)
    {
        mexErrMsgTxt(e.what());
    }

    // Finish the function
    return;
}
