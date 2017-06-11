/* This file is part of MULTILOUVAIN, a program to find network partitions
*
*  Copyright (C) 2017 Carlo Nicolini <carlo.nicolini@iit.it>
*
*  MULTILOUVAIN is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 3 of the License, or (at your option) any later version.
*
*  Alternatively, you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  MULTILOUVAIN is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License and a copy of the GNU General Public License along with
*  MULTILOUVAIN. If not, see <http://www.gnu.org/licenses/>.
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
#include "DCSurpriseVertexPartition.h"
#include "igraph_utils.h"
#include "SurpriseVertexPartition.h"

#include <Eigen/Core>


using namespace std;

void exit_with_help()
{
    std::printf(
                "Usage: run_multilouvain [options] graph_file\n"
                "Version 0.3 20 March 2017\n"
                "graph_file the file containing the graph. Accepted formats are adjacency matrices in csv format or"
                "\nedges list (the ncol format), additionally with a third column with edge weights\n"
                "\n"
                "Options:\n"
                "-q [quality]\n"
                "   0 Asymptotical Surprise\n"
                "   1 Significance\n"
                "   2 Reichardt-Bornholdt Erdos Renyi model\n"
                "   3 Reichardt-Bornholdt Configuration Model\n"
                "   4 Constant Potts Model\n"
                "   5 Newman Modularity\n"
                "   6 Degree Corrected Surprise\n"
                "-c [consider_comms]:\n"
                "		1: ALL_COMMS. Looks for improvements in all communities. Slower but better results.\n"
                "		2: ALL_NEIGH_COMMS Looks for improvements just in neighboring communities. Faster and still good results.\n"
                "		3: RAND_COMMS Choose a random communitiy to look for improvements.\n"
                "		4: RAND_NEIGH_COMMS Choose a random neighbor community. Fastest option but poor results.\n"
                "-v [report_level] ERROR=0, WARNING=1, INFO=2, DEBUG=3, DEBUG1=4, DEBUG2=5, DEBUG3=6, DEBUG4=7\n"
                "-s [seed] specify the random seed, default time(0)\n"
                "-b [bool] wheter to start with initial random cluster or every node in its community\n"
                "-r [repetitions], number of , default=1\n"
                "-p [print solution]\n"
                "\n"
                );
    exit(1);
}

/**
 * @brief parse_command_line
 * @param argc
 * @param argv
 * @param input_file_name
 */
LouvainParams parse_command_line(int argc, char **argv)
{
    LouvainParams params;
    int i=1;
    for(i=1; i<argc; i++)
    {
        if(argv[i][0] != '-')
            break;
        if(++i>=argc)
            exit_with_help();
        switch(argv[i-1][1])
        {
        case 'v':
        case 'V':
        {
            params.verbosity_level = atoi(argv[i]);
            if (params.verbosity_level>7)
                params.verbosity_level=7;
            break;
        }
        case 's':
        case 'S':
        {
            params.rand_seed = atoi(argv[i]);
            break;
        }
        case 'q':
        case 'Q':
        {
            if (atoi(argv[i])<QualitySurprise || atoi(argv[i])>QualityDCSurprise)
            {
                exit_with_help();
            }
            else
            {
                params.quality = static_cast<QualityFunction>(atoi(argv[i]));
            }
            break;
        }
        case 'r':
        case 'R':
        {
            params.random_order = atoi(argv[i]);
            break;
        }
        case 'c':
        case 'C':
        {
            params.consider_comms = atoi(argv[i]);
            if ( params.consider_comms <0 || params.consider_comms > 4)
            {
                exit_with_help();
            }
            break;
        }
        default:
            exit_with_help();
        }
    }

    // Determine filenames
    if(i>=argc)
        exit_with_help();

    char input[1024];
    strcpy(input, argv[i]);
    params.filename = std::string(input);

    std::ifstream is(input);
    if (!is.good())
    {
        std::cerr << std::string("File \"" + params.filename + "\" not found") << std::endl;
        exit_with_help();
    }
    return params;
}

/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[])
{
    try
    {
        LouvainParams params = parse_command_line(argc,argv);
        // Set the random seed on the current time in microseconds, if not specified
        if (params.rand_seed == -1)
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
            srand(params.rand_seed); // initialize random seed from parameters
        }

        Eigen::MatrixXd W = read_adj_matrix<Eigen::MatrixXd>(params.filename,' ');

        Graph *G = init(W.data(),W.rows(),W.cols());
        if (params.verbosity_level > 0)
            cout << "|V|=" << G->vcount() << " |E|=" << G->ecount() << endl;
        MutableVertexPartition *partition=NULL;
        switch ( params.quality )
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
            partition = new CPMVertexPartition(G, params.cpmgamma);
            break;
        }
        case QualityModularity:
        {
            partition = new ModularityVertexPartition(G);
            break;
        }
        case QualityDCSurprise:
        {
            partition = new DCSurpriseVertexPartition(G);
            break;
        }
        default:
            exit_with_help();
        }

        Optimiser *opt = new Optimiser;

        double qual = opt->optimize_partition(partition);
        cout << qual << endl;
        cout << partition->membership() << endl;

        delete opt;
        delete partition;
        G->dispose();
        delete G;
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
        exit_with_help();
    }

    return 0;
}

