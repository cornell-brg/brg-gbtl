/*
 * GraphBLAS Template Library (GBTL), Version 3.0
 *
 * Copyright 2020 Carnegie Mellon University, Battelle Memorial Institute, and
 * Authors.
 *
 * THIS MATERIAL WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY AN AGENCY OF
 * THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES GOVERNMENT NOR THE
 * UNITED STATES DEPARTMENT OF ENERGY, NOR THE UNITED STATES DEPARTMENT OF
 * DEFENSE, NOR CARNEGIE MELLON UNIVERSITY, NOR BATTELLE, NOR ANY OF THEIR
 * EMPLOYEES, NOR ANY JURISDICTION OR ORGANIZATION THAT HAS COOPERATED IN THE
 * DEVELOPMENT OF THESE MATERIALS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
 * OR USEFULNESS OR ANY INFORMATION, APPARATUS, PRODUCT, SOFTWARE, OR PROCESS
 * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED
 * RIGHTS.
 *
 * Released under a BSD-style license, please see LICENSE file or contact
 * permission@sei.cmu.edu for full terms.
 *
 * [DISTRIBUTION STATEMENT A] This material has been approved for public release
 * and unlimited distribution.  Please see Copyright notice for non-US
 * Government use and distribution.
 *
 * DM20-0442
 */

#include <iostream>
#include <fstream>
#include <chrono>

#define GRAPHBLAS_DEBUG 1

#include <graphblas/graphblas.hpp>
#include <algorithms/cluster_louvain.hpp>
#include "Timer.hpp"
#include "read_edge_list.hpp"

//****************************************************************************
int main(int argc, char **argv)
{
#if 0
    if (argc < 2)
    {
        std::cerr << "ERROR: too few arguments." << std::endl;
        exit(1);
    }

    // Read the edgelist and create the tuple arrays
    grb::IndexArrayType iA;
    grb::IndexArrayType jA;
    std::vector<double> weights;
    grb::IndexType NUM_NODES(read_triples<double>(argv[1], iA, jA, weights));

    using MatType = grb::Matrix<double>;

    MatType A(NUM_NODES, NUM_NODES);
    A.build(iA.begin(), jA.begin(), weights.begin(), iA.size());

    std::cout << "#Nodes = " << NUM_NODES << std::endl;
    std::cout << "#Edges = " << A.nvals() << std::endl;

    std::cout << "Running louvain clustering (2 ways)..." << std::endl;

    Timer<std::chrono::steady_clock> my_timer;

    // Perform clustering with 2 different algorithms
    //===================
    my_timer.start();
    auto cluster_matrix = algorithms::louvain_cluster(A);
    my_timer.stop();

    std::cout << "Elapsed time: " << my_timer.elapsed() << " msec." << std::endl;
    auto cluster_assignments =
        algorithms::get_louvain_cluster_assignments(cluster_matrix);
    print_vector(std::cout, cluster_assignments, "cluster assignments");

    //===================
    my_timer.start();
    auto cluster2_matrix = algorithms::louvain_cluster_masked(A);
    my_timer.stop();

    std::cout << "Elapsed time: " << my_timer.elapsed() << " msec." << std::endl;
    auto cluster2_assignments =
        algorithms::get_louvain_cluster_assignments(cluster2_matrix);
    print_vector(std::cout, cluster2_assignments, "cluster (masked)  assignments");
#endif
    return 0;
}
