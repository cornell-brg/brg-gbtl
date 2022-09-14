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
#include <filesystem>
#include <fstream>
#include <cassert>

#include <graphblas/graphblas.hpp>
#include <algorithms/bfs.hpp>

grb::IndexType const num_nodes = 34;
grb::IndexArrayType i = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,
    2,2,2,2,2,2,2,2,2,2,
    3,3,3,3,3,3,
    4,4,4,
    5,5,5,5,
    6,6,6,6,
    7,7,7,7,
    8,8,8,8,8,
    9,9,
    10,10,10,
    11,
    12,12,
    13,13,13,13,13,
    14,14,
    15,15,
    16,16,
    17,17,
    18,18,
    19,19,19,
    20,20,
    21,21,
    22,22,
    23,23,23,23,23,
    24,24,24,
    25,25,25,
    26,26,
    27,27,27,27,
    28,28,28,
    29,29,29,29,
    30,30,30,30,
    31,31,31,31,31,31,
    32,32,32,32,32,32,32,32,32,32,32,32,
    33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33};

grb::IndexArrayType j = {
    1,2,3,4,5,6,7,8,10,11,12,13,19,21,23,31,
    0,2,3,7,13,17,19,21,30,
    0,1,3,7,8,9,13,27,28,32,
    0,1,2,7,12,13,
    0,6,10,
    0,6,10,16,
    0,4,5,16,
    0,1,2,3,
    0,2,30,32,33,
    2,33,
    0,4,5,
    0,
    0,3,
    0,1,2,3,33,
    32,33,
    32,33,
    5,6,
    0,1,
    32,33,
    0,1,33,
    32,33,
    0,1,
    32,33,
    25,27,29,32,33,
    25,27,31,
    23,24,31,
    29,33,
    2,23,24,33,
    2,31,33,
    23,26,32,33,
    1,8,32,33,
    0,24,25,28,32,33,
    2,8,14,15,18,20,22,23,29,30,31,33,
    8,9,13,14,15,18,19,20,22,23,26,27,28,29,30,31,32};


//****************************************************************************
using ScalarT = int32_t;

int main(int argc, char* argv[])
{
    grb::IndexType       nnodes;
    grb::IndexType       nedges;
    grb::IndexType       data_type = 0;
    grb::IndexArrayType  src_arr;
    grb::IndexArrayType  dst_arr;
    std::vector<ScalarT> weights;

    // read mtx from an input file
    if (argc == 2) {
        if (std::filesystem::path(argv[1]).extension().string() != std::string(".list")) {
            std::cerr << "Unsupported file format" << std::endl;
            return 1;
        }

        std::cout << "Reading input graph from an input file " << argv[1] << std::endl;

        std::string fname(argv[1]);
        std::ifstream in;
        in.open(fname, std::ios::in | std::ios::binary);

        // check open file for write
        if (!in.is_open()) {
            std::cerr << "Error in open file " << fname << std::endl;
            return 1;
        }

        uint32_t uint_tmp;

        in.read((char*) &uint_tmp, sizeof(uint_tmp));
        assert(in);
        nnodes = static_cast<grb::IndexType>(uint_tmp);
        std::cout << "nnodes = " << nnodes << std::endl;

        in.read((char*) &uint_tmp, sizeof(uint_tmp));
        assert(in);
        nedges = static_cast<grb::IndexType>(uint_tmp);
        std::cout << "nedges = " << nedges << std::endl;

        in.read((char*) &uint_tmp, sizeof(uint_tmp));
        assert(in);
        data_type = static_cast<grb::IndexType>(uint_tmp);
        std::cout << "data_type = " << data_type
                  << " (0-binary, 1-integer, 2-float)" << std::endl;

        // src arr
        src_arr.resize(nedges);
        for (size_t i = 0; i < nedges; ++i) {
	          in.read((char*) &uint_tmp, sizeof(uint_tmp));
            assert(in);
            src_arr[i] = static_cast<grb::IndexType>(uint_tmp);
        }

        // dst arr
        dst_arr.resize(nedges);
        for (size_t i = 0; i < nedges; ++i) {
	          in.read((char*) &uint_tmp, sizeof(uint_tmp));
            assert(in);
            dst_arr[i] = static_cast<grb::IndexType>(uint_tmp);
        }

        // close
        in.close();
    } else {
        std::cout << "Using the default matrix" << std::endl;
        nnodes  = num_nodes;
        src_arr = i;
        dst_arr = j;
        nedges  = src_arr.size();
    }

    // set weights to 1s for this BFS_level since actual weights don't matter here
    weights.resize(nedges, 1);

    // TODO Assignment from Initalizer list.
    grb::Matrix<ScalarT> G_karate(nnodes, nnodes);

    G_karate.build(src_arr.begin(), dst_arr.begin(), weights.begin(), src_arr.size());
    std::cout << "Graph: " << std::endl;
    //grb::print_matrix(std::cout, G_karate);

    std::cout << "\n\nRunning bfs_level_masked_v2 ..." << std::endl;
    grb::Vector<grb::IndexType> levels1(nnodes);
    grb::Vector<ScalarT> root(nnodes);
    root.setElement(grb::IndexType(0), 1);
    algorithms::bfs_level_masked_v2(G_karate, root, levels1);
    std::cout << "levels:" << std::endl;
    grb::print_vector(std::cout, levels1);

    //std::cout << "\n\nRunning bfs_level (using mxv) ..." << std::endl;
    //grb::Vector<grb::IndexType> levels(nnodes);
    //algorithms::bfs_level(G_karate, grb::IndexType(0), levels);
    //std::cout << "levels:" << std::endl;
    //grb::print_vector(std::cout, levels);

//    // Trying the row vector approach
//    grb::Matrix<ScalarT>  root(1, nnodes);
//    // pick an arbitrary root:
//    root.setElement(0, 0, 1);
//
//    grb::Matrix<ScalarT> levels1(1, nnodes);
//
//    algorithms::bfs_level(G_karate, root, levels1);
//
//    std::cout << "bfs_level output" << std::endl;
//    std::cout << "root:" << std::endl;
//    grb::print_matrix(std::cout, root);
//    std::cout << "levels (using mxm):" << std::endl;
//    grb::print_matrix(std::cout, levels1);

//    grb::Matrix<ScalarT> levels(1, nnodes);
//    algorithms::batch_bfs_level_masked(G_karate, root, levels);
//
//    std::cout << "Graph: " << std::endl;
//    grb::print_matrix(std::cout, G_karate);
//    std::cout << std::endl;
//    std::cout << "root:" << std::endl;
//    grb::print_matrix(std::cout, root);
//    std::cout << "levels:" << std::endl;
//    grb::print_matrix(std::cout, levels);

    return 0;
}
