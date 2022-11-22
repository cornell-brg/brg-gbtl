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

//=========================================================================
// bfs_level_ref_demo.cpp
//=========================================================================
// Doing BFS-level and generate output in a binary file

#include <iostream>
#include <filesystem>
#include <fstream>
#include <cassert>

#include <graphblas/graphblas.hpp>
#include <algorithms/bfs.hpp>
#include <graphblas/gem5_helpers.hpp>

using ScalarT = uint32_t;

int main(int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Wrong command line arguments <inp_graph> <root_node> <out_graph>" << std::endl;
    return 1;
  }

  if (std::filesystem::path(argv[1]).extension().string() != std::string(".list")) {
    std::cerr << "Unsupported file format: " << argv[1] << std::endl;
    return 1;
  }

  grb::IndexType       nnodes;
  grb::IndexType       nedges;
  grb::IndexType       data_type;
  grb::IndexType       root_node;
  grb::IndexArrayType  src_arr;
  grb::IndexArrayType  dst_arr;
  std::vector<ScalarT> weights;

  try {
    root_node = std::stoi(std::string(argv[2]));
  } catch (...) {
    std::cerr << "Invalid root node ID: " << argv[2] << std::endl;
  }

  std::cout << "Reading input graph from " << argv[1] << std::endl;
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

  if (root_node >= nnodes) {
    std::cerr << "Invalid root node: " << root_node << std::endl;
    return 1;
  }

  std::cout << "root_node = " << root_node << std::endl;

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

  // set weights to 1s for this BFS_level since actual weights don't matter here
  weights.resize(nedges, 1);

  // TODO Assignment from Initalizer list.
  grb::Matrix<ScalarT> G_karate(nnodes, nnodes);
  G_karate.build(src_arr.begin(), dst_arr.begin(), weights.begin(), src_arr.size());

  // print the matrix
  //grb::print_matrix(std::cout, G_karate);

  grb::Vector<ScalarT> root(nnodes);
  root.setElement(root_node, 1);
  grb::Vector<grb::IndexType> levels(nnodes);

  // do BFS-level (masked)
  algorithms::bfs_level_masked(G_karate, root, levels);

  // dump the output to a file
  //std::cout << "levels:" << std::endl;
  //grb::print_vector(std::cout, levels);

  std::ofstream out(argv[3], std::ios::out | std::ios::binary);
  uint_tmp = nnodes;
  out.write((char *) &uint_tmp, sizeof(uint_tmp));
  for (auto i = 0; i < nnodes; ++i) {
    try {
      uint_tmp = levels.extractElement(i);
    } catch (...) {
      uint_tmp = 0;
    }
    out.write((char *) &uint_tmp, sizeof(uint_tmp));
  }

  return 0;
}
