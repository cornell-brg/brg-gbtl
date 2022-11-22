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
// bfs_level_demo.cpp
//=========================================================================
// Doing BFS-level
// Input graphs are expected in CSR or transposed CSR format

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
    std::cerr << "Wrong command line arguments <inp_graph> <root_node> <ref_graph>" << std::endl;
    return 1;
  }

  bool is_transposed_csr = false;
  if (std::filesystem::path(argv[1]).extension().string() != std::string(".csr")) {
    is_transposed_csr = false;
  } else if (std::filesystem::path(argv[1]).extension().string() != std::string(".tcsr")) {
    is_transposed_csr = true;
  } else {
    std::cerr << "Unsupported file format: " << argv[1] << std::endl;
    return 1;
  }

  grb::IndexType nnodes;
  grb::IndexType nedges;
  grb::IndexType data_type;
  grb::IndexType root_node = 0;

  // process root node
  try {
    root_node = std::stoi(std::string(argv[2]));
  } catch (...) {
    std::cerr << "Invalid root node: " << argv[2] << std::endl;
  }

  // process input graph's info
  std::string fname(argv[1]);
  std::ifstream in_file;
  in_file.open(fname, std::ios::in | std::ios::binary);

  if (!in_file.is_open()) {
    std::cerr << "Error in opening file " << fname << std::endl;
    return 1;
  }

  // read nnodes, nedges, data type
  grb::IndexType uint_tmp;
  in_file.read((char*) &uint_tmp, sizeof(grb::IndexType));
  assert(in_file);
  nnodes = static_cast<grb::IndexType>(uint_tmp);
  std::cout << "nnodes = " << nnodes << std::endl;

  in_file.read((char*) &uint_tmp, sizeof(grb::IndexType));
  assert(in_file);
  nedges = static_cast<grb::IndexType>(uint_tmp);
  std::cout << "nedges = " << nedges << std::endl;

  in_file.read((char*) &uint_tmp, sizeof(grb::IndexType));
  assert(in_file);
  data_type = static_cast<grb::IndexType>(uint_tmp);
  std::cout << "data_type = " << data_type
            << " (0-binary, 1-integer, 2-float)" << std::endl;

  if (root_node >= nnodes) {
      std::cerr << "Invalid root node: " << root_node << std::endl;
      return 1;
  }

  std::cout << "root_node = " << root_node << std::endl;

  // read row_ptr
  grb::IndexArrayType row_ptr_arr(nnodes + 1);
  for (grb::IndexType i = 0; i < nnodes + 1; ++i) {
    in_file.read((char*) &uint_tmp, sizeof(grb::IndexType));
    assert(in_file);
    row_ptr_arr[i] = static_cast<grb::IndexType>(uint_tmp);
  }

  // read col_idx
  grb::IndexArrayType col_idx_arr(nedges);
  for (grb::IndexType i = 0; i < nedges; ++i) {
    in_file.read((char*) &uint_tmp, sizeof(grb::IndexType));
    assert(in_file);
    col_idx_arr[i] = static_cast<grb::IndexType>(uint_tmp);
  }

  // read mtx_dat
  // for now, just all 1s
  grb::IndexArrayType mtx_dat_arr(nedges, 1);

  // close the file
  in_file.close();

  // create a matrix
  grb::Matrix<ScalarT> G_karate(nnodes, nnodes);
  G_karate.build_from_csr(row_ptr_arr, col_idx_arr, mtx_dat_arr, is_transposed_csr);

  // print the matrix
  //grb::print_matrix(std::cout, G_karate);

  grb::Vector<grb::IndexType> levels(nnodes);
  grb::Vector<ScalarT> root(nnodes);
  root.setElement(root_node, 1);

  ///** Switch to detailed CPU */
  //gem5::switch_cpus(true);

  ///** Turn on vector engine */
  //gem5::vstart();

  algorithms::bfs_level_masked_v2(G_karate, root, levels);

  ///** Turn off vector engine */
  //gem5::vend();

  // check the output
  //std::cout << "levels:" << std::endl;
  //grb::print_vector(std::cout, levels);
  std::ifstream ref_file;
  ref_file.open(argv[3]);

  if (!ref_file.is_open()) {
    std::cerr << "Error in opening file " << argv[3] << std::endl;
    return 1;
  }

  ref_file.read((char*) &uint_tmp, sizeof(grb::IndexType));
  assert(ref_file);
  auto ref_nnodes = static_cast<grb::IndexType>(uint_tmp);
  if (ref_nnodes != nnodes) {
    std::cerr << "[FAILED] nnodes != ref_nnodes: " << nnodes << " != " << ref_nnodes << std::endl;
    return 1;
  }

  for (auto i = 0; i < nnodes; ++i) {
    ref_file.read((char*) &uint_tmp, sizeof(uint_tmp));
    assert(ref_file);
    ScalarT ref_val = static_cast<ScalarT>(uint_tmp);
    ScalarT val;
    try {
      val = levels.extractElement(i);
    } catch (...) {
      val = 0;
    }
    if (ref_val != val) {
      std::cerr << "[FAILED] i = " << i << ", ref_val != val: " << ref_val << " != " << val << std::endl;
      return 1;
    }
  }

  std::cout << "[passed]" << std::endl;
  return 0;
}
