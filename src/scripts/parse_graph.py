#!/usr/bin/env python

import argparse
import sys
import os
import struct

graph_dict = {
  # tiny
  'cage5'            : ('directed',   'float'),     #   37, 233   (toy-ish)
  # small
  'wiki-Vote'        : ('directed',   'binary'),    #   8K, 100K  (social)
  ## large
  #'soc-Pokec'        : ('directed',   'binary'),    # 1.6M, 30.0M (social)
  #'belgium_osm'      : ('undirected', 'binary'),    # 1.4M,  1.5M (road)
  #'kron_g500-logn17' : ('undirected', 'integer'),   # 0.1M, 10.0M (kron synthetic)
  #'web-Google'       : ('directed',   'binary'),    # 0.9M,  5.1M (web)
  #'delaunay_n19'     : ('undirected', 'binary'),    # 0.5M,  3.1M (delaunay)
  #'ASIC_680ks'       : ('directed',   'float'),     # 0.7M,  1.7M (circuit sim)
}

def main():

  # args
  parser = argparse.ArgumentParser(description='Parse raw matrix into lil and csr format')
  parser.add_argument('--input', required=True)
  parser.add_argument('--output', required=True)
  args = parser.parse_args()

  # check if all graph files exist
  input_dir  = os.path.abspath(args.input)
  output_dir = os.path.abspath(args.output)
  assert os.path.exists(input_dir),  '{} does not exist'(input_dir)
  assert os.path.exists(output_dir), '{} does not exist'(output_dir)

  for graph in graph_dict.keys():
    graph_file = os.path.join(input_dir, graph, '{}.mtx'.format(graph))
    if not os.path.exists(graph_file):
      print('File not found: {}'.format(graph_file))
      continue

    print('Processing {}'.format(graph))

    direction = graph_dict[graph][0]
    data_type = graph_dict[graph][1]

    #
    # Read input graph file
    #
    with open(graph_file, 'r') as fp:
      lines = fp.readlines()
      count = 0

      edge_list = []

      for line in lines:
        # skip comment lines
        if line.startswith('%'):
          continue

        vals = [ s for s in line.split(' ') ]

        # read header line
        if count == 0:
          assert (len(vals) == 3) and (vals[0] == vals[1]), \
                      'Wrong format: {nnodes} {nnodes} {nedges}'
          nnodes = int(vals[0])
          nedges = int(vals[2])

        # read data lines
        else:
          assert (data_type == 'binary' and len(vals) == 2) or (len(vals) == 3)

          src_id = int(vals[0]) - 1
          dst_id = int(vals[1]) - 1

          # assuming node indices start from 0, substract 1 from each node index
          if   data_type == 'binary':
            val = 1
          elif data_type == 'integer':
            val = int(vals[2])
          elif data_type == 'float':
            val = float(vals[2])
          else:
            assert False, 'Wrong data type {}'.format(data_type)

          edge_list.append((src_id, dst_id, val))
          if direction == 'undirected':
            edge_list.append((dst_id, src_id, val))

        count += 1

    # reset number of edges
    nedges = len(edge_list)

    print('nnodes = {}, nedges = {}, direction = {}, data_type = {}'.format(
                    nnodes, nedges, direction, data_type))

    # sort all edges by rows and then cols
    sorted_edge_list = sorted(edge_list, key = lambda x: (x[0], x[1]))

    ##for e in sorted_edge_list:
    ##  print('({}, {})'.format(e[0], e[1]))

    # dump to a binary file in the following format
    # <nnodes> <nedges> <is_weighted>
    # <list_of_src_node> <list_of_dst_node> <list_of_weight>

    print('Writing to {}.list ...'.format(graph))
    with open(os.path.join(output_dir, '{}.list'.format(graph)), 'wb') as fp:
      fp.write(struct.pack('<I', nnodes))
      fp.write(struct.pack('<I', nedges))

      if   data_type == 'binary':
        fp.write(struct.pack('<I', 0))
      elif data_type == 'integer':
        fp.write(struct.pack('<I', 1))
      elif data_type == 'float':
        fp.write(struct.pack('<I', 2))
      else:
        assert False, 'Wrong data type {}'.format(data_type)

      src_list = [ x[0] for x in sorted_edge_list ]
      for s in src_list:
        fp.write(struct.pack('<I', s))

      dst_list = [ x[1] for x in sorted_edge_list ]
      for d in dst_list:
        fp.write(struct.pack('<I', d))

      weight_list = [ x[2] for x in sorted_edge_list ]
      for w in weight_list:
        if data_type == 'binary' or data_type == 'integer':
          fp.write(struct.pack('<i', w))
        elif data_type == 'float':
          fp.write(struct.pack('<f', w))

    #
    # construct a CSR format for this matrix
    #

    row_ptr_arr = []
    col_idx_arr = dst_list
    val_arr     = weight_list

    cur_row_ptr = 0
    for row in range(nnodes):
      row_ptr_arr.append(cur_row_ptr)

      # skip all values in the same row
      while (cur_row_ptr < nedges) and (src_list[cur_row_ptr] == row):
        cur_row_ptr += 1

    # append the last ptr to the row_ptr_arr
    row_ptr_arr.append(nedges)

    # dump to a binary file in the following format
    # <nnodes> <nedges> <is_weighted>
    # <row_ptr_arr>
    # <col_idx_arr>
    # <val_arr>

    print('Writing to {}.csr ...'.format(graph))
    with open(os.path.join(output_dir, '{}.csr'.format(graph)), 'wb') as fp:
      fp.write(struct.pack('<I', nnodes))
      fp.write(struct.pack('<I', nedges))

      if   data_type == 'binary':
        fp.write(struct.pack('<I', 0))
      elif data_type == 'integer':
        fp.write(struct.pack('<I', 1))
      elif data_type == 'float':
        fp.write(struct.pack('<I', 2))
      else:
        assert False, 'Wrong data type {}'.format(data_type)

      for r in row_ptr_arr:
        fp.write(struct.pack('<I', r))

      for c in col_idx_arr:
        fp.write(struct.pack('<I', c))

      for v in val_arr:
        if data_type == 'binary' or data_type == 'integer':
          fp.write(struct.pack('<i', v))
        elif data_type == 'float':
          fp.write(struct.pack('<f', v))

if __name__ == '__main__':
  main()
