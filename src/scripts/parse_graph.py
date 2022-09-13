#!/usr/bin/env python

import argparse
import sys
import os
import struct

def main():

  # args
  parser = argparse.ArgumentParser(description='Parse raw matrix into lil and csr format')
  parser.add_argument('--input', required=True)
  args = parser.parse_args()

  raw_mtx_file = args.input
  assert os.path.exists(raw_mtx_file), '{} does not exist'(raw_mtx_file)
  raw_mtx_file = os.path.abspath(raw_mtx_file)

  # matrix's name
  mtx_name = raw_mtx_file.split('/')[-1].split('.')[0]
  print('Matrix name: {}'.format(mtx_name))

  # output dir
  out_dir_path = os.path.dirname(raw_mtx_file)
  print('Output directory: {}'.format(out_dir_path))

  with open(raw_mtx_file, 'r') as fp:
    lines = fp.readlines()
    count = 0

    is_weighted = None
    edge_list   = []

    for line in lines:
      # skip comment lines
      if line.startswith('%'):
        continue

      vals = [ int(s) for s in line.split(' ') ]

      # read header line
      if count == 0:
        assert (len(vals) == 3) and (vals[0] == vals[1]), \
                    'Wrong format: {nnodes} {nnodes} {nedges}'
        nnodes = vals[0]
        nedges = vals[2]
      # read data lines
      else:
        assert (len(vals) == 2) or (len(vals) == 3)
        if is_weighted is None:
          is_weighted = True if len(vals) == 3 else False

        # assuming node indices start from 0, substract 1 from each node index
        if is_weighted:
          edge_list.append((vals[0] - 1, vals[1] - 1, vals[2]))
        else:
          edge_list.append((vals[0] - 1, vals[1] - 1, 1))

      count += 1

  assert nedges == (count - 1), "nedges={} count={}".format(nedges, count)

  print('nnodes = {}, nedges = {}, is_weighted = {}'.format(nnodes, nedges, is_weighted))

  # sort all edges by rows and then cols
  sorted_edge_list = sorted(edge_list, key = lambda x: (x[0], x[1]))

  #for e in sorted_edge_list:
  #  print('({}, {})'.format(e[0], e[1]))

  # dump to a binary file in the following format
  # <nnodes> <nedges> <is_weighted>
  # <list_of_src_node> <list_of_dst_node> <list_of_weight>

  with open(out_dir_path + '/' + mtx_name + '.list', 'wb') as fp:
    fp.write(struct.pack('<I', nnodes))
    fp.write(struct.pack('<I', nedges))
    fp.write(struct.pack('<I', 1 if is_weighted else 0))

    src_list = [ x[0] for x in sorted_edge_list ]
    for s in src_list:
      fp.write(struct.pack('<I', s))

    dst_list = [ x[1] for x in sorted_edge_list ]
    for d in dst_list:
      fp.write(struct.pack('<I', d))

    weight_list = [ x[2] for x in sorted_edge_list ]
    for w in weight_list:
      fp.write(struct.pack('<i', w))

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

  with open(out_dir_path + '/' + mtx_name + '.csr', 'wb') as fp:
    fp.write(struct.pack('<I', nnodes))
    fp.write(struct.pack('<I', nedges))
    fp.write(struct.pack('<I', 1 if is_weighted else 0))

    for r in row_ptr_arr:
      fp.write(struct.pack('<I', r))

    for c in col_idx_arr:
      fp.write(struct.pack('<I', c))

    for v in val_arr:
      fp.write(struct.pack('<i', v))

if __name__ == '__main__':
  main()
