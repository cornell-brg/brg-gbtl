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

#pragma once

#include <iostream>
#include <vector>
#include <typeinfo>
#include <stdexcept>
#include <algorithm>

#include <graphblas/graphblas.hpp>

//****************************************************************************

namespace grb
{
    namespace backend
    {

        template<typename ScalarT>
        class CsrSparseMatrix
        {
        public:
            using ScalarType = ScalarT;
            using ElementType = std::tuple<IndexType, ScalarT>;
            using RowType = std::vector<ElementType>;

            using ScalarArrayType     = std::vector<ScalarType>;
            using ScalarArrayIterator = typename std::vector<ScalarType>::iterator;
            using IndexArrayIterator  = typename std::vector<IndexType>::iterator;

            // Constructor
            CsrSparseMatrix(IndexType num_rows,
                            IndexType num_cols)
                : m_num_rows(num_rows),
                  m_num_cols(num_cols),
                  m_nvals(0),
                  m_is_transposed(false)
            {
            }

            // Constructor - copy
            CsrSparseMatrix(CsrSparseMatrix<ScalarT> const &rhs)
                : m_num_rows(rhs.m_num_rows),
                  m_num_cols(rhs.m_num_cols),
                  m_nvals(rhs.m_nvals),
                  m_is_transposed(false)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::CsrSparseMatrix() (copy) INTERNAL ERROR");
            }

            // Constructor - dense from dense matrix
            CsrSparseMatrix(std::vector<std::vector<ScalarT>> const &val)
                : m_num_rows(val.size()),
                  m_num_cols(val[0].size()),
                  m_is_transposed(false)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::CsrSparseMatrix() (dense from dense) INTERNAL ERROR");
            }

            // Constructor - sparse from dense matrix, removing specifed implied zeros
            CsrSparseMatrix(std::vector<std::vector<ScalarT>> const &val,
                            ScalarT zero)
                : m_num_rows(val.size()),
                  m_num_cols(val[0].size()),
                  m_is_transposed(false)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::CsrSparseMatrix() (sparse from dense) INTERNAL ERROR");
            }

            // Destructor
            ~CsrSparseMatrix()
            {}

            // Assignment (currently restricted to same dimensions)
            CsrSparseMatrix<ScalarT> &operator=(CsrSparseMatrix<ScalarT> const &rhs)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::operator=() INTERNAL ERROR");
                return *this;
            }

            // EQUALITY OPERATORS
            /**
             * @brief Equality testing for LilMatrix.
             * @param rhs The right hand side of the equality operation.
             * @return If this LilMatrix and rhs are identical.
             */
            bool operator==(CsrSparseMatrix<ScalarT> const &rhs) const
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::operator==() INTERNAL ERROR");
                return false;
            }

            /**
             * @brief Inequality testing for LilMatrix.
             * @param rhs The right hand side of the inequality operation.
             * @return If this LilMatrix and rhs are not identical.
             */
            bool operator!=(CsrSparseMatrix<ScalarT> const &rhs) const
            {
                return !(*this == rhs);
            }

            template<typename RAIteratorI,
                     typename RAIteratorJ,
                     typename RAIteratorV,
                     typename DupT>
            void build(RAIteratorI  i_it,
                       RAIteratorJ  j_it,
                       RAIteratorV  v_it,
                       IndexType    n,
                       DupT         dup)
            {
                if (n == 0) {
                    // no non-zero value, just return
                    return;
                }

                // sanity check if row indices are sorted
                if (!std::is_sorted(i_it, i_it + n)) {
                    throw grb::InvalidIndexException(
                          "CsrSparseMatrix::build() INTERNAL ERROR");
                }

                m_col_idx_arr.assign(j_it, j_it + n);
                m_mtx_dat_arr.assign(v_it, v_it + n);

                m_row_ptr_arr.resize(m_num_rows + 1);
                IndexType cur_row_ptr = 0;
                for (IndexType r = 0; r < m_num_rows; ++r) {
                    m_row_ptr_arr[r] = cur_row_ptr;

                    // skip all values in the same row
                    // sanity check: column indices must be increasing
                    IndexType col_idx;
                    uint64_t count = 0;

                    while ((cur_row_ptr < n) && (*(i_it + cur_row_ptr) == r)) {
                        if (count > 0) {
                            if (m_col_idx_arr[cur_row_ptr] <= col_idx) {
                                throw grb::InvalidIndexException(
                                      "CsrSparseMatrix::build() INTERNAL ERROR");
                            }
                        }

                        col_idx = m_col_idx_arr[cur_row_ptr];
                        count++;
                        cur_row_ptr++;
                    }
                }

                // set number of non-zeros
                m_nvals = m_mtx_dat_arr.size();
                m_row_ptr_arr[m_num_rows] = m_nvals;
            }

            template<typename RowPtrIterator,
                     typename ColIdxIterator,
                     typename MtxDatIterator>
            void build_from_csr(RowPtrIterator row_ptr_it,
                                ColIdxIterator col_idx_it,
                                MtxDatIterator mtx_dat_it,
                                IndexType      nnodes,
                                IndexType      nedges,
                                bool           is_transposed)
            {
                if (nedges == 0) {
                    // no non-zero value, just return
                    return;
                }

                m_row_ptr_arr.assign(row_ptr_it, row_ptr_it + (nnodes + 1));
                m_col_idx_arr.assign(col_idx_it, col_idx_it + nedges);
                m_mtx_dat_arr.assign(mtx_dat_it, mtx_dat_it + nedges);

                m_is_transposed = is_transposed;
                m_nvals         = nedges;
            }

            void clear()
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::clear() INTERNAL ERROR");
            }

            IndexType nrows() const { return m_num_rows; }
            IndexType ncols() const { return m_num_cols; }
            IndexType nvals() const { return m_nvals; }

            /**
             * @brief Resize the matrix dimensions (smaller or larger)
             *
             * @param[in]  new_num_rows  New number of rows (zero is invalid)
             * @param[in]  new_num_cols  New number of columns (zero is invalid)
             *
             */
            void resize(IndexType new_num_rows, IndexType new_num_cols)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::resize() INTERNAL ERROR");
            }

            bool hasElement(IndexType irow, IndexType icol) const
            {
                if (irow >= m_num_rows || icol >= m_num_cols)
                {
                    throw IndexOutOfBoundsException(
                        "get_value_at: index out of bounds");
                }

                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::hasElement() INTERNAL ERROR");

                return false;
            }

            // Get value at index
            ScalarT extractElement(IndexType irow, IndexType icol) const
            {
                if (irow >= m_num_rows || icol >= m_num_cols)
                {
                    throw IndexOutOfBoundsException(
                        "extractElement: index out of bounds");
                }

                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::extractElement() INTERNAL ERROR");
            }

            // Set value at index
            void setElement(IndexType irow, IndexType icol, ScalarT const &val)
            {
                if (irow >= m_num_rows || icol >= m_num_cols)
                {
                    throw IndexOutOfBoundsException("setElement: index out of bounds");
                }

                // FIXME @Tuan: check if the matrix has not been built yet!
                //IndexType col_start = m_row_ptr_arr[irow];
                //IndexType col_end   = (irow < m_num_rows - 1) ?
                //                          m_row_ptr_arr[irow + 1] : m_nvals;

                //auto i = col_start;
                //for (; i < col_end; ++i) {
                //    printf("%d %d", i, m_col_idx_arr[i]);
                //    if (icol == m_col_idx_arr[i]) {
                //        // found an existing element, just update it
                //        m_mtx_dat_arr[i] = val;
                //        return;
                //    } else if (icol < m_col_idx_arr[i]) {
                //        // found a place to insert the new element
                //        break;
                //    }
                //}
                //printf("\n");

                //// insert the value right before the i-th element
                //m_col_idx_arr.insert(m_col_idx_arr.begin() + i, icol);
                //m_mtx_dat_arr.insert(m_mtx_dat_arr.begin() + i, val);
                //m_nvals++;

                //// shift row_ptr after irow by 1
                //for (i = irow + 1; i < m_num_rows; ++i)
                //    m_row_ptr_arr[i]++;

                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::setElement() INTERNAL ERROR");
            }

            // Set value at index + 'merge' with any existing value
            // according to the BinaryOp passed.
            template <typename BinaryOpT>
            void setElement(IndexType irow, IndexType icol, ScalarT const &val,
                            BinaryOpT merge)
            {
                if (irow >= m_num_rows || icol >= m_num_cols)
                {
                    throw IndexOutOfBoundsException(
                        "setElement(merge): index out of bounds");
                }

                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::setElement() INTERNAL ERROR");
            }

            void removeElement(IndexType irow, IndexType icol)
            {
                if (irow >= m_num_rows || icol >= m_num_cols)
                {
                    throw IndexOutOfBoundsException("removeElement: index out of bounds");
                }

                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::removeElement() INTERNAL ERROR");
            }

            void recomputeNvals()
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::recomputeNvals() INTERNAL ERROR");
            }

            // TODO: add error checking on dimensions?
            void swap(CsrSparseMatrix<ScalarT> &rhs)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::swap() INTERNAL ERROR");
            }

            // Row access
            // Warning if you use this non-const row accessor then you should
            // call recomputeNvals() at some point to fix it
            RowType &operator[](IndexType row_index)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::operator[]() INTERNAL ERROR");
                return m_dummy_row;
            }

            RowType const &operator[](IndexType row_index) const
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::operator[]() INTERNAL ERROR");
                return m_dummy_row;
            }

            IndexType getNvals(IndexType row_idx) const
            {
                if (row_idx >= m_num_rows)
                {
                    throw IndexOutOfBoundsException("getNvals: index out of bounds");
                }

                return m_row_ptr_arr[row_idx + 1] - m_row_ptr_arr[row_idx];
            }

            const IndexType* getColIdxArr(IndexType row_idx) const
            {
                if (row_idx >= m_num_rows)
                {
                    throw IndexOutOfBoundsException("getNvals: index out of bounds");
                }

                return m_col_idx_arr.data() + m_row_ptr_arr[row_idx];
            }

            const ScalarType* getValArr(IndexType row_idx) const
            {
                if (row_idx >= m_num_rows)
                {
                    throw IndexOutOfBoundsException("getNvals: index out of bounds");
                }

                return m_mtx_dat_arr.data() + m_row_ptr_arr[row_idx];
            }

            // Allow casting
            template <typename OtherScalarT>
            void setRow(
                IndexType row_index,
                std::vector<std::tuple<IndexType, OtherScalarT> > const &row_data)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::setRow() INTERNAL ERROR");
            }

            // When not casting vector swap used...should we use move semantics?
            void setRow(
                IndexType row_index,
                std::vector<std::tuple<IndexType, ScalarT> > &&row_data)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::setRow() INTERNAL ERROR");
            }


            // Allow casting. TODO Do we need one that does not need casting?
            // mergeRow with no accumulator is same as setRow
            template <typename OtherScalarT, typename AccumT>
            void mergeRow(
                IndexType row_index,
                std::vector<std::tuple<IndexType, OtherScalarT> > &row_data,
                NoAccumulate const &op)
            {
                setRow(row_index, row_data);
            }


            // Allow casting. TODO Do we need one that does not need casting?
            template <typename OtherScalarT, typename AccumT>
            void mergeRow(
                IndexType row_index,
                std::vector<std::tuple<IndexType, OtherScalarT> > &row_data,
                AccumT const &op)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::mergeRow() INTERNAL ERROR");
            }

            /// @deprecated Only needed for 4.3.7.3 assign: column variant"
            /// @todo need move semantics.
            using ColType = std::vector<std::tuple<IndexType, ScalarT> >;
            ColType getCol(IndexType col_index) const
            {
                std::vector<std::tuple<IndexType, ScalarT> > data;

                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::getCol() INTERNAL ERROR");

                return data;  // hopefully compiles to a move
            }

            /// @deprecated Only needed for 4.3.7.3 assign: column variant"
            /// @note col_data must be in increasing index order
            /// @todo this could be vastly improved.
            template <typename OtherScalarT>
            void setCol(
                IndexType col_index,
                std::vector<std::tuple<IndexType, OtherScalarT> > const &col_data)
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::setCol() INTERNAL ERROR");
            }

            template<typename RAIteratorIT,
                     typename RAIteratorJT,
                     typename RAIteratorVT>
            void extractTuples(RAIteratorIT        row_it,
                               RAIteratorJT        col_it,
                               RAIteratorVT        values) const
            {
                /** TODO */
                throw grb::NotImplementedException(
                      "CsrSparseMatrix::extractTuples() INTERNAL ERROR");
            }

            // output specific to the storage layout of this type of matrix
            void printInfo(std::ostream &os) const
            {
                os << "Optimized RVV Serial Backend: ";
                os << "backend::CsrSparseMatrix<" << typeid(ScalarT).name() << "> ";
                os << "(" << m_num_rows << " x " << m_num_cols << "), nvals = "
                   << nvals() << std::endl;

                if (m_num_rows == 0 && m_num_cols == 0)
                    return;

                // Used to print data in storage format instead of like a matrix
                #ifdef GRB_MATRIX_PRINT_RAW_STORAGE
                    os << "row_ptr: ";
                    for (auto r : m_row_ptr_arr) {
                        os << r << " ";
                    }
                    os << std::endl;

                    os << "col_idx: ";
                    for (auto c : m_col_idx_arr) {
                        os << c << " ";
                    }
                    os << std::endl;

                    os << "values : ";
                    for (auto v : m_mtx_dat_arr) {
                        os << v << " ";
                    }
                    os << std::endl;
                #else
                    IndexType row_ptr_idx = 0;
                    IndexType val_idx     = 0;

                    for (IndexType row_idx = 0; row_idx < m_num_rows; ++row_idx)
                    {
                        // We like to start with a little whitespace indent
                        os << ((row_idx == 0) ? "  [[" : "   [");

                        IndexType min_col_idx = m_row_ptr_arr[row_idx];
                        IndexType max_col_idx = m_row_ptr_arr[row_idx + 1];

                        for (IndexType col_idx = 0; col_idx < m_num_cols; ++col_idx) {
                            if (min_col_idx < max_col_idx &&
                                col_idx == m_col_idx_arr[min_col_idx])
                            {
                                // non-zero element
                                os << ((col_idx > 0) ? ", " : " ");
                                os << m_mtx_dat_arr[val_idx];
                                min_col_idx++;
                            } else {
                                // zero element
                                os << ((col_idx > 0) ? ", " : " ");
                                os << " ";
                            }
                        }

                        os << ((row_idx == m_num_rows - 1 ) ? "]]" : "]\n");
                    }
                #endif
            }

            friend std::ostream &operator<<(std::ostream             &os,
                                            CsrSparseMatrix<ScalarT> const &mat)
            {
                mat.printInfo(os);
                return os;
            }

        private:
            IndexType m_num_rows;
            IndexType m_num_cols;
            IndexType m_nvals;

            // FIXME: dumy RowType
            RowType m_dummy_row;

            IndexArrayType  m_col_idx_arr;
            IndexArrayType  m_row_ptr_arr;
            ScalarArrayType m_mtx_dat_arr;

            bool            m_is_transposed;
        };

    } // namespace backend

} // namespace grb
