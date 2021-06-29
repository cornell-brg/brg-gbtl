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
#include <numeric>
#include <map>

#include <graphblas/graphblas.hpp>
#include <graphblas/indices.hpp>

//****************************************************************************

namespace grb
{
    namespace backend
    {

        /**
         * @brief Class representing a sparse vector using an index vector
         * and a values vector. Based on frontier class in GKC 1.0
         */

        // Sort or no sort? 
        // When to keep things synchronized?
        // How to keep updates?

        template<typename ScalarT>
        class GKCSparseVector
        {
        public:
            using ScalarType = ScalarT;

            // Constructor with some default values
            /**
             * @brief construct an empty sparse vector with space for 
             * num_vals values.
             * The sparse vector can be structure only, with weights set
             * to false, so only vertex IDs are stored.
             */
            GKCSparseVector(IndexType num_vals,
                      bool weighted = true)
                : m_num_vals(num_vals),
                  m_num_stored_vals(0),
                  m_weighted(weighted)
            {
                if (m_num_vals <= 0){
                    throw InvalidValueException();
                }
                m_indices.resize(m_num_vals);
                if (m_weighted)
                    m_weights.resize(m_num_vals);
            }

            // Create GKC Sparse Vector with default weight/value
            /* THIS IS NOT A SENSIBLE USE OF THIS VECTOR. A DENSE VECTOR SHOULD INSTEAD BE USED. */
            GKCSparseVector(IndexType num_vals, ScalarT const &value)
                : m_num_vals(num_vals),
                  m_num_stored_vals(num_vals),
                  m_weighted(true),
                  m_weights(num_vals, value),
                  m_indices(num_vals)
            {
                if (m_num_vals <= 0){
                    throw InvalidValueException();
                }
                // Fill with ascending indices
                for (auto idx = 0; idx < m_num_vals; idx++){
                    m_indices[idx] = idx;
                }
                // Is iota as good (and parallelizable) as the above loop?
                //std::iota(m_indices.begin(), m_indices.end(), 0);
            }

            // Constructor - copy
            GKCSparseVector(GKCSparseVector<ScalarT> const &rhs)
                : m_num_vals(rhs.m_num_vals),
                  m_weighted(rhs.m_weighted),
                  m_indices(rhs.m_indices),
                  m_weights(rhs.m_weights),
                  m_num_stored_vals(rhs.m_num_stored_vals)
            {
            }

            // Should the non-filtering dense constructor be supplied?
            // It makes little sense to create a dense vector with a 
            // Sparse vector class.
            // Constructor - dense from dense vector
            GKCSparseVector(std::vector<ScalarT> const &val)
            {
                m_num_vals = val.size();
                // Dense has both structure and weights.
                m_weighted = true;
                m_indices.resize(m_num_vals);
                m_weights.resize(m_num_vals);
                m_num_stored_vals = 0;
                
                /// @todo: parallelize?
                // Copy values from dense vector into sparse
                for (auto idx = 0; idx < m_num_vals; idx++){
                    m_indices[m_num_stored_vals] = idx;
                    m_weights[m_num_stored_vals] = val[idx];
                    m_num_stored_vals ++;
                }
            }

            // Constructor - sparse from dense vector, removing specifed implied zeros
            GKCSparseVector(std::vector<ScalarT> const &val,
                            ScalarT zero)
            {
                m_num_vals = val.size();
                // Dense has both structure and weights.
                m_weighted = true;
                m_indices.resize(m_num_vals);
                m_weights.resize(m_num_vals);
                m_num_stored_vals = 0;
                
                /// @todo: parallelize?
                // Copy values from dense vector into sparse
                // While ignoring 'zero' elements
                for (auto idx = 0; idx < m_num_vals; idx++){
                    if (val[idx] != zero)
                    {
                        m_indices[m_num_stored_vals] = idx;
                        m_weights[m_num_stored_vals] = val[idx];
                        m_num_stored_vals ++;
                    }
                }
            }

            // Build Constructor - parse sparse coordinate data and construct vector.
            // Similar to build method, but baked into the constructor.
            // Does not use addElement.
            /// @todo: how do we limit the size of this vector? Right now it's 
            // implicit in the data provided by i_it.
            template <typename RAIteratorI,
                      typename RAIteratorV,
                      typename BinaryOpT = grb::Second<ScalarType> >
            GKCSparseVector(RAIteratorI i_it,
                            RAIteratorV v_it,
                            IndexType n,
                            BinaryOpT dup = BinaryOpT())
            {
                /// @todo require random access iterators
                // scan the data to determine num_vals
                /// @todo: OMP max reduction
                IndexType max_idx = 0;
                for (size_t i = 0; i < n; i++){
                    max_idx = std::max(*(i_it+i), max_idx);
                }

                /// @todo: should we always allocate to the max number of 
                /// vertices? Or in this case just to n that we know we may need?
                m_num_vals = max_idx+1;
                m_num_stored_vals = 0;

                // allocate memory
                m_indices.resize(m_num_vals);
                /// @todo how to detect if graph is weighted?
                m_weighted = true;
                m_weights.resize(m_num_vals);

                // Copy data from iterators into vector
                std::map<IndexType, IndexType> already_inserted;
                for (IndexType idx = 0; idx < n; idx++)
                {
                    IndexType vidx = *(i_it+idx);
                    auto found_itr = already_inserted.find(vidx);
                    if (found_itr == already_inserted.end())
                    {
                        // Index not recognized, add new entry
                        already_inserted[vidx] = m_num_stored_vals;
                        m_indices[m_num_stored_vals] = vidx;
                        m_weights[m_num_stored_vals] = *(v_it + idx);
                        m_num_stored_vals++;
                    }
                    else  
                    {
                        // Already have a value, so merge the weight.
                        auto old_idx = found_itr->second;
                        m_weights[old_idx] = dup(m_weights[old_idx], *(v_it + idx));
                    }
                }
            }

            // Destructor
            ~GKCSparseVector()
            {}

            // Copy assignment (currently restricted to same dimensions)
            GKCSparseVector<ScalarT> &operator=(GKCSparseVector<ScalarT> const &rhs)
            {
                if (this != &rhs) 
                {
                    if (m_num_vals != rhs.m_num_vals)
                    {
                        throw DimensionException("Dimensions of vectors do not match.");
                    }    
                    m_num_vals = rhs.m_num_vals;
                    m_num_stored_vals = rhs.m_num_stored_vals;
                    m_weighted = rhs.m_weighted;
                    m_indices = rhs.m_indices;
                    m_weights = rhs.m_weights;
                }
                return *this;
            }

            // EQUALITY OPERATORS
            /**
             * @brief Equality testing for GKC Vector.
             * @param rhs The right hand side of the equality operation.
             * @return If this GKC Vector and rhs are identical.
             */
            bool operator==(GKCSparseVector<ScalarT> const &rhs) const
            {
                return ((m_num_vals == rhs.m_num_vals) &&
                        (m_weighted == rhs.m_weighted) &&
                        (m_num_stored_vals == rhs.m_num_stored_vals) &&
                        (m_indices == rhs.m_indices) && 
                        (m_weights == rhs.m_weights));
            }

            /**
             * @brief Inequality testing for GKC Vector.
             * @param rhs The right hand side of the inequality operation.
             * @return If this GKC Vector and rhs are not identical.
             */
            
            bool operator!=(GKCSparseVector<ScalarT> const &rhs) const
            {
                return !(*this == rhs);
            }
            
            template<typename RAIteratorIT,
                     typename RAIteratorVT,
                     typename BinaryOpT = grb::Second<ScalarType> >
            void build(RAIteratorIT  i_it,
                       RAIteratorVT  v_it,
                       IndexType     n,
                       BinaryOpT     dup = BinaryOpT())
            {
                /// @todo require random access iterators
                // scan the data to determine num_edges, num_rows, num_cols
                m_num_vals = n;
                m_num_stored_vals = 0;

                // allocate memory
                m_indices.resize(m_num_vals);
                /// @todo how to detect if graph is weighted?
                m_weighted = true;
                m_weights.resize(m_num_vals);

                // Copy data from iterators into vector
                std::map<IndexType, IndexType> already_inserted;
                for (IndexType idx = 0; idx < n; idx++)
                {
                    IndexType vidx = *(i_it+idx);
                    if (already_inserted[vidx] < idx + 1){
                        // Already have a value, so merge the weight.
                        auto old_idx = already_inserted[vidx] - 1;
                        m_weights[old_idx] = dup(m_weights[old_idx], *(v_it + idx));
                    }
                    else  // Index not recognized, add new entry
                    {
                        already_inserted[vidx] = idx + 1;
                        m_indices[idx] = vidx;
                        m_weights[idx] = *(v_it + idx);
                        m_num_stored_vals++;
                    }
                }
            }
            
            void clear()
            {
                /// @todo make atomic? transactional?
                m_num_stored_vals = 0;
                /// @todo clear or resize to 0?
                m_indices.clear();
                m_weights.clear();
            }

            IndexType size() const { return m_num_vals; }
            IndexType nvals() const { return m_num_stored_vals; }

            /**
             * @brief Resize the vector dimensions (smaller or larger)
             *
             * @param[in]  new_num_rows  New number of rows (zero is invalid)
             * @param[in]  new_num_cols  New number of columns (zero is invalid)
             *
             */
            void resize(IndexType new_size)
            {
                // Invalid values check by frontend
                if (new_size == 0)
                    throw DimensionException("Cannot resize vector to have zero dimension.");
                if (new_size < m_num_vals)
                {
                    m_num_vals = new_size;
                    // Linear scan and shift
                    IndexType new_num_stored = 0;
                    for (size_t i = 0; i < m_num_stored_vals; i++)
                    {
                        size_t idx = m_indices[i];
                        if (idx < new_size)
                        {
                            m_indices[new_num_stored] = idx;
                            if (m_weighted) {
                                m_weights[new_num_stored] = m_weights[i];
                            }
                            new_num_stored++;
                        }
                    }
                    m_num_stored_vals = new_num_stored;
                }
                else // Increase size of vectors
                {
                    m_num_vals = new_size;
                    m_indices.resize(new_size);
                    if (m_weighted) {
                        m_weights.resize(new_size);
                    }
                } 
            }

            bool hasElement(IndexType index) const
            {
                if (index >= m_num_vals)
                {
                    throw IndexOutOfBoundsException();
                }
                for (size_t idx = 0; idx < m_num_stored_vals; idx++)
                {
                    auto vidx = m_indices[idx];
                    if (vidx == index) return true;
                }
                return false;
            }

            ScalarT extractElement(IndexType index) const
            {
                /// @todo mark: if vector is sorted, 
                /// use binary search.
                /// Need to add a 'sorted' flag that is reset if 
                // vector is modified.
                if (index >= m_num_vals)
                {
                    throw IndexOutOfBoundsException();
                }
                for (size_t idx = 0; idx < m_num_stored_vals; idx++)
                {
                    auto vidx = m_indices[idx];
                    if (vidx == index)
                    {
                        if (m_weighted)
                        {
                            return m_weights[idx];
                        }
                        else // What to do if no weights?
                        {
                            return (ScalarT)1;
                        }
                    }
                }
                // No value found; throw error
                throw NoValueException();
            }

            void setElement(IndexType index, ScalarT const &new_val)
            {
                if (index >= m_num_vals)
                {
                    throw IndexOutOfBoundsException();
                }
                for (size_t idx = 0; idx < m_num_stored_vals; idx++)
                {
                    auto vidx = m_indices[idx];
                    if (vidx == index)
                    {
                        if (m_weighted){
                            m_weights[idx] = new_val;
                        }
                        return;
                    }
                }
                // No value found; insert it:
                m_indices[m_num_stored_vals] = index;
                if (m_weighted){
                    m_weights[m_num_stored_vals] = new_val;
                }
                m_num_stored_vals++;
            }

            void removeElement(IndexType index)
            {
                if (index >= m_num_vals)
                {
                    throw IndexOutOfBoundsException();
                }
                // Step 1: find element
                if (index > m_num_vals)
                {
                    throw IndexOutOfBoundsException();
                }
                for (size_t idx = 0; idx < m_num_stored_vals; idx++)
                {
                    auto vidx = m_indices[idx];
                    if (vidx == index)
                    {
                        // Step 2: vector doesn't need to remain sorted, 
                        // so just replace with last element.
                        // NOT THREAD SAFE!
                        if (idx < m_num_stored_vals - 1){
                            m_indices[idx] = m_indices[m_num_stored_vals - 1];
                            if (m_weighted)
                            {
                                m_weights[idx] = m_weights[m_num_stored_vals-1];
                            }
                            m_num_stored_vals--;
                        }
                        return;
                    }
                }
                // No value found; throw error
                throw NoValueException();

                
            }

            template<typename RAIteratorIT,
                     typename RAIteratorVT>
            void extractTuples(RAIteratorIT        i_it,
                               RAIteratorVT        v_it) const
            {
                throw NotImplementedException();
            }

            // Note: this has to be const because changes to it could break 
            // the weights vector. 
            // further, using end() on the vector may iterate too far since
            // the vector is resized for the maximum number of vertices, 
            // not the amount currently stored.
            /// @todo mark: add iterator for GKC Sparse Vector
            const std::vector<IndexType>& getIndices()
            {
                return m_indices;
            }

            const std::vector<ScalarT>& getWeights()
            {
                if (m_weighted){
                    return m_weights;
                }
                throw NoValueException();
            }

            // output specific to the storage layout of this type of matrix
            void printInfo(std::ostream &os) const
            {
                os << "GKC Backend: ";
                os << "backend::GKCSparseVector<" << typeid(ScalarT).name() << ">";
                os << ", size  = " << m_num_vals;
                os << ", nvals = " << m_num_stored_vals << std::endl;

                os << "[";
                if (m_num_stored_vals > 0){
                    os << m_indices[0];
                    if (m_weighted){
                        os << ":" << m_weights[0];
                    } 
                }
                for (IndexType idx = 1; idx < m_num_stored_vals; ++idx)
                {
                    os << ", " << m_indices[idx];
                    if (m_weighted) {
                        os << ":" << m_weights[idx];
                    }
                }
                os << "]";
                os << std::endl;
            }

            friend std::ostream &operator<<(std::ostream             &os,
                                            GKCSparseVector<ScalarT> const &mat)
            {
                mat.printInfo(os);
                return os;
            }
            
        private:
            IndexType m_num_vals;
            IndexType m_num_stored_vals;
            bool m_weighted;

            // Two array compressed sparse vector
            std::vector<IndexType> m_indices;
            std::vector<ScalarType> m_weights;
        };

    } // namespace backend

} // namespace grb
