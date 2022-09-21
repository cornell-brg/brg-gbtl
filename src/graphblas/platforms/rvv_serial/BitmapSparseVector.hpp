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
#include <numeric>

namespace grb
{
    namespace backend
    {
        /**
         * @brief Class representing a sparse vector by using a bitmap + dense vector
         */
        template<typename ScalarT>
        class BitmapSparseVector
        {
        public:
            using ScalarType = ScalarT;

            /**
             * @brief Construct an empty sparse vector with given size
             *
             * @param[in] nsize  Size of vector.
             */
            BitmapSparseVector(IndexType nsize)
                : m_size(nsize),
                  m_vals(nsize),
#ifndef ARCH_RVV
                  m_bitmap(nsize, false)
#else
                  m_bitmap(nsize, 0)
#endif
            {
                if (nsize == 0)
                {
                    throw InvalidValueException();
                }
            }

            BitmapSparseVector(IndexType nsize, ScalarT const &value)
                : m_size(nsize),
                  m_vals(nsize, value),
#ifndef ARCH_RVV
                  m_bitmap(nsize, false)
#else
                  m_bitmap(nsize, 0)
#endif
            {
            }

            /**
             * @brief Construct from a dense vector.
             *
             * @param[in]  rhs  The dense vector to assign to this BitmapSparseVector.
             *                  Size is implied by the vector.
             * @return *this.
             */
            BitmapSparseVector(std::vector<ScalarT> const &rhs)
                : m_size(rhs.size()),
                  m_vals(rhs),
#ifndef ARCH_RVV
                  m_bitmap(rhs.size(), true)
#else
                  m_bitmap(rhs.size(), 1)
#endif
            {
                if (rhs.size() == 0)
                {
                    throw InvalidValueException();
                }
            }

            /**
             * @brief Construct a sparse vector from a dense array and zero val.
             *
             * @param[in]  rhs  The dense vector to assign to this BitmapSparseVector.
             *                  Size is implied by the vector.
             * @param[in]  zero An values in the rhs equal to this value will result
             *                  in an implied zero in the resulting sparse vector
             * @return *this.
             */
            BitmapSparseVector(std::vector<ScalarT> const &rhs,
                               ScalarT const              &zero)
                : m_size(rhs.size()),
                  m_vals(rhs.size()),
#ifndef ARCH_RVV
                  m_bitmap(rhs.size(), false)
#else
                  m_bitmap(rhs.size(), 0)
#endif
            {
                if (rhs.size() == 0)
                {
                    throw InvalidValueException();
                }

                for (IndexType idx = 0; idx < rhs.size(); ++idx)
                {
                    if (rhs[idx] != zero)
                    {
                        m_vals[idx] = rhs[idx];
#ifndef ARCH_RVV
                        m_bitmap[idx] = true;
#else
                        m_bitmap[idx] = 1;
#endif
                    }
                }
            }

            /**
             * @brief Construct from index and value arrays.
             * @deprecated Use vectorBuild method
             */
            BitmapSparseVector(
                IndexType                     nsize,
                std::vector<IndexType> const &indices,
                std::vector<ScalarT>   const &values)
                : m_size(nsize),
                  m_vals(nsize),
#ifndef ARCH_RVV
                  m_bitmap(nsize, false)
#else
                  m_bitmap(nsize, 0)
#endif
            {
                /// @todo check for same size indices and values
                for (IndexType idx = 0; idx < indices.size(); ++idx)
                {
                    IndexType i = indices[idx];
                    if (i >= m_size)
                    {
                        throw DimensionException();  // Should this be IndexOutOfBounds?
                    }

                    m_vals[i] = values[idx];
#ifndef ARCH_RVV
                    m_bitmap[i] = true;
#else
                    m_bitmap[i] = 1;
#endif
                }
            }

            /**
             * @brief Copy constructor for BitmapSparseVector.
             *
             * @param[in] rhs  The BitmapSparseVector to copy construct this
             *                 BitmapSparseVector from.
             */
            BitmapSparseVector(BitmapSparseVector<ScalarT> const &rhs)
                : m_size(rhs.m_size),
                  m_vals(rhs.m_vals),
                  m_bitmap(rhs.m_bitmap)
            {
            }

            ~BitmapSparseVector() {}

            /**
             * @brief Copy assignment.
             *
             * @param[in] rhs  The BitmapSparseVector to assign to this
             *
             * @return *this.
             */
            BitmapSparseVector<ScalarT>& operator=(
                BitmapSparseVector<ScalarT> const &rhs)
            {
                if (this != &rhs)
                {
                    if (m_size != rhs.m_size)
                    {
                        throw DimensionException();
                    }

                    m_vals = rhs.m_vals;
                    m_bitmap = rhs.m_bitmap;
                }
                return *this;
            }

            /**
             * @brief Assignment from a dense vector.
             *
             * @param[in]  rhs  The dense vector to assign to this BitmapSparseVector.
             *
             * @return *this.
             */
            BitmapSparseVector<ScalarT>& operator=(std::vector<ScalarT> const &rhs)
            {
                if (rhs.size() != m_size)
                {
                    throw DimensionException();
                }
                for (IndexType idx = 0; idx < rhs.size(); ++idx)
                {
                    m_vals[idx] = rhs[idx];
#ifndef ARCH_RVV
                    m_bitmap[idx] = true;
#else
                    m_bitmap[idx] = 1;
#endif
                }
                return *this;
            }

            // EQUALITY OPERATORS
            /**
             * @brief Equality testing for BitmapSparseVector.
             * @param rhs The right hand side of the equality operation.
             * @return If this BitmapSparseVector and rhs are identical.
             */
            bool operator==(BitmapSparseVector<ScalarT> const &rhs) const
            {
                if (m_size != rhs.m_size)
                {
                    return false;
                }

                for (IndexType i = 0; i < m_size; ++i)
                {
                    if (m_bitmap[i] != rhs.m_bitmap[i])
                    {
                        return false;
                    }
                    if (m_bitmap[i])
                    {
                        if (m_vals[i] != rhs.m_vals[i])
                        {
                            return false;
                        }
                    }
                }

                return true;
            }

            /**
             * @brief Inequality testing for BitmapSparseVector.
             * @param rhs The right hand side of the inequality operation.
             * @return If this BitmapSparseVector and rhs are not identical.
             */
            bool operator!=(BitmapSparseVector<ScalarT> const &rhs) const
            {
                return !(*this == rhs);
            }

        public:
            // METHODS

            void clear()
            {
                //m_vals.clear();
#ifndef ARCH_RVV
                m_bitmap.assign(m_size, false);
#else
                m_bitmap.assign(m_size, 0);
#endif
            }

            IndexType size() const
            {
                return m_size;
            }

            IndexType nvals() const
            {
                // count how many non-zero elements
                // @Tuan: TODO considering maintaining a variable m_nvals
                // instead of re-counting how many non-zeros everytime
#ifndef ARCH_RVV
                return std::count(m_bitmap.begin(), m_bitmap.end(), true);
#else
                size_t max_vlen = vsetvl_e32m1(m_size);
                size_t vlen     = 0;
                vuint32m1_t partial_sum_v = vmv_v_x(static_cast<uint32_t>(0), max_vlen);

                for (auto i = 0; i < m_size; i += vlen) {
                    vlen = vsetvl_e32m1(m_size - i);
                    partial_sum_v = vadd_vv(partial_sum_v,
                                            vle_v(m_bitmap.data() + i, vlen),
                                            vlen);
                }

                vuint32m1_t red_sum_v =
                        vredsum_vs_u32m1_u32m1(
                                red_sum_v,
                                partial_sum_v,
                                vmv_v_x(static_cast<uint32_t>(0), vlen),
                                max_vlen);
                return vmv_x_s_u32m1_u32(red_sum_v);
#endif
            }

            /**
             * @brief Resize the vector (smaller or larger)
             *
             * @param[in]  new_size  New number of elements (zero is invalid)
             *
             */
            void resize(IndexType new_size)
            {
                // Check in the frontend
                //if (nsize == 0)
                //   throw InvalidValueException();

                if (new_size < m_size)
                {
                    m_size = new_size;
                    m_bitmap.resize(new_size);
                    m_vals.resize(new_size);
                }
                else if (new_size > m_size)
                {
                    m_vals.resize(new_size);
#ifndef ARCH_RVV
                    m_bitmap.resize(new_size, false);
#else
                    m_bitmap.resize(new_size, 0);
#endif
                    m_size = new_size;
                }
            }

            /**
             *
             */
            template<typename RAIteratorIT,
                     typename RAIteratorVT,
                     typename BinaryOpT = grb::Second<ScalarType> >
            void build(RAIteratorIT  i_it,
                       RAIteratorVT  v_it,
                       IndexType     nvals,
                       BinaryOpT     dup = BinaryOpT())
            {
                std::vector<ScalarType> vals(m_size);
#ifndef ARCH_RVV
                std::vector<bool> bitmap(m_size);
#else
                std::vector<uint32_t> bitmap(m_size);
#endif

                /// @todo check for same size indices and values
                for (IndexType idx = 0; idx < nvals; ++idx)
                {
                    IndexType i = i_it[idx];
                    if (i >= m_size)
                    {
                        throw IndexOutOfBoundsException();
                    }

#ifndef ARCH_RVV
                    if (bitmap[i] == true)
#else
                    if (bitmap[i] == 1)
#endif
                    {
                        vals[i] = dup(vals[i], v_it[idx]);
                    }
                    else
                    {
                        vals[i] = v_it[idx];
#ifndef ARCH_RVV
                        bitmap[i] = true;
#else
                        bitmap[i] = 1;
#endif
                    }
                }

                m_vals.swap(vals);
                m_bitmap.swap(bitmap);
            }

            bool hasElement(IndexType index) const
            {
                if (index >= m_size)
                {
                    throw IndexOutOfBoundsException();
                }

#ifndef ARCH_RVV
                return m_bitmap[index];
#else
                return m_bitmap[index] == 1;
#endif
            }

            bool hasElementNoCheck(IndexType index) const
            {
#ifndef ARCH_RVV
                return m_bitmap[index];
#else
                return m_bitmap[index] == 1;
#endif
            }

#ifdef ARCH_RVV
            vbool32_t hasElementNoCheck(const RVVIndexType& index_v, size_t vlen) const
            {
                auto bitmap_v = vlxe_v(m_bitmap.data(), index_v, vlen);
                return vmseq_vx_u32m1_b32(bitmap_v, 1, vlen);
            }

            vbool32_t hasElementNoCheck(IndexType start_index, size_t vlen) const
            {
                auto bitmap_v = vle_v(m_bitmap.data() + start_index, vlen);
                return vmseq_vx_u32m1_b32(bitmap_v, 1, vlen);
            }
#endif

            /**
             * @brief Access the elements of this BitmapSparseVector given index.
             *
             * Function provided to access the elements of this BitmapSparseVector
             * given the index.
             *
             * @param[in] index  Position to access.
             *
             * @return The element of this BitmapSparseVector at the given row and
             *         column.
             */
            ScalarT extractElement(IndexType index) const
            {
                if (index >= m_size)
                {
                    throw IndexOutOfBoundsException();
                }

#ifndef ARCH_RVV
                if (m_bitmap[index] == false)
#else
                if (m_bitmap[index] == 0)
#endif
                {
                    throw NoValueException();
                }

                return m_vals[index];
            }

            ScalarT extractElementNoCheck(IndexType index) const
            {
                return m_vals[index];
            }

#ifdef ARCH_RVV
            auto extractElementNoCheck(const RVVIndexType& index_v, size_t vlen) const
            {
                return grb::vlxe_v(m_vals.data(), index_v, vlen);
            }

            auto extractElementNoCheck(IndexType start_index, size_t vlen) const
            {
                return grb::vle_v(m_vals.data() + start_index, vlen);
            }
#endif

            /// @todo Not certain about this implementation
            void setElement(IndexType      index,
                            ScalarT const &new_val)
            {
                if (index >= m_size)
                {
                    throw IndexOutOfBoundsException();
                }
                m_vals[index] = new_val;
#ifndef ARCH_RVV
                if (m_bitmap[index] == false)
                {
                    m_bitmap[index] = true;
                }
#else
                if (m_bitmap[index] == 0)
                {
                    m_bitmap[index] = 1;
                }
#endif
            }

            /// @todo Not certain about this implementation
            void setElementNoCheck(IndexType      index,
                                   ScalarT const &new_val)
            {
                m_vals[index] = new_val;
#ifndef ARCH_RVV
                if (m_bitmap[index] == false)
                {
                    m_bitmap[index] = true;
                }
#else
                if (m_bitmap[index] == 0)
                {
                    m_bitmap[index] = 1;
                }
#endif
            }

#ifdef ARCH_RVV
            template<typename VectorT>
            void setElementNoCheck(const RVVIndexType& index_vec,
                                   const VectorT&     new_val_vec,
                                   size_t             vlen)
            {
                if constexpr(std::is_same<ScalarT, uint32_t>::value) {
                    static_assert(std::is_same<VectorT, vuint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
                    static_assert(std::is_same<VectorT, vint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, float>::value) {
                    static_assert(std::is_same<VectorT, vfloat32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else {
                    static_assert(grb::always_false<ScalarT>, "vle_v: Unsupported type");
                }

                vsxe_v(m_vals.data(),
                       index_vec,
                       new_val_vec,
                       vlen);
                vsxe_v(m_bitmap.data(),
                       index_vec,
                       vmv_v_x(static_cast<grb::IndexType>(1), vlen),
                       vlen);
            }

            template<typename VectorT>
            void setElementNoCheck(IndexType          start_index,
                                   const VectorT&     new_val_vec,
                                   size_t             vlen)
            {
                if constexpr(std::is_same<ScalarT, uint32_t>::value) {
                    static_assert(std::is_same<VectorT, vuint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
                    static_assert(std::is_same<VectorT, vint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, float>::value) {
                    static_assert(std::is_same<VectorT, vfloat32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else {
                    static_assert(grb::always_false<ScalarT>, "vle_v: Unsupported type");
                }

                vse_v(m_vals.data() + start_index,
                      new_val_vec,
                      vlen);
                vse_v(m_bitmap.data() + start_index,
                      vmv_v_x(static_cast<grb::IndexType>(1), vlen),
                      vlen);
            }

            template<typename VectorT>
            void setElementNoCheck(const RVVIndexType& index_vec,
                                   const VectorT&     new_val_vec,
                                   const vbool32_t&   mask_vec,
                                   size_t             vlen)
            {
                if constexpr(std::is_same<ScalarT, uint32_t>::value) {
                    static_assert(std::is_same<VectorT, vuint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
                    static_assert(std::is_same<VectorT, vint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, float>::value) {
                    static_assert(std::is_same<VectorT, vfloat32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else {
                    static_assert(grb::always_false<ScalarT>, "vle_v: Unsupported type");
                }

                vsxe_v_m(m_vals.data(),
                         index_vec,
                         new_val_vec,
                         mask_vec,
                         vlen);
                vsxe_v_m(m_bitmap.data(),
                         index_vec,
                         vmv_v_x(static_cast<grb::IndexType>(1), vlen),
                         mask_vec,
                         vlen);
            }

            template<typename VectorT>
            void setElementNoCheck(IndexType          start_index,
                                   const VectorT&     new_val_vec,
                                   const vbool32_t&   mask_vec,
                                   size_t             vlen)
            {
                if constexpr(std::is_same<ScalarT, uint32_t>::value) {
                    static_assert(std::is_same<VectorT, vuint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
                    static_assert(std::is_same<VectorT, vint32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else if constexpr(std::is_same<ScalarT, float>::value) {
                    static_assert(std::is_same<VectorT, vfloat32m1_t>::value,
                                  "BitmapSparseVector::setElementNoCheck: Mismatched types");
                } else {
                    static_assert(grb::always_false<ScalarT>, "vle_v: Unsupported type");
                }

                vse_v_m(m_vals.data() + start_index,
                        new_val_vec,
                        mask_vec,
                        vlen);
                vse_v_m(m_bitmap.data() + start_index,
                        vmv_v_x(static_cast<grb::IndexType>(1), vlen),
                        mask_vec,
                        vlen);
            }
#endif

            void removeElement(IndexType index)
            {
                if (index >= m_size)
                {
                    throw IndexOutOfBoundsException();
                }
#ifndef ARCH_RVV
                if (m_bitmap[index] == true)
                {
                    m_bitmap[index] = false;
                }
#else
                if (m_bitmap[index] == 1)
                {
                    m_bitmap[index] = 0;
                }
#endif
            }

            void removeElementNoCheck(IndexType index)
            {
#ifndef ARCH_RVV
                if (m_bitmap[index] == true)
                {
                    m_bitmap[index] = false;
                }
#else
                if (m_bitmap[index] == 1)
                {
                    m_bitmap[index] = 0;
                }
#endif
            }

#ifdef ARCH_RVV
            void removeElementNoCheck(const RVVIndexType& index_vec,
                                      const vbool32_t&    mask_vec,
                                      size_t              vlen)
            {
                vsxe_v_m(m_bitmap.data(),
                         index_vec,
                         vmv_v_x(static_cast<grb::IndexType>(0), vlen),
                         mask_vec,
                         vlen);
            }

            void removeElementNoCheck(IndexType          start_index,
                                      const vbool32_t&   mask_vec,
                                      size_t             vlen)
            {
                vse_v_m(m_bitmap.data() + start_index,
                        vmv_v_x(static_cast<grb::IndexType>(0), vlen),
                        mask_vec,
                        vlen);
            }
#endif

            template<typename RAIteratorIT,
                     typename RAIteratorVT>
            void extractTuples(RAIteratorIT        i_it,
                               RAIteratorVT        v_it) const
            {
                for (IndexType idx = 0; idx < m_size; ++idx)
                {
                    if (m_bitmap[idx])
                    {
                        *i_it = idx;         ++i_it;
                        *v_it = m_vals[idx]; ++v_it;
                    }
                }
            }

            void extractTuples(IndexArrayType        &indices,
                               std::vector<ScalarT>  &values) const
            {
                extractTuples(indices.begin(), values.begin());
            }

            // output specific to the storage layout of this type of matrix
            void printInfo(std::ostream &os) const
            {
                os << "Optimized Sequential Backend: ";
                os << "backend::BitmapSparseVector<" << typeid(ScalarT).name() << ">";
                os << ", size  = " << m_size;
                os << ", nvals = " << nvals() << std::endl;

                os << "[";
                if (m_bitmap[0]) os << m_vals[0]; else os << "-";
                for (IndexType idx = 1; idx < m_size; ++idx)
                {
                    if (m_bitmap[idx]) os << ", " << m_vals[idx]; else os << ", -";
                }
                os << "]";
            }

            friend std::ostream &operator<<(std::ostream             &os,
                                            BitmapSparseVector<ScalarT> const &mat)
            {
                mat.printInfo(os);
                return os;
            }

        public:
#ifndef ARCH_RVV
            std::vector<bool> const &get_bitmap() const
#else
            std::vector<uint32_t> const &get_bitmap() const
#endif
            {
                return m_bitmap;
            }

            std::vector<ScalarT> const &get_vals() const { return m_vals; }

            std::vector<std::tuple<IndexType,ScalarT> > getContents() const
            {
                std::vector<std::tuple<IndexType,ScalarT> > contents;
                for (IndexType idx = 0; idx < m_size; ++idx)
                {
                    if (m_bitmap[idx])
                    {
                        contents.emplace_back(idx, m_vals[idx]);
                    }
                }
                return contents;
            }

            template <typename OtherScalarT>
            void setContents(
                std::vector<std::tuple<IndexType,OtherScalarT> > const &contents)
            {
                clear();
                for (auto&& [idx, val] : contents)
                {
#ifndef ARCH_RVV
                    m_bitmap[idx] = true;
#else
                    m_bitmap[idx] = 1;
#endif
                    m_vals[idx]   = static_cast<ScalarT>(val);
                }
            }

        private:
            IndexType             m_size;
            std::vector<ScalarT>  m_vals;
#ifndef ARCH_RVV
            std::vector<bool>     m_bitmap;
#else
            std::vector<uint32_t> m_bitmap;
#endif
        };
    } // backend
} // grb
