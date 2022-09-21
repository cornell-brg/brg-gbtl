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

#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <iostream>
#include <graphblas/algebra.hpp>

#include "sparse_helpers.hpp"

#ifdef ARCH_RVV
#include <riscv_vector.h>
#include "graphblas/rvv_helpers.hpp"
#endif

//****************************************************************************

namespace grb
{
    namespace backend
    {
        //********************************************************************
        /// Implementation of 4.3.2 vxm: u * A
        //********************************************************************
        template<typename WVectorT,
                 typename MaskT,
                 typename AccumT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AScalarT>
        inline void vxm(WVectorT                        &w,
                        MaskT                     const &mask,
                        AccumT                    const &accum,
                        SemiringT                        op,
                        UVectorT                  const &u,
                        CsrSparseMatrix<AScalarT> const &A,
                        OutputControlEnum                outp)
        {
            GRB_LOG_VERBOSE("w<M,z> := u +.* A");
            GRB_PROFILE("w<M,z> := u +.* A -- u.size=" + std::to_string(u.size()) + ", u.nvals=" + std::to_string(u.nvals()) + ", mask.nvals=" + std::to_string(mask.nvals()));

            // @Tuan: add some static check for input types
            static_assert(std::is_same<AScalarT, bool>::value == false,
                          "vxm() does not support AScalarT as bool");

            // =================================================================
            // Use axpy approach with the semi-ring.
            using TScalarType = typename SemiringT::result_type;
            BitmapSparseVector<TScalarType> t(w.size());

            if ((A.nvals() > 0) && (u.nvals() > 0))
            {
                for (IndexType row_idx = 0; row_idx < u.size(); ++row_idx)
                {
                    IndexType        row_nvals   = A.getNvals(row_idx);
                    const IndexType* col_idx_arr = A.getColIdxArr(row_idx);
                    const AScalarT*  val_arr     = A.getValArr(row_idx);

                    if (u.hasElementNoCheck(row_idx) && row_nvals != 0) {
                        // @Tuan: doing an axpy(u[row_idx], A[row_idx]) -> t
                        auto a = u.extractElementNoCheck(row_idx);
#ifndef ARCH_RVV
                        for (IndexType i = 0; i < row_nvals; ++i)
                        {
                            auto j   = col_idx_arr[i];
                            auto b_j = val_arr[i];
                            auto t_j = op.mult(a, b_j);

                            if (t.hasElementNoCheck(j))
                                t.setElementNoCheck(j, op.add(t.extractElementNoCheck(j), t_j));
                            else
                                t.setElementNoCheck(j, t_j);
                        }
#else
                        size_t vlen = 0;
                        for (IndexType i = 0; i < row_nvals; i += vlen) {
                            vlen = vsetvl_e32m1(row_nvals - i);

                            auto j_vec   = vle_v(col_idx_arr + i, vlen);
                            auto b_j_vec = vle_v(val_arr + i, vlen);
                            auto t_j_vec = op.mult(grb::vmv_v_x(a, vlen), b_j_vec, vlen);

                            // accumulate in t_j_vec
                            t_j_vec = op.add(vlxe_v(t.get_vals().data(), j_vec, vlen), t_j_vec, vlen);

                            // update t
                            t.setElementNoCheck(j_vec, t_j_vec, vlen);
                        }
#endif
                    }
                }
            }

            // =================================================================
            // Accumulate into final output, w, considering mask and replace/merge
            //using ZScalarType = typename std::conditional_t<
            //    std::is_same_v<AccumT, NoAccumulate>,
            //    TScalarType,
            //    decltype(accum(std::declval<typename WVectorT::ScalarType>(),
            //                   std::declval<TScalarType>()))>;
            opt_accum_with_opt_mask_1D(w, mask, accum, t, outp);
        }

        //**********************************************************************
        //**********************************************************************
        //**********************************************************************

        //**********************************************************************
        /// Implementation of 4.3.2 vxm: u * A'
        //**********************************************************************
        template<typename WVectorT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AMatrixT>
        inline void vxm_dot_nomask_noaccum_noalias(
            WVectorT                      &w,
            NoMask                  const &mask,
            NoAccumulate            const &accum,
            SemiringT                      op,
            UVectorT                const &u,
            TransposeView<AMatrixT> const &AT,
            OutputControlEnum              outp)
        {
            auto const &A(AT.m_mat);
            w.clear();  // ERROR if u and w are same vector

            // =================================================================
            // Do the basic dot-product work with the semi-ring.

            throw grb::NotImplementedException(
                  std::string("INTERNAL ERROR at ") + __FILE__ + ":" + std::to_string(__LINE__));
        }

        //**********************************************************************
        template<typename WVectorT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AMatrixT>
        inline void vxm(WVectorT                      &w,
                        NoMask                  const &mask,
                        NoAccumulate            const &accum,
                        SemiringT                      op,
                        UVectorT                const &u,
                        TransposeView<AMatrixT> const &AT,
                        OutputControlEnum              outp)
        {
            //mxv_dot_nomask_noaccum(w, op, A, u);
            GRB_LOG_VERBOSE("w := u +.* A'");
            auto const &A(AT.m_mat);

            throw grb::NotImplementedException(
                  std::string("INTERNAL ERROR at ") + __FILE__ + ":" + std::to_string(__LINE__));
        }

        //**********************************************************************
        /// Implementation for 4.3.2 mxv: w + u * A'
        //**********************************************************************
        template<typename WVectorT,
                 typename AccumT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AMatrixT>
        inline void vxm(WVectorT           &w,
                        NoMask       const &mask,
                        AccumT       const &accum,
                        SemiringT           op,
                        UVectorT     const &u,
                        TransposeView<AMatrixT> const &AT,
                        OutputControlEnum  outp)
        {
            //mxv_dot_nomask_accum(w, accum, op, A, u);
            GRB_LOG_VERBOSE("w := w + (u +.* A')");
            auto const &A(AT.m_mat);

            // =================================================================
            // Do the basic dot-product work with the semi-ring.
            using TScalarType = typename SemiringT::result_type;

            throw grb::NotImplementedException(
                  std::string("INTERNAL ERROR at ") + __FILE__ + ":" + std::to_string(__LINE__));
        }

        //**********************************************************************
        /// Implementation for 4.3.2 mxv: <m,r>(u * A')
        //**********************************************************************
        template<typename WVectorT,
                 typename MaskT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AMatrixT>
        inline void vxm_dot_mask_noaccum_noalias(
            WVectorT                      &w,
            MaskT                   const &mask,
            NoAccumulate            const &accum,
            SemiringT                      op,
            UVectorT                const &u,
            TransposeView<AMatrixT> const &AT,
            OutputControlEnum              outp)
        {
            auto const &A(AT.m_mat);

            // =================================================================
            // Do the basic dot-product work with the semi-ring.
            using TScalarType = typename SemiringT::result_type;
            std::vector<std::tuple<IndexType, TScalarType> > t;

            throw grb::NotImplementedException(
                  std::string("INTERNAL ERROR at ") + __FILE__ + ":" + std::to_string(__LINE__));
        }

        //**********************************************************************
        template<typename WVectorT,
                 typename MaskT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AMatrixT>
        inline void vxm(WVectorT                      &w,
                        MaskT                   const &mask,
                        NoAccumulate            const &accum,
                        SemiringT                      op,
                        UVectorT                const &u,
                        TransposeView<AMatrixT> const &AT,
                        OutputControlEnum              outp)
        {
            //mxv_dot_mask_noaccum(w, mask, op, A, u, outp);
            GRB_LOG_VERBOSE("w<M,r> := (u +.* A')");

            throw grb::NotImplementedException(
                  std::string("INTERNAL ERROR at ") + __FILE__ + ":" + std::to_string(__LINE__));
        }

        //**********************************************************************
        /// Implementation for 4.3.2 mxv: <m,r>(w + u * A')
        //**********************************************************************
        template<typename WVectorT,
                 typename MaskT,
                 typename AccumT,
                 typename SemiringT,
                 typename UVectorT,
                 typename AMatrixT>
        inline void vxm(WVectorT          &w,
                        MaskT       const &mask,
                        AccumT      const &accum,
                        SemiringT          op,
                        UVectorT    const &u,
                        TransposeView<AMatrixT> const &AT,
                        OutputControlEnum  outp)
        {
            //mxv_dot_mask_accum(w, mask, accum, op, A, u, outp);
            GRB_LOG_VERBOSE("w<M,r> := w + (u +.* A')");
            auto const &A(AT.m_mat);

            // =================================================================
            // Do the basic dot-product work with the semi-ring.
            using TScalarType = typename SemiringT::result_type;

            // =================================================================
            // Accumulate into Z
            using ZScalarType =
                decltype(accum(std::declval<typename WVectorT::ScalarType>(),
                               std::declval<TScalarType>()));
            std::vector<std::tuple<IndexType, ZScalarType> > z;

            throw grb::NotImplementedException(
                  std::string("INTERNAL ERROR at ") + __FILE__ + ":" + std::to_string(__LINE__));
        }
    } // backend
} // grb
