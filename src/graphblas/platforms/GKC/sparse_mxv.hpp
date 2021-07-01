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


//****************************************************************************

namespace grb
{
    namespace backend
    {
        //**********************************************************************
        /// Implementation for 4.3.3 mxv: A * u
        //**********************************************************************
        template<typename WVectorT,
                 typename MaskT,
                 typename AccumT,
                 typename SemiringT,
                 typename AMatrixT,
                 typename UVectorT>
        inline void mxv(WVectorT          &w,
                        MaskT       const &mask,
                        AccumT      const &accum,
                        SemiringT          op,
                        AMatrixT    const &A,
                        UVectorT    const &u,
                        OutputControlEnum  outp)
        {
            GRB_LOG_VERBOSE("w<M,z> := A +.* u");

            // =================================================================
            // Do the basic dot-product work with the semi-ring.
            using TScalarType = typename SemiringT::result_type;
            std::vector<std::tuple<IndexType, TScalarType> > t;

            if ((A.nvals() > 0) && (u.nvals() > 0))
            {
                auto u_contents(u.getContents());
                for (IndexType row_idx = 0; row_idx < w.size(); ++row_idx)
                {
                    if (!A[row_idx].empty())
                    {
                        TScalarType t_val;
                        /// @note In mxv_timing_test, if I reverse u_contents and
                        /// A[row_idx], the performance improves by a factor of 2.
                        /// But I cannot reorder in case op is not commutative.
                        ///
                        /// I have added dot_rev() helper that reverses the two
                        /// vectors but keeps the order correct for op.
                        ///
                        /// I suspect this is strictly data dependent performance
                        if (dot_rev(t_val, A[row_idx], u_contents, op))
                        {
                            t.emplace_back(row_idx, t_val);
                        }
                    }
                }
            }

            // =================================================================
            // Accumulate into Z
            using ZScalarType = typename std::conditional_t<
                std::is_same_v<AccumT, NoAccumulate>,
                TScalarType,
                decltype(accum(std::declval<typename WVectorT::ScalarType>(),
                               std::declval<TScalarType>()))>;

            std::vector<std::tuple<IndexType, ZScalarType> > z;
            ewise_or_opt_accum_1D(z, w, t, accum);

            // =================================================================
            // Copy Z into the final output, w, considering mask and replace/merge
            write_with_opt_mask_1D(w, z, mask, outp);
        }

        //**********************************************************************
        //**********************************************************************
        //**********************************************************************

        //**********************************************************************
        /// Implementation of 4.3.3 mxv: A' * u
        //**********************************************************************
        template<typename WVectorT,
                 typename MaskT,
                 typename AccumT,
                 typename SemiringT,
                 typename AMatrixT,
                 typename UVectorT>
        inline void mxv(WVectorT                      &w,
                        MaskT                   const &mask,
                        AccumT                  const &accum,
                        SemiringT                      op,
                        TransposeView<AMatrixT> const &AT,
                        UVectorT                const &u,
                        OutputControlEnum              outp)
        {
            GRB_LOG_VERBOSE("w<M,z> := A' +.* u");
            auto const &A(AT.m_mat);

            // =================================================================
            // Use axpy approach with the semi-ring.
            using TScalarType = typename SemiringT::result_type;
            std::vector<std::tuple<IndexType, TScalarType> > t;

            if ((A.nvals() > 0) && (u.nvals() > 0))
            {
                for (IndexType row_idx = 0; row_idx < u.size(); ++row_idx)
                {
                    if (u.hasElement(row_idx) && !A[row_idx].empty())
                    {
                        axpy(t, op, u.extractElement(row_idx), A[row_idx]);
                    }
                }
            }

            // =================================================================
            // Accumulate into Z
            using ZScalarType = typename std::conditional_t<
                std::is_same_v<AccumT, NoAccumulate>,
                TScalarType,
                decltype(accum(std::declval<typename WVectorT::ScalarType>(),
                               std::declval<TScalarType>()))>;

            std::vector<std::tuple<IndexType, ZScalarType> > z;
            ewise_or_opt_accum_1D(z, w, t, accum);

            // =================================================================
            // Copy Z into the final output, w, considering mask and replace/merge
            write_with_opt_mask_1D(w, z, mask, outp);
        }
        
        //**********************************************************************
        /// Implementation for mxv with GKC Matrix and GKC Sparse Vector: A * u
        //**********************************************************************
        // w, mask, and u are vectors. A is a matrix.
        template <typename AccumT,
                  typename MaskT,
                  typename SemiringT,
                  typename ScalarT>
        inline void mxv(GKCSparseVector<ScalarT> &w,
                        MaskT const &mask,
                        AccumT const &accum,
                        SemiringT op,
                        GKCMatrix<ScalarT> const &A,
                        GKCSparseVector<ScalarT> const &u,
                        OutputControlEnum outp)
        {
            GRB_LOG_VERBOSE("w<M,z> := A +.* u");

            // =================================================================
            // Do the basic dot-product work with the semi-ring.
            using TScalarType = typename SemiringT::result_type;

            GKCSparseVector<TScalarType> t(u.nvals());
            u.sortSelf();
            //std::vector<std::tuple<IndexType, TScalarType> > t;

            if ((A.nvals() > 0) && (u.nvals() > 0))
            {
                //auto u_contents(u.getContents());
                for (IndexType row_idx = 0; row_idx < w.size(); ++row_idx) 
                {
                    auto AIst = A.idxBegin(row_idx);
                    auto AInd = A.idxEnd(row_idx);
                    auto AWst = A.wgtBegin(row_idx);
                    auto AWnd = A.wgtEnd(row_idx);
                    auto UIst = u.getIndices().begin();
                    auto UInd = u.getIndices().end(); 
                    auto UWst = u.getWeights().begin();
                    auto UWnd = u.getWeights().end(); 
                    // Do dot product here, into t
                    bool value_set(false);
                    TScalarType sum;
                    //if (AIst < AInd && UIst < UInd) std::cout << "row idx: " << row_idx << std::endl;
                    while (AIst < AInd && UIst < UInd)
                    {                            
                        if (*AIst == *UIst)
                        {
                            //std::cout << "AIdx: " << *AIst << " UIdx: " << *UIst << std::endl;
                            //std::cout << "AWgt: " << *AWst<< " UWgt: " << *UWst << std::endl;
                            if (value_set)
                            {
                                sum = op.add(sum, op.mult(*AWst, *UWst));
                            }
                            else 
                            {                                
                                sum = op.mult(*AWst, *UWst);
                                value_set = true;
                            }
                            AIst++; AWst++;
                            UIst++; UWst++;
                        }
                        else if (*AIst < *UIst)
                        {                                
                            AIst++; AWst++;
                        }
                        else 
                        {
                            UIst++; UWst++;
                        }
                        //if (value_set) std::cout << "Sum: " << sum << " ";
                    }
                    /// @todo: outputs control enum, masking, etc...
                    if (value_set) w.setElement(row_idx, sum);
                    //std::cout << std::endl;
                }
            }
            w.sortSelf();

            // =================================================================
            // Accumulate into Z
            //!!using ZScalarType = typename std::conditional_t<
            //!!    std::is_same_v<AccumT, NoAccumulate>,
            //!!    TScalarType,
            //!!    decltype(accum(std::declval<typename WVectorT::ScalarType>(),
            //!!                   std::declval<TScalarType>()))>;

            //!!std::vector<std::tuple<IndexType, ZScalarType> > z;
            //!!ewise_or_opt_accum_1D(z, w, t, accum);

            // =================================================================
            // Copy Z into the final output, w, considering mask and replace/merge
            //!!write_with_opt_mask_1D(w, z, mask, outp);
        }

        
        //**********************************************************************
        //**********************************************************************
        //**********************************************************************

        //**********************************************************************
        /// Implementation of mxv for GKC Matrix and Sparse Vector: A' * u
        //**********************************************************************
        template <typename MaskT,
                  typename AccumT,
                  typename SemiringT,
                  typename ScalarT>
        inline void mxv(GKCSparseVector<ScalarT> &w,
                        MaskT const &mask,
                        AccumT const &accum,
                        SemiringT op,
                        TransposeView<GKCMatrix<ScalarT>> const &AT,
                        GKCSparseVector<ScalarT> const &u,
                        OutputControlEnum outp)
        {
            GRB_LOG_VERBOSE("w<M,z> := A' +.* u");
            auto const &A(AT.m_mat);

            // =================================================================
            // Use axpy approach with the semi-ring.
            using TScalarType = typename SemiringT::result_type;
            GKCSparseVector<TScalarType> t();
            
            if (outp == REPLACE) w.clear();

            if ((A.nvals() > 0) && (u.nvals() > 0))
            {
                auto UIst = u.getIndices().begin();
                auto UInd = u.getIndices().end();
                auto UWst = u.getWeights().begin();
                for ( ; UIst < UInd; UIst++, UWst++)
                {
                    bool value_set(false);
                    TScalarType res;
                    auto AIst = A.idxBegin(*UIst); 
                    auto AWst = A.wgtBegin(*UIst);
                    for ( ; AIst < A.idxEnd(*UIst); AIst++, AWst++)
                    {
                        res = op.mult(*UWst, *AWst);
                        w.mergeSetElement(*AIst, res, op);
                        //std::cout << "Wrote val: " << res << " at idx: " << *AIst << std::endl;
                        //std::cout << "Current val is " << w.extractElement(*AIst) << std::endl;
                    }
                }
            }

            // =================================================================
            // Accumulate into Z
            //using ZScalarType = typename std::conditional_t<
            //    std::is_same_v<AccumT, NoAccumulate>,
            //    TScalarType,
            //    decltype(accum(std::declval<typename WVectorT::ScalarType>(),
            //                   std::declval<TScalarType>()))>;

            //std::vector<std::tuple<IndexType, ZScalarType> > z;
            //ewise_or_opt_accum_1D(z, w, t, accum);

            // =================================================================
            // Copy Z into the final output, w, considering mask and replace/merge
            //write_with_opt_mask_1D(w, z, mask, outp);
        }
    } // backend
} // grb
