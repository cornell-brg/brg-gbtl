#ifdef ARCH_RVV

#pragma once

#include <riscv_vector.h>

#include <graphblas/types.hpp>

namespace grb
{

//---------------------------------------------------------------------------
// vector arithmetics
//---------------------------------------------------------------------------

template<typename VectorT>
inline auto vadd_vv(const VectorT& lhs, const VectorT& rhs, size_t vlen)
{
    if constexpr(std::is_same<VectorT, vuint32m1_t>::value) {
        return vadd_vv_u32m1(lhs, rhs, vlen);
    } else if constexpr(std::is_same<VectorT, vint32m1_t>::value) {
        return vadd_vv_i32m1(lhs, rhs, vlen);
    } else if constexpr(std::is_same<VectorT, vfloat32m1_t>::value) {
        return vfadd_vv_f32m1(lhs, rhs, vlen);
    } else {
        static_assert(grb::always_false<VectorT>, "vadd_vv: Unsupported type");
    }
}

template<typename VectorT>
inline auto vmul_vv(const VectorT& lhs, const VectorT& rhs, size_t vlen)
{
    if constexpr(std::is_same<VectorT, vuint32m1_t>::value) {
        return vmul_vv_u32m1(lhs, rhs, vlen);
    } else if constexpr(std::is_same<VectorT, vint32m1_t>::value) {
        return vmul_vv_i32m1(lhs, rhs, vlen);
    } else if constexpr(std::is_same<VectorT, vfloat32m1_t>::value) {
        return vfmul_vv_f32m1(lhs, rhs, vlen);
    } else {
        static_assert(grb::always_false<VectorT>, "vmul_vv: Unsupported type");
    }
}

template<typename VectorT>
inline auto vor_vv(const VectorT& lhs, const VectorT& rhs, size_t vlen)
{
    if constexpr(std::is_same<VectorT, vuint32m1_t>::value) {
        return vor_vv_u32m1(lhs, rhs, vlen);
    } else if constexpr(std::is_same<VectorT, vint32m1_t>::value) {
        return vor_vv_i32m1(lhs, rhs, vlen);
    } else {
        static_assert(grb::always_false<VectorT>, "vor_vv: Unsupported type");
    }
}

template<typename VectorT>
inline auto vand_vv(const VectorT& lhs, const VectorT& rhs, size_t vlen)
{
    if constexpr(std::is_same<VectorT, vuint32m1_t>::value) {
        return vand_vv_u32m1(lhs, rhs, vlen);
    } else if constexpr(std::is_same<VectorT, vint32m1_t>::value) {
        return vand_vv_i32m1(lhs, rhs, vlen);
    } else {
        static_assert(grb::always_false<VectorT>, "vand_vv: Unsupported type");
    }
}

//---------------------------------------------------------------------------
// move
//---------------------------------------------------------------------------

template<typename ScalarT>
inline auto vmv_v_x(const ScalarT& scalar, size_t vlen)
{
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        return vmv_v_x_u32m1(scalar, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        return vmv_v_x_i32m1(scalar, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        return vfmv_v_f_f32m1(scalar, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vmv_v_x: Unsupported type");
    }
}

//---------------------------------------------------------------------------
// vector memory
//---------------------------------------------------------------------------

template<typename ScalarT>
inline auto vle_v(const ScalarT* addr, size_t vlen)
{
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        return vle32_v_u32m1(addr, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        return vle32_v_i32m1(addr, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        return vle32_v_f32m1(addr, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vle_v: Unsupported type");
    }
}

template<typename ScalarT>
inline auto vlxe_v(const ScalarT* addr, const vuint32m1_t& index_vec, size_t vlen)
{
    auto index_vec_byte = vmul_vx_u32m1(index_vec, sizeof(uint32_t), vlen);
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        return vloxei32_v_u32m1(addr, index_vec_byte, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        return vloxei32_v_i32m1(addr, index_vec_byte, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        return vloxei32_v_f32m1(addr, index_vec_byte, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vlxe_v: Unsupported type");
    }
}

template<typename ScalarT, typename VectorT>
inline void vse_v(ScalarT* addr, const VectorT& val_vec, size_t vlen)
{
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        static_assert(std::is_same<VectorT, vuint32m1_t>::value, "vse_v: Mismatched type");
        vse32_v_u32m1(addr, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        static_assert(std::is_same<VectorT, vint32m1_t>::value, "vse_v: Mismatched type");
        vse32_v_i32m1(addr, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        static_assert(std::is_same<VectorT, vfloat32m1_t>::value, "vse_v: Mismatched type");
        vse32_v_f32m1(addr, val_vec, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vse_v: Unsupported type");
    }
}

template<typename ScalarT, typename VectorT>
inline void vse_v_m(ScalarT* addr, const VectorT& val_vec, const vbool32_t& mask_vec, size_t vlen)
{
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        static_assert(std::is_same<VectorT, vuint32m1_t>::value, "vse_v_m: Mismatched type");
        vse32_v_u32m1_m(mask_vec, addr, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        static_assert(std::is_same<VectorT, vint32m1_t>::value, "vse_v_m: Mismatched type");
        vse32_v_i32m1_m(mask_vec, addr, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        static_assert(std::is_same<VectorT, vfloat32m1_t>::value, "vse_v_m: Mismatched type");
        vse32_v_f32m1_m(mask_vec, addr, val_vec, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vse_v_m: Unsupported type");
    }
}

template<typename ScalarT, typename VectorT>
inline void vsxe_v(ScalarT* addr, const vuint32m1_t& index_vec,
                   const VectorT& val_vec, size_t vlen)
{
    auto index_vec_byte = vmul_vx_u32m1(index_vec, sizeof(uint32_t), vlen);
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        static_assert(std::is_same<VectorT, vuint32m1_t>::value, "vsxe_v: Mismatched type");
        vsoxei32_v_u32m1(addr, index_vec_byte, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        static_assert(std::is_same<VectorT, vint32m1_t>::value, "vsxe_v: Mismatched type");
        vsoxei32_v_i32m1(addr,index_vec_byte, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        static_assert(std::is_same<VectorT, vfloat32m1_t>::value, "vsxe_v: Mismatched type");
        vsoxei32_v_f32m1(addr,index_vec_byte, val_vec, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vsxe_v: Unsupported type");
    }
}

template<typename ScalarT, typename VectorT>
inline void vsxe_v_m(const ScalarT* addr, const vuint32m1_t& index_vec,
                     const VectorT& val_vec, const vbool32_t& mask_vec, size_t vlen)
{
    auto index_vec_byte = vmul_vx_u32m1(index_vec, sizeof(uint32_t), vlen);
    if constexpr(std::is_same<ScalarT, uint32_t>::value) {
        static_assert(std::is_same<VectorT, vuint32m1_t>::value, "vsxe_v_m: Mismatched type");
        vsoxei32_v_u32m1_m(mask_vec, addr, index_vec_byte, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, int32_t>::value) {
        static_assert(std::is_same<VectorT, vint32m1_t>::value, "vsxe_v_m: Mismatched type");
        vsoxei32_v_i32m1_m(mask_vec, addr,index_vec_byte, val_vec, vlen);
    } else if constexpr(std::is_same<ScalarT, float>::value) {
        static_assert(std::is_same<VectorT, vfloat32m1_t>::value, "vsxe_v_m: Mismatched type");
        vsoxei32_v_f32m1_m(mask_vec, addr,index_vec_byte, val_vec, vlen);
    } else {
        static_assert(grb::always_false<ScalarT>, "vsxe_v_m: Unsupported type");
    }
}

} // end namespace grb

#endif // ARCH_RVV
