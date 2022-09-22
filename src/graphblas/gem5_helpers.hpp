#pragma once

namespace gem5
{

inline void toggle_stats(bool on)
{
#if defined(ARCH_RVV) && defined(GEM5)
  __asm__ volatile ("csrw 0x7C1, %0;"
                    :
                    : "r" (on)
                    :);
  // toggle activity stats as well
  __asm__ volatile ("csrw 0x800, %0;"
                    :
                    : "r" (on)
                    :);
#endif
}

inline void switch_cpus(bool to_detailed)
{
#if defined(ARCH_RVV) && defined(GEM5)
  __asm__ volatile ("csrw 0x7C2, %0;"
                    :
                    : "r" (to_detailed)
                    :);
#endif
}

inline void vmfence()
{
#if defined(ARCH_RVV) && defined(GEM5)
  __asm__ volatile ("csrrwi zero, 0x83F, 0x0 \n"
                    :::);
#endif
}

inline void vfence()
{
#if defined(ARCH_RVV) && defined(GEM5)
  __asm__ volatile ("csrrwi zero, 0x83E, 0x0 \n"
                    :::);
#endif
}

inline bool has_vector_support()
{
#if defined(ARCH_RVV) && defined(GEM5)
  uint64_t csr_value = 0;
  __asm__ volatile ("csrr %[csr_value], 0xCC0"
                    : [csr_value] "=r" (csr_value)
                    :
                    :);
  return csr_value == 1;
#else
  return false;
#endif
}

inline void vstart()
{
#if defined(ARCH_RVV) && defined(GEM5)
  if (gem5::has_vector_support()) {
    uint64_t vs_mask = 0x600;
    asm volatile ("csrrs x0, ustatus, %[vs_mask]; \n"
                    :
                    : [vs_mask] "r" (vs_mask)
                    :);

  }
#endif
}

inline void vend()
{
#if defined(ARCH_RVV) && defined(GEM5)
  if (gem5::has_vector_support()) {
    uint64_t vs_mask = 0x600;
    gem5::vmfence();
    gem5::vfence();
    asm volatile ("csrrc x0, ustatus, %[vs_mask]; \n"
                    :
                    : [vs_mask] "r" (vs_mask)
                    :);
  }
#endif
}

}
