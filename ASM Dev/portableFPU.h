/** 
 * bits to set the floating point control word register
 *
 * Sections 4.9, 8.1.4, 10.2.2 and 11.5 in 
 * IA-32 Intel Architecture Software Developer's Manual
 *   Volume 1: Basic Architecture
 *
 * http://www.intel.com/design/pentium4/manuals/245471.htm
 *http://www.stereopsis.com/sree/fpu2006.html
 * http://www.geisswerks.com/ryan/FAQS/fpu.html
 *
 * windows has _controlfp() but it takes different parameters
 *
 * 0 : IM invalid operation mask
 * 1 : DM denormalized operand mask
 * 2 : ZM divide by zero mask
 * 3 : OM overflow mask
 * 4 : UM underflow mask
 * 5 : PM precision, inexact mask
 * 6,7 : reserved
 * 8,9 : PC precision control
 * 10,11 : RC rounding control
 *
 * precision control:
 * 00 : single precision
 * 01 : reserved
 * 10 : double precision
 * 11 : extended precision
 *
 * rounding control:
 * 00 = Round to nearest whole number. (default)
 * 01 = Round down, toward -infinity.
 * 10 = Round up, toward +infinity.
 * 11 = Round toward zero (truncate).
 */
#ifndef _PORTABLEFPU_H_INCLUDED
#define _PORTABLEFPU_H_INCLUDED
#endif

#ifndef _PSTDINT_H_INCLUDED
#include "pstdint.h"
#endif

#define __FPU_CW_EXCEPTION_MASK__   (0x003f)
#define __FPU_CW_INVALID__          (0x0001)
#define __FPU_CW_DENORMAL__         (0x0002)
#define __FPU_CW_ZERODIVIDE__       (0x0004)
#define __FPU_CW_OVERFLOW__         (0x0008)
#define __FPU_CW_UNDERFLOW__        (0x0010)
#define __FPU_CW_INEXACT__          (0x0020)

#define __FPU_CW_PREC_MASK__        (0x0300)
#define __FPU_CW_PREC_SINGLE__      (0x0000)
#define __FPU_CW_PREC_DOUBLE__      (0x0200)
#define __FPU_CW_PREC_EXTENDED__    (0x0300)

#define __FPU_CW_ROUND_MASK__       (0x0c00)
#define __FPU_CW_ROUND_NEAR__       (0x0000)
#define __FPU_CW_ROUND_DOWN__       (0x0400)
#define __FPU_CW_ROUND_UP__         (0x0800)
#define __FPU_CW_ROUND_CHOP__       (0x0c00)

#define __FPU_CW_MASK_ALL__         (0x1f3f)


#define __SSE_CW_FLUSHZERO__        (0x8000)
    
#define __SSE_CW_ROUND_MASK__       (0x6000)
#define __SSE_CW_ROUND_NEAR__       (0x0000)
#define __SSE_CW_ROUND_DOWN__       (0x2000)
#define __SSE_CW_ROUND_UP__         (0x4000)
#define __SSE_CW_ROUND_CHOP__       (0x6000)

#define __SSE_CW_EXCEPTION_MASK__   (0x1f80)
#define __SSE_CW_PRECISION__        (0x1000)
#define __SSE_CW_UNDERFLOW__        (0x0800)
#define __SSE_CW_OVERFLOW__         (0x0400)
#define __SSE_CW_DIVIDEZERO__       (0x0200)
#define __SSE_CW_DENORMAL__         (0x0100)
#define __SSE_CW_INVALID__          (0x0080)
// not on all SSE machines
// #define __SSE_CW_DENORMALZERO__     (0x0040)

#define __SSE_CW_MASK_ALL__         (0xffc0)

#define __MOD_FPU_CW_DEFAULT__ (__FPU_CW_EXCEPTION_MASK__ + __FPU_CW_PREC_DOUBLE__ + __FPU_CW_ROUND_CHOP__)
#define __MOD_SSE_CW_DEFAULT__ (__SSE_CW_EXCEPTION_MASK__ + __SSE_CW_ROUND_CHOP__ + __SSE_CW_FLUSHZERO__)


#ifdef __WIN32__
#define ___INLINE __inline
#elif defined(__GNUC__)
#define ___INLINE __inline__
#endif

#ifdef USE_SSE
___INLINE uint32_t getSSEStateX86(void);
___INLINE void setSSEModDefault(uint32_t control);
___INLINE void modifySSEStateX86(uint32_t control, uint32_t mask);
#endif // USE_SSE

___INLINE void setRoundingMode(uint32_t round);
___INLINE uint32_t getFPUStateX86(void);
___INLINE void setFPUStateX86(uint32_t control);
___INLINE void setFPUModDefault(void);
___INLINE void assertFPUModDefault(void);
___INLINE void modifyFPUStateX86(const uint32_t control, const uint32_t mask);
___INLINE int32_t FastFtol(const float a);
