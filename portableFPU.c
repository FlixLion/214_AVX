//assume for now that we are running on an x86
// #ifdef __i386__

#ifndef __WIN32__
    #ifdef WIN32
    #define __WIN32__
    #endif
    #ifdef _WIN32
    #define __WIN32__
    #endif
    #ifdef __WIN32
    #define __WIN32__
    #endif
    #ifdef _WIN64
    #define __WIN32__
    #endif
    #ifdef WIN64
    #define __WIN32__
    #endif
    #ifdef _WINDOWS
    #define __WIN32__
    #endif
#endif

#ifndef _PSTDINT_H_INCLUDED
#include "pstdint.h"
#endif

#ifndef _PORTABLEFPU_H_INCLUDED
#include "portableFPU.h"
#endif

#ifdef USE_SSE
___INLINE unsigned int getSSEStateX86(void)
{
#ifdef __WIN32__
    return _FPU_SETCW(control);
#else
    return _mm_getcsr();
#endif
}

___INLINE void setSSEStateX86(unsigned int control)
{
#ifdef __WIN32__
    _FPU_GETCW(control);
#else
    _mm_setcsr(control);
#endif
}

___INLINE void modifySSEStateX86(uint32_t control, uint32_t mask);
{
    unsigned int oldControl = getFPUStateX86();
    unsigned int newControl = ((oldControl & (~mask)) | (control & mask));
    setFPUStateX86(newControl);
}

___INLINE void setSSEModDefault(void)
{
    modifySSEStateX86(__MOD_SSE_CW_DEFAULT__, __SSE_CW_MASK_ALL__)
}    
#endif // USE_SSE

___INLINE void setRoundingMode(uint32_t round)
{
    //ASSERT(round < 4);
    uint32_t mask = 0x3;
    uint32_t fpuControl = getFPUStateX86();
    fpuControl &= ~(mask << 10);
    fpuControl |= round << 10;
    setFPUStateX86(fpuControl);

#ifdef USE_SSE
    uint32_t sseControl = getSSEStateX86();
    sseControl &= ~(mask << 13);
    sseControl |= round << 13;
    setSSEStateX86(sseControl);
#endif // USE_SSE
}


uint32_t getFPUStateX86(void)
{
    uint32_t control = 0;
#if defined(_MSVC)
    __asm fnstcw control;
#elif defined(__GNUG__)
    __asm__ __volatile__ ("fnstcw %0" : "=m" (control));
#endif
    return control;
}



/* set status */
___INLINE void setFPUStateX86(uint32_t control)
{
#if defined(_MSVC)
    __asm fldcw control;
#elif defined(__GNUG__)
    __asm__ __volatile__ ("fldcw %0" : : "m" (control));
#endif
}


___INLINE void modifyFPUStateX86(const uint32_t control, const uint32_t mask)
{
    unsigned int oldControl = getFPUStateX86();
    unsigned int newControl = ((oldControl & (~mask)) | (control & mask));
    setFPUStateX86(newControl);
}


___INLINE void setFPUModDefault(void)
{
    modifyFPUStateX86(__MOD_FPU_CW_DEFAULT__, __FPU_CW_MASK_ALL__);
    assertFPUModDefault();
}


___INLINE void assertFPUModDefault(void)
{
#if defined(__GNUG__)
    assert((getFPUStateX86() & (__FPU_CW_MASK_ALL__)) == 
            __MOD_FPU_CW_DEFAULT__);
#endif
}

// taken from http://www.stereopsis.com/FPU.html
// this assumes the CPU is in double precision mode
___INLINE int32_t FastFtol(const float a)
{
    static int32_t    b;
    
#if defined(_MSVC)
    __asm fld a
    __asm fistp b
#elif defined(__GNUG__)
    // use AT&T inline assembly style, document that
    // we use memory as output (=m) and input (m)
    __asm__ __volatile__ (
        "flds %1        \n\t"
        "fistpl %0      \n\t"
        : "=m" (b)
        : "m" (a));
#endif
    return b;
}


//////////////////////////////////////////////////////////////////////

void test()
{
    float testFloats[] = { 1.8, 1.8, 1.1, -1.8, -1.1 };
    int32_t testInt;
    uint32_t i;
    uint32_t oldControl = getFPUStateX86();
    for (i = 0; i < 5; i++) {
        float curTest = testFloats[i];

        setFPUStateX86(oldControl);

        testInt = (int32_t) curTest;
        printf("ANSI slow, default FPU, float %.5f int %d\n",
                curTest, testInt);
        testInt = FastFtol(curTest);
        printf("fast, default FPU, float %.5f int %d\n",
                curTest, testInt);

        setFPUModDefault();
        
        testInt = (int32_t) curTest;
        printf("ANSI slow, modified FPU, float %.5f int %d\n",
                curTest, testInt);
        testInt = FastFtol(curTest);
        printf("fast, modified FPU, float %.5f int %d\n\n",
                curTest, testInt);
    }
}



