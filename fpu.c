/* David Kirkby 21st August 2010
Licenced under the GPL version 2 or at your option any later version.

Set the FPU's precision control to 53 bits (double precision) instead of
the default 64-bits (extended precision). I've commented it fairly
liberally, with the hope it's helpful if anyone needs to edit it.

Note, the extended precision uses 80 bits in total,
of which 64 are for the mantissa.

Double precsion uses 64 bits in total, but ony 53 are
for the mantissa.

The precision is set by bits 8 and 9 of the Control Word in the floating
point processor. The Control word has 16 bits.

Data taken from The 80387 Programmer's reference Manual, Intel,
1987.

00 = 24-bits (single precision)
01 = reserved (or at least it was at the time the 387 was released)
10 = 53-bits (double precision)
11 = 64-bits (extended precision).

FLDCW is an x86 instruction to "Load the Control Word"
FNSTCW is an x88 instruction to "Store FPU Control Word"
It does so without checking for pending unmasked floating-point
exceptions. (A similar FSTCW checks for them first).
*/

/*
 * Rrefreshed and revisited for Debian on behalf of the Debian Science Team
 * by Jerome Benoit <calculus@rezozer.net>, 2014-10-07.
 */

#if defined(ISOC99_FENV)
#include <fenv.h>
#elif defined(FPUCONTROLH)
#include <fpu_control.h>
#define ADHOC__FPU_OR_MASK_EXTENDED (((fpu_control_t)(0x1) << 8) | ((fpu_control_t)(0x1) << 9))
#define ADHOC__FPU_AND_MASK_DOUBLE (~((fpu_control_t)(0x1) << 8))
#elif defined(x86)
#define _SET_FPU_CONTROL_WORD(x) asm volatile ("fldcw %0": :"m" (x));
#define _READ_FPU_CONTROL_WORD(x) asm volatile ("fnstcw %0":"=m" (x));
#else
#endif

void fpu_53bits()
{

#if defined(ISOC99_FENV)
    fesetprec(FE_DBLPREC);
#elif defined(FPUCONTROLH)
    fpu_control_t fpu_control=_FPU_DEFAULT;
    _FPU_GETCW(fpu_control);
    fpu_control|=ADHOC__FPU_OR_MASK_EXTENDED;
    fpu_control&=ADHOC__FPU_AND_MASK_DOUBLE;
		_FPU_SETCW(fpu_control);
#elif defined(x86)
    /* The control word is 16 bits, numbered 0 to 15 */
    volatile unsigned short control_word;

    _READ_FPU_CONTROL_WORD(control_word); /* Read the FPU control word */
    control_word=control_word & 0xfeff; /* Set bit 8 = 0 */
    control_word=control_word | 0x200; /* Set bit 9 = 1 */
    _SET_FPU_CONTROL_WORD(control_word); /* Force double-precision, 53-bit mantissa */
#endif

}
