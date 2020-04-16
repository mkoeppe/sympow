/*
 * Copyright Jerome G. M. Benoit <jgmbenoit@rezozer.net> 2020
 *
 */
#if !(__STDC_VERSION__ >= 201112L)
#error "require C11 or higher"
#endif

#include <sys/types.h>
#include <sys/param.h>
#include <stdalign.h>
#include <stdio.h>

int main() {
	const char * ADHOC_ARCH_ENDIAN="";
	const char * ADHOC_ARCH_SIZE_ALIGN="";
	int is_bigendian=-1;

/* inspired by AC_C_BIGENDIAN */
#ifdef __APPLE_CC__
	is_bigendian=0;
#elif (defined BYTE_ORDER && defined BIG_ENDIAN && defined LITTLE_ENDIAN && \
		BYTE_ORDER && BIG_ENDIAN && LITTLE_ENDIAN)
#	if BYTE_ORDER == LITTLE_ENDIAN
		is_bigendian=0;
# elif BYTE_ORDER == BIG_ENDIAN
		is_bigendian=1;
# else
#	error "unexpected event"
#	endif
#elif (defined _LITTLE_ENDIAN || defined _BIG_ENDIAN)
#	ifndef _BIG_ENDIAN
		is_bigendian=0;
# else
		is_bigendian=1;
#	endif
#else
	/* From Harbison&Steele, section 6.1.2, ``Byte Ordering'' */
	union {
		long int l;
		char c[sizeof(long int)];
		} u;
	u.l=1;

	if (u.c[0] == 1) {
		is_bigendian=0;
		}
	else if (u.c[sizeof(long int)-1] == 1) {
		is_bigendian=1;
		}
#endif

	if (!(is_bigendian < 0)) {
		if (is_bigendian == 0) {
			ADHOC_ARCH_ENDIAN="le";
			}
		else if (0 < is_bigendian) {
			ADHOC_ARCH_ENDIAN="be";
			}
		if (sizeof(void *) == 4) {
			if (alignof(double) == 4) {
				ADHOC_ARCH_SIZE_ALIGN="32d4";
				}
			else {
				ADHOC_ARCH_SIZE_ALIGN="32d8";
				}
			}
		else { /* (sizeof(void *) != 4) */
			ADHOC_ARCH_SIZE_ALIGN="64";
			}
		}

	printf("%s%s\n",ADHOC_ARCH_ENDIAN,ADHOC_ARCH_SIZE_ALIGN);

	return 0; }
