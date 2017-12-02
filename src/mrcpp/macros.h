/*
 * 
 * \breif
 */

#pragma once

#define  MAX(X,Y)   (((X) > (Y)) ? (X):(Y))
#define  MIN(X,Y)   (((X) < (Y)) ? (X):(Y))
#define  SQUARE(X)  ((X) * (X))
#define  CUBE(X)    ((X) * (X) * (X))
#define CLEAR_ARRAY(X,N) for (int _i = 0; _i < N; _i++) X[_i] = 0;

#define SET_BIT(A,N) A |= (1 << N)
#define CLEAR_BIT(A,N) A &= ~(1 << N)
#define TOGGLE_BIT(A,N) A ^= (1 << N)
#define TEST_BIT(A,N) {if (A & (1 << N)) return true; return false;}
#define IS_ODD(A) ( A & 1 )
#define IS_EVEN(A) (!( A & 1 ))
#define IS_EQUAL(A,B) (fabs(A - B) < MachineZero)

#define SET_BITS(A,N) A |= (N)
#define CLEAR_BITS(A,N) A &= ~(N)
#define TOGGLE_BITS(A,N) A ^= (N)

/* Binary constant generator macro
 * By Tom Torfs - donated to the public domain
 * All macro's evaluate to compile-time constants
 */

 /* turn a numeric literal into a hex constant
 (avoids problems with leading zeroes)
 8-bit constants max value 0x11111111, always fits in unsigned long
 */
#define HEX__(n) 0x##n##LU

/* 8-bit conversion function */
#define B8__(x) ((x&0x0000000FLU)?1:0)      \
			   +((x&0x000000F0LU)?2:0)      \
			   +((x&0x00000F00LU)?4:0)      \
			   +((x&0x0000F000LU)?8:0)      \
			   +((x&0x000F0000LU)?16:0)     \
			   +((x&0x00F00000LU)?32:0)     \
			   +((x&0x0F000000LU)?64:0)     \
			   +((x&0xF0000000LU)?128:0)

/* *** user macros *** */

 /* for upto 8-bit binary constants */
#define B8(d) ((unsigned char)B8__(HEX__(d)))

/* for upto 16-bit binary constants, MSB first */
#define B16(dmsb,dlsb) (((unsigned short)B8(dmsb)<<8)     \
						+ B8(dlsb))

/* for upto 32-bit binary constants, MSB first */
#define B32(dmsb,db2,db3,dlsb) (((unsigned long)B8(dmsb)<<24)      \
								  + ((unsigned long)B8(db2)<<16) \
								  + ((unsigned long)B8(db3)<<8)    \
								  + B8(dlsb))

/* Sample usage:
 B8(01010101) = 85
 B16(10101010,01010101) = 43605
 B32(10000000,11111111,10101010,01010101) = 2164238933
 */

