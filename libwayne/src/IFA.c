/* FA accesses a 1D FORTRAN array using FORTRAN semantics, ie, from 1
 * If you compiled your FORTRAN code using -r8, then integers take up
 * 64 bits, although half of them are not used. (Repeat after me:
 * "I love FORTRAN, honest I do!")
 */

#ifdef sgi
#define R8 1    /* 1 for no, 2 for yes if you compiled f77 -r8 */
#else
#define R8 1
#endif

#define IFA(name, index) (*((name) + R8*((index)-1)))   /* integer array */
#define RFA(name, index) (*((name) + (index)-1))        /* Real array */

