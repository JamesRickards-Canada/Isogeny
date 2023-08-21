/*Supersingular elliptic curve graph*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "isogeny.h"

/*STATIC DECLARATIONS*/

/*MAIN BODY*/


/*SECTION 1: BASIC METHODS*/

/*Returns a finite field element that is the j-value of a supersingular elliptic curve over F_p.*/
GEN
getss(GEN p)
{
  pari_sp av = avma;
  forprime_t T;
  forprime_init(&T, stoi(3), NULL);
  GEN q;
  while ( (q = forprime_next(&T)) ) {
	if (mod4(q) != 3) continue;/*Need q==3 mod 4*/
	if (kronecker(negi(q), p) == -1) break;
  }
  GEN hpol = polclass(negi(q), 0, 0);/*Hilbert class polynomial in the variable x.*/
  GEN g = ffgen(mkvec2(p, gen_2), 1);/*Use the variable y for F_p^2.*/
  GEN rts = FFX_roots(hpol, g);
  if (poldegree(hpol, 0) != lg(rts) - 1) pari_err(e_MISC, "Not enough roots, the Hilbert polynomial should split over F_p^2.");
  return gerepilecopy(av, gel(rts, 1));
}

