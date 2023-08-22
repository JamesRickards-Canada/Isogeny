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

/*CURRENTLY ONLY WORKS FOR l=2*/

/*Returns the modular polynomial for the given l. The format is [v1, v2, ..., vn], which represents pol1*y^0+...+poln+y^n-1 and vi=[a1,...,am] means poli=a1x^0+...+amx^(m-1).*/
GEN
getmodpol(GEN l)
{
  pari_sp av = avma;
  if (!gequal(l, gen_2)) pari_err(e_IMPL, "Not yet implemented.");
  return gerepilecopy(av, mkvec4(mkvec4s(-157464000000000, 8748000000, -162000, 1), mkvec3s(8748000000, 40773375, 1488), mkvec3s(-162000, 1488, -1), mkvec(gen_1)));
}

/*Returns the neighbours of the supersingular j-invariant j in the l-isogeny graph. pol should be the output of getmodpol, or can be passed as NULL to compute it.*/
GEN
ss_nbrs(GEN jval, GEN l, GEN pol)
{
  pari_sp av = avma;
  GEN fone = FF_1(jval);/*1 as a t_FFELT*/
  if (!pol) pol = getmodpol(l);
  long sl = itos(l), i;
  GEN f = cgetg(sl + 4, t_POL);/*Stores the polynomial.*/
  for (i = 1; i <= sl + 2; i++) {
	GEN curcoef = gel(pol, i);/*Polynomial for the current coefficient.*/
	long j = lg(curcoef) - 1;
	GEN c = FF_Z_mul(fone, gel(curcoef, j));/*Convert last element to being in the finite field.*/
	j--;
	for ( ;j > 0; j--) {
	  c = FF_mul(c, jval);
	  c = FF_Z_add(c, gel(curcoef, j));/*Update the evaulation.*/
	}
	gel(f, i + 1) = c;
  }
  setvarn(f, 0);/*Make sure it is in the variable x.*/
  GEN fact = FFX_factor(f, jval);/*Factor it since we need the multiplicities.*/
  GEN allrts = cgetg(sl + 2, t_VEC);
  long nbfact = nbrows(fact), ind = 1;
  for (i = 1; i <= nbfact; i++) {
	GEN rt = FF_neg(polcoef(gcoeff(fact, i, 1), 0, 0));/*The root, we are assuming the polynomial factors (which it does).*/
	long j, po = itos(gcoeff(fact, i, 2));
	for (j = 1; j <= po; j++) {
	  gel(allrts, ind) = rt;
	  ind++;
	}
  }
  return gerepilecopy(av, allrts);
}


