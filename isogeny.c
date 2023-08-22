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

/*Returns the number of supersingular elliptic curves over F_p.*/
long
ss_count(GEN p)
{
  long r, pm1 = itos(p) - 1;
  long d = sdivss_rem(pm1, 12, &r);/*d=floor(p-1)/12*/
  if (r % 4 == 2) d++;/*p==3(4), add 1 for j=1728*/
  if (r % 3 == 1) d++;/*p==2(3), add 1 for j=0*/
  return d;
}

/*Returns the supersingular isogeny graph. The output is [v, G], where v is the vector of possible j-invariants (as finite field elements), and G is the vector of Vecsmall of indices of where the ith element of v has directed arrows towards.*/
GEN
ss_graph(GEN p, GEN l)
{
  pari_sp av = avma;
  GEN jval = getss(p);
  GEN pol = getmodpol(l);
  long nss = ss_count(p), vind = 1, i, lp2 = itos(l) + 2;
  GEN v = cgetg(nss + 1, t_VEC);/*Tracks the j-invariants*/
  GEN G = cgetg(nss + 1, t_VEC);/*Tracks the indices that they go to*/
  for (i = 1; i <= nss; i++) gel(G, i) = cgetg(lp2, t_VECSMALL);/*l+1-regular*/
  gel(v, 1) = jval;
  hashtable locs;/*Tracks the found j-invariants and their location in v.*/
  hash_init_GEN(&locs, nss, &FF_equal, 1);/*Initialize the hash.*/
  hash_insert(&locs, (void *)jval, (void *)1);/*First j-value has index 1.*/
  long maxdepth = 100;/*Maximal depth, to start.*/
  GEN swaps = const_vecsmall(maxdepth, 0);/*Tracks the sequence of swaps.*/
  GEN depthseqlocs = const_vecsmall(maxdepth, 0);/*Tracks the indices in v of the depth sequence.*/
  depthseqlocs[1] = 1;
  GEN depthseq = cgetg(maxdepth + 1, t_VEC);/*Tracks the sequence of j-values found.*/
  gel(depthseq, 1) = jval;
  gel(depthseq, 2) = ss_nbrs(jval, l, pol);/*The starting neighbours.*/
  long ind = 2;/*Tracks the depth we are working on.*/
  while (ind > 1) {/*We are finding the neighbours of the next j-value.*/
	long cind = ++swaps[ind];/*Increment the swapping index. We insert this (if required) and update the previous term in G to go there. If we did insert this anew, we also compute its neighbours and move deeper; otherwise we already did this before.*/
	if (cind == lp2) {/*Overflowed, go back.*/
	  swaps[ind] = 0;
	  ind--;
	  continue;
	}
	GEN curj = gmael(depthseq, ind, cind);/*current j value.*/
	hashentry *entry = hash_search(&locs, (void *)curj);
	if (entry) {/*h appeared before, so we have already done it.*/
	  long where = (long) entry -> val;/*Index it was inserted.*/
	  gel(G, depthseqlocs[ind - 1])[cind] = where;/*Update the previous location value in G.*/
	  continue;
	}
	vind++;
	gel(v, vind) = curj;/*Insert it.*/
	hash_insert(&locs, (void *)curj, (void *)vind);
	gel(G, depthseqlocs[ind - 1])[cind] = vind;/*Update the previous location value in G.*/
	depthseqlocs[ind] = vind;
	ind++;/*We move on.*/
	if (ind == maxdepth) {/*We are going too deep, must increase our allocations*/
      long newdepth = maxdepth << 1;/*Double it.*/
	  swaps = vecsmall_lengthen(swaps, newdepth);
	  for (i = maxdepth + 1; i <= newdepth; i++) swaps[i] = 0;
	  depthseqlocs = vecsmall_lengthen(depthseqlocs, newdepth);
	  depthseq = vec_lengthen(depthseq, newdepth);
      maxdepth = newdepth;
    }
	gel(depthseq, ind) = ss_nbrs(curj, l, pol);
  }
  hash_destroy(&locs);/*Done with the hashtable*/
  for (i = 1; i <= nss; i++) vecsmall_sort(gel(G, i));
  return gerepilecopy(av, mkvec2(v, G));
}

/*Returns the adjacency matrix for the supersingular isogeny graph. Can pass either p and l or the output of ss_graph.*/
GEN
ss_graphadjmat(GEN p, GEN l)
{
  pari_sp av = avma;
  GEN gdat;
  if (!l) gdat = p;
  else gdat = ss_graph(p, l);
  GEN G = gel(gdat, 2);
  long lenG = lg(G) - 1, i, j;
  long lp2 = lg(gel(G, 1));
  GEN M = zeromatcopy(lenG, lenG);
  for (i = 1; i <= lenG; i++) {
	GEN row = gel(G, i);
	for (j = 1; j < lp2; j++) gcoeff(M, i, row[j]) = addis(gcoeff(M, i, row[j]), 1);
  }
  return gerepilecopy(av, M);
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


