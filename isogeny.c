/*Supersingular elliptic curve graph*/

/*INCLUSIONS*/
#include <pari.h>
#include "isogeny.h"

/*STATIC DECLARATIONS*/

/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/
static GEN ssl_graph_i(GEN p, GEN l);
static GEN ssl_graph_givenjvals(GEN p, GEN l, GEN jvals);
static GEN ssl_graph_L(GEN p, GEN L);

/*SECTION 2: MODULAR POLYNOMIALS*/
static GEN getmodpol(GEN l);


/*MAIN BODY*/


/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/

/*Returns a finite field element that is the j-value of a supersingular elliptic curve over F_p.*/
GEN
getssl(GEN p)
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

/*Returns the number of supersingular elliptic curves over F_p.*/
long
ssl_count(GEN p)
{
  long r, pm1 = itos(p) - 1;
  long d = sdivss_rem(pm1, 12, &r);/*d=floor(p-1)/12*/
  if (r % 4 == 2) d++;/*p==3(4), add 1 for j=1728*/
  if (r % 3 == 1) d++;/*p==2(3), add 1 for j=0*/
  return d;
}

/*Returns the supersingular isogeny graph. The output is [v, G], where v is the vector of possible j-invariants (as finite field elements), and G is the vector of Vecsmall of indices of where the ith element of v has directed arrows towards.*/
GEN
ssl_graph(GEN p, GEN l, GEN jvals)
{
  pari_sp av = avma;
  if (typ(l) == t_VEC) return gerepilecopy(av, ssl_graph_L(p, l));/*Vector of L's*/
  if (jvals) return gerepilecopy(av, ssl_graph_givenjvals(p, l, jvals));/*Given jvals*/
  return gerepilecopy(av, ssl_graph_i(p, l));/*Generic, single l and no precomputed jvals.*/
}

/*ssl_graph but no garbage collection.*/
static GEN
ssl_graph_i(GEN p, GEN l)
{
  GEN jval = getssl(p);
  GEN pol = getmodpol(l);
  long nssl = ssl_count(p), vind = 1, lgnbrs = ssl_regularity(l) + 1, i;
  for (i = 1; i <= lgnbrs; i++) {/*Reduce pol modulo p*/
    long lgterm = lg(gel(pol, i)), j;
    for (j = 1; j < lgterm; j++) gmael(pol, i, j) = modii(gmael(pol, i, j), p);
  }
  GEN v = cgetg(nssl + 1, t_VEC);/*Tracks the j-invariants*/
  GEN G = cgetg(nssl + 1, t_VEC);/*Tracks the indices that they go to*/
  for (i = 1; i <= nssl; i++) gel(G, i) = cgetg(lgnbrs, t_VECSMALL);/*(lgnbrs-1)-regular*/
  gel(v, 1) = jval;
  hashtable locs;/*Tracks the found j-invariants and their location in v.*/
  hash_init_GEN(&locs, nssl, &FF_equal, 1);/*Initialize the hash.*/
  hash_insert(&locs, (void *)jval, (void *)1);/*First j-value has index 1.*/
  long maxdepth = 100;/*Maximal depth, to start.*/
  GEN swaps = const_vecsmall(maxdepth, 0);/*Tracks the sequence of swaps.*/
  GEN depthseqlocs = const_vecsmall(maxdepth, 0);/*Tracks the indices in v of the depth sequence.*/
  depthseqlocs[1] = 1;
  GEN depthseq = cgetg(maxdepth + 1, t_VEC);/*Tracks the sequence of j-values found.*/
  gel(depthseq, 1) = jval;
  gel(depthseq, 2) = ssl_nbrs(jval, l, pol);/*The starting neighbours.*/
  long ind = 2;/*Tracks the depth we are working on.*/
  while (ind > 1) {/*We are finding the neighbours of the next j-value.*/
    long cind = ++swaps[ind];/*Increment the swapping index. We insert this (if required) and update the previous term in G to go there. If we did insert this anew, we also compute its neighbours and move deeper; otherwise we already did this before.*/
    if (cind == lgnbrs) {/*Overflowed, go back.*/
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
    gel(depthseq, ind) = ssl_nbrs(curj, l, pol);
  }
  hash_destroy(&locs);/*Done with the hashtable*/
  for (i = 1; i <= nssl; i++) vecsmall_sort(gel(G, i));
  return mkvec2(v, G);
}

/*Given all the supersingular j-values, this gives the supersingular isogeny graph. No garbage collection, not gerepile safe*/
static GEN
ssl_graph_givenjvals(GEN p, GEN l, GEN jvals)
{
  long lgnbrs = ssl_regularity(l) + 1, lj, i, j;
  GEN pol = getmodpol(l);
  for (i = 1; i <= lgnbrs; i++) {/*Reduce pol modulo p*/
    long lgterm = lg(gel(pol, i));
    for (j = 1; j < lgterm; j++) gmael(pol, i, j) = modii(gmael(pol, i, j), p);
  }
  GEN G = cgetg_copy(jvals, &lj);
  hashtable locs;/*Tracks the found j-invariants and their location in v.*/
  hash_init_GEN(&locs, lj, &FF_equal, 1);/*Initialize the hash.*/
  for (i = 1; i < lj; i++) {
	hash_insert(&locs, (void *)gel(jvals, i), (void *)i);/*Insert the indices*/
	gel(G, i) = cgetg(lgnbrs, t_VECSMALL);/*Make the entries.*/
  }
  for (i = 1; i < lj; i++) {
    GEN nbrs = ssl_nbrs(gel(jvals, i), l, pol);
    for (j = 1; j < lgnbrs; j++) {
      hashentry *entry = hash_search(&locs, (void *)gel(nbrs, j));
      gel(G, i)[j] = (long) entry -> val;
    }
  }
  hash_destroy(&locs);/*Done with the hashtable*/
  for (i = 1; i < lj; i++) vecsmall_sort(gel(G, i));
  return mkvec2(jvals, G);
}

/*Does ssl_graph for a vector L.*/
static GEN
ssl_graph_L(GEN p, GEN L)
{
  long lL, i, lgnbrs;
  GEN nbrvec = cgetg_copy(L, &lL);/*Stores the neighbours for each l in L.*/
  if (lL == 1) pari_err_TYPE("L must be non-empty", L);
  GEN G0 = ssl_graph_i(p, gel(L, 1));
  GEN jvals = gel(G0, 1);
  gel(nbrvec, 1) = gel(G0, 2);
  lgnbrs = lg(gmael(G0, 2, 1)) - 1;/*Tracks the total number of neighbours*/
  for (i = 2; i < lL; i++) {
    gel(nbrvec, i) = gel(ssl_graph_givenjvals(p, gel(L, i), jvals), 2);
    lgnbrs += (lg(gmael(nbrvec, i, 1)) - 1);
  }
  lgnbrs++;
  long ljvals;
  GEN allnbrs = cgetg_copy(jvals, &ljvals);
  for (i = 1; i < ljvals; i++) {/*concatenate them all.*/
    long j, k, ind = 0;
    gel(allnbrs, i) = cgetg(lgnbrs, t_VECSMALL);
    for (j = 1; j < lL; j++) {
      long lnb = lg(gmael(nbrvec, j, i));
      for (k = 1; k < lnb; k++) {
      ind++;
      gel(allnbrs, i)[ind] = gmael(nbrvec, j, i)[k];
      }
    }
    vecsmall_sort(gel(allnbrs, i));
  }
  return mkvec2(jvals, allnbrs);
}

/*Writes the adjacency matrix to a file readable as a csr_array in scipy.sparse. Can input p=ssl_graph(p, l) if desired. Can also pass the list of jvals, in case you want to keep the same ordering as another graph. l can be a vector. The file is stored in scipy_adj/p_l.dat. If l is a vector, there are dashes between the l values.*/
void
ssl_graph_scipy(GEN p, GEN l, GEN jvals)
{
  pari_sp av = avma;
  GEN G;
  if (typ(p) == t_VEC) G = gel(p, 2);/*Graph supplied*/
  else if (typ(l) == t_VEC) G = gel(ssl_graph_L(p, l), 2);/*Vector of L's*/
  else if (jvals) G = gel(ssl_graph_givenjvals(p, l, jvals), 2);/*Supplied j values*/
  else G = gel(ssl_graph_i(p, l), 2);/*Regular input*/
  long nrows = lg(G) - 1;
  long reg = lg(gel(G, 1)) - 1;
  long maxnonzero = nrows * reg + 1, i, j;
  GEN rows = vecsmalltrunc_init(maxnonzero);/*Converting to csr format.*/
  GEN cols = vecsmalltrunc_init(maxnonzero);
  GEN data = vecsmalltrunc_init(maxnonzero);
  for (i = 1; i <= nrows; i++) {
	GEN r = gel(G, i);
	long ct = 1;
	for (j = 2; j <= reg; j++) {
	  if (r[j] != r[j - 1]) {
      vecsmalltrunc_append(rows, i);
      vecsmalltrunc_append(cols, r[j - 1]);
      vecsmalltrunc_append(data, ct);
      ct = 1;
	  }
	  else ct++;
	}
	vecsmalltrunc_append(rows, i);
	vecsmalltrunc_append(cols, r[reg]);
	vecsmalltrunc_append(data, ct);
  }
  if (!pari_is_dir("scipy_adj")) {/*Make the directory if it doesn't exist.*/
    int s = system("mkdir -p scipy_adj");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY scipy_adj");
  }
  char *fname;
  if (typ(l) == t_VEC) {
    fname = stack_sprintf("./scipy_adj/%Pd_%Pd", p, gel(l, 1));
    long lL = lg(l);
    for (i = 2; i < lL; i++) fname = stack_sprintf("%s-%Pd", fname, gel(l, i));
    fname = stack_sprintf("%s.dat", fname);
  }
  else fname = stack_sprintf("./scipy_adj/%Pd_%Pd.dat", p, l);
  FILE *f = fopen(fname, "w");
  long ndata = lg(rows) - 1;
  for (i = 1; i < ndata; i++) pari_fprintf(f, "%ld ", data[i]);
  pari_fprintf(f, "%ld\n", data[ndata]);
  for (i = 1; i < ndata; i++) pari_fprintf(f, "%ld ", rows[i] - 1);/*Shift by 1 for python convention.*/
  pari_fprintf(f, "%ld\n", rows[ndata] - 1);/*Shift by 1 for python convention.*/
  for (i = 1; i < ndata; i++) pari_fprintf(f, "%ld ", cols[i] - 1);/*Shift by 1 for python convention.*/
  pari_fprintf(f, "%ld", cols[ndata] - 1);/*Shift by 1 for python convention.*/
  fclose(f);
  set_avma(av);
}

/*Returns the adjacency matrix for the supersingular isogeny graph. Can pass p as the output of ssl_graph and l=jvals=NULL if desired. Can also pass in the vector of j-values.*/
GEN
ssl_graphadjmat(GEN p, GEN l, GEN jvals)
{
  pari_sp av = avma;
  GEN G;
  if (!l) G = gel(p, 2);/*p=ssl_graph(p, l)*/
  else if (typ(l) == t_VEC) G = gel(ssl_graph_L(p, l), 2);/*l is a vector*/
  else if (jvals) G = gel(ssl_graph_givenjvals(p, l, jvals), 2);/*jvals supplied*/
  else G = gel(ssl_graph_i(p, l), 2);/*Generic input.*/
  long lenG = lg(G) - 1, i, j;
  long lnbrs = lg(gel(G, 1));
  GEN M = zeromatcopy(lenG, lenG);
  for (i = 1; i <= lenG; i++) {
    GEN row = gel(G, i);
    for (j = 1; j < lnbrs; j++) gcoeff(M, i, row[j]) = addis(gcoeff(M, i, row[j]), 1);
  }
  return gerepilecopy(av, M);
}

/*Returns the neighbours of the supersingular j-invariant j in the l-isogeny graph. pol should be the output of getmodpol, or can be passed as NULL to compute it.*/
GEN
ssl_nbrs(GEN jval, GEN l, GEN pol)
{
  pari_sp av = avma;
  GEN fone = FF_1(jval);/*1 as a t_FFELT*/
  if (!pol) pol = getmodpol(l);
  long lgpol = lg(pol), i;
  GEN f = cgetg(lgpol + 1, t_POL);/*Stores the polynomial.*/
  for (i = 1; i < lgpol; i++) {
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
  GEN allrts = cgetg(lgpol - 1, t_VEC);
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

/*Returns the regularity of the l-isogeny graph, i.e. product of (p+1)*p^(e-1)*/
long
ssl_regularity(GEN l)
{
  pari_sp av = avma;
  GEN fact = factoru(itou(l));
  long i, s = 1, lfact = lg(gel(fact, 1));
  for (i = 1; i < lfact; i++) {
    long p = gel(fact, 1)[i];
    s *= (p + 1);
    s *= upowuu(p, gel(fact, 2)[i] - 1);
  }
  return gc_long(av, s);
}


/*SECTION 2: MODULAR POLYNOMIALS*/

/*Returns the modular polynomial for the given l. The format is [v1, v2, ..., vn], which represents pol1*y^0+...+poln+y^n-1 and vi=[a1,...,am] means poli=a1x^0+...+amx^(m-1).*/
static GEN
getmodpol(GEN l)
{
  char *fname = stack_sprintf("./modpolys/%Pd.dat", l);
  GEN rdat = gp_readvec_file(fname);/*Read the raw data.*/
  return rdat;
}

/*Assuming that there is a presaved file for the modular polynomial Phi_n(x, y), this returns said polynomial.*/
GEN
modpol(GEN n)
{
  pari_sp av = avma;
  char *fname = stack_sprintf("./modpolys/%Pd.dat", n);
  GEN rdat = gp_readvec_file(fname);/*Read the raw data.*/
  long i, j, lrdat = lg(rdat);
  GEN f = gen_0;
  for (i = 1; i < lrdat; i++) {/*y^(i-1) term. This is not the best way of doing it but doesn't matter.*/
    GEN v = gel(rdat, i);
    long lv = lg(v);
    GEN t = gen_0;
    for (j = 1; j < lv; j++) t = gadd(t, gmul(gel(v, j), pol_xn(j - 1, 0)));
    f = gadd(f, gmul(t, pol_xn(i - 1, 1)));
  }
  return gerepilecopy(av, f);
}

/*Takes in the raw data and makes a file amenable to use in getmodpol.*/
void
modpol_processraw(GEN n)
{
  pari_sp av = avma;
  char *fname = stack_sprintf("./modpolys/%Pd_raw.dat", n);
  GEN rdat = gp_readvec_file(fname);/*Read the raw data.*/
  FILE *f = fopen(stack_sprintf("./modpolys/%Pd.dat", n), "w");
  long i, lrdat = lg(rdat);
  GEN pol = gen_0;/*Stores the modular polynomial.*/
  for (i = 1; i < lrdat; i++) {
    GEN pair = gmael(rdat, i, 1);
    long d1 = itos(gel(pair, 1)), d2 = itos(gel(pair, 2));
    GEN toadd = gmul(pol_xn(d1, 0), pol_xn(d2, 1));
    if (d1 != d2) toadd = gadd(toadd, gmul(pol_xn(d1, 1), pol_xn(d2, 0)));/*Symmetrize it*/
    pol = gadd(pol, gmul(toadd, gmael(rdat, i, 2)));
  }
  long deg = poldegree(pol, 0);
  for (i = 0; i <= deg; i++) {
    GEN term = polcoef(pol, i, 0);/*Term x^i.*/
    long termdeg = poldegree(term, 1), j;
    GEN v = cgetg(termdeg + 2, t_VEC);
    for(j = 0; j <= termdeg; j++) gel(v, j + 1) = polcoef(term, j, 1);/*Make the vector*/
    pari_fprintf(f, "%Ps\n", v);
  }
  fclose(f);
  set_avma(av);
}


