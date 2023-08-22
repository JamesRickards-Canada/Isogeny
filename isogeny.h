/*isogeny.c*/

/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/
GEN getss(GEN p);
long ss_count(GEN p);
GEN ss_graph(GEN p, GEN l);
GEN ss_graphadjmat(GEN p, GEN l);
GEN ss_nbrs(GEN jval, GEN l, GEN pol);

/*SECTION 2: MODULAR POLYNOMIALS*/
GEN modpol(GEN n);
void modpol_processraw(GEN n);