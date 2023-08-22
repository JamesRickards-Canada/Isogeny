/*isogeny.c*/

/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/
GEN getssl(GEN p);
long ssl_count(GEN p);
GEN ssl_graph(GEN p, GEN l);
GEN ssl_graphadjmat(GEN p, GEN l);
GEN ssl_nbrs(GEN jval, GEN l, GEN pol);

/*SECTION 2: MODULAR POLYNOMIALS*/
GEN modpol(GEN n);
void modpol_processraw(GEN n);