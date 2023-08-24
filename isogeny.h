/*isogeny.c*/

/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/
GEN getssl(GEN p);
long ssl_count(GEN p);
GEN ssl_graph(GEN p, GEN l, GEN jvals);
void ssl_graph_scipy(GEN p, GEN l, GEN jvals);
GEN ssl_graphadjmat(GEN p, GEN l, GEN jvals);
GEN ssl_nbrs(GEN jval, GEN l, GEN pol);
long ssl_regularity(GEN l);

/*SECTION 2: MODULAR POLYNOMIALS*/
GEN modpol(GEN n);
void modpol_processraw(GEN n);