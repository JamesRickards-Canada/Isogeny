parigp_version=version();
isogeny_library=strprintf("./libisogeny-%d-%d.so", parigp_version[1], parigp_version[2]);

/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/
	install(getss,"G",,isogeny_library);
	addhelp(getss,"getss(p): returns a finite field element of the j-value of a possible supersingular elliptic curve over F_p.");
	install(ss_count,"lG");
	addhelp(ss_count,"ss_count(p): returns the number of supersingular elliptic curves over F_p.");
	install(ss_graph,"GG");
	addhelp(ss_graph,"ss_graph(p, l): Returns the supersingular isogeny graph G(p, l). The output is [v, G], where v is the vector of possible j-invariants (as finite field elements), and G is the vector of Vecsmall of indices of where the ith element of v has directed arrows towards.");
	install(ss_graphadjmat,"GDG");
	addhelp(ss_graphadjmat,"ss_graphadjmat(p, {l}): Returns the adjacency matrix of the supersingular isogeny graph G(p, l); can also pass in p=ss_graph(p, l) only.");
	install(ss_nbrs,"GGDG");
	addhelp(ss_nbrs,"ss_nbrs(jval, l, {pol}): returns the j-values of the neighbours of jval in the l-isogeny graph.");

/*SECTION 2: MODULAR POLYNOMIALS*/
	install(modpol,"G");
	addhelp(modpol,"modpol(n): assuming that there is a presaved file for the modular polynomial Phi_n(x, y), this returns said polynomial. See modpol_processraw for how to create this file.");
	install(modpol_processraw,"vG");
	addhelp(modpol_processraw,"modpol_processraw(n): takes the raw data file for the modular polynomial Phi_n(x, y), and makes a file accessible by this package for this polynomial (for use in ss_graph or modpol). To get the raw data file, download one of the files from Drew Sutherland's website (https://math.mit.edu/~drew/ClassicalModPolys.html), and save it as n_raw.dat in the folder modpolys. Then go through and make every line enclosed by square brackets, and replace every space by a comma. Once this is done, you can call this method to build the file. ");

default(parisize, "4096M");\\Must come at the end