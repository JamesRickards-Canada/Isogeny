print("\n\nType '?supersingular' for help.\n\n");
addhelp(supersingular, "This package can be used to compute directed supersingular isogeny graphs.\n Installed methods:\ngetssl, ssl_count, ssl_graph, ssl_graph_scipy, ssl_graphadjmat, ssl_nbrs, ssl_regularity\nmodpol, modpol_processraw");
parigp_version=version();
isogeny_library=strprintf("./libisogeny-%d-%d.so", parigp_version[1], parigp_version[2]);

/*SECTION 1: SUPERSINGULAR ELLIPTIC CURVES*/
	install(getssl,"G",,isogeny_library);
	addhelp(getssl,"getssl(p): returns a finite field element of the j-value of a possible supersingular elliptic curve over F_p.");
	install(ssl_count,"lG");
	addhelp(ssl_count,"ssl_count(p): returns the number of supersingular elliptic curves over F_p.");
	install(ssl_graph,"GGDG");
	addhelp(ssl_graph,"ssl_graph(p, l, {jvals}): returns the supersingular isogeny graph G(p, l). The output is [v, G], where v is the vector of possible j-invariants (as finite field elements), and G is the vector of Vecsmall of indices of where the ith element of v has directed arrows towards. If all the supersingular j-invariants are already known, can pass this in as a vector, and we will retain the same ordering. Can also pass l as a vector of l's, to get the capital L isogeny graph.");
	install(ssl_graph_scipy,"vGGDG");
	addhelp(ssl_graph_scipy,"ssl_graph_scipy(p, l, {jvals}): writes the adjacency matrix of G(p, l) to scipy_adj/p_l.dat, in a format readable by csr_array in scipy.sparse. See eigenvalues.py to easily load this. Can also pass in p=ssl_graph(p, l); alternatively, can supply the list of j-values (useful when comparing with different l's). You can supply l as a vector, in which case the l values are separated by dashes. The file format is three lines of space separated integers, with the first being the data, the second being the rows, and the third being the columns of the corresponding sparse matrix.");
	install(ssl_graphadjmat,"GDGDG");
	addhelp(ssl_graphadjmat,"ssl_graphadjmat(p, {l}, {jvals}): returns the adjacency matrix of the supersingular isogeny graph G(p, l); can also pass in p=ssl_graph(p, l) only, l can be passed as a vector, and can also pass in the vector of all supersingular j-values.");
	install(ssl_nbrs,"GGDG");
	addhelp(ssl_nbrs,"ssl_nbrs(jval, l, {pol}): returns the j-values of the neighbours of jval in the l-isogeny graph.");
	install(ssl_regularity,"lG");
	addhelp(ssl_regularity,"ssl_regularity(l): returns the regularity of a l-isogeny graph, i.e. the sum of the divisors of l.");

/*SECTION 2: MODULAR POLYNOMIALS*/
	install(modpol,"G");
	addhelp(modpol,"modpol(n): assuming that there is a presaved file for the modular polynomial Phi_n(x, y), this returns said polynomial. See modpol_processraw for how to create this file.");
	install(modpol_processraw,"vG");
	addhelp(modpol_processraw,"modpol_processraw(n): takes the raw data file for the modular polynomial Phi_n(x, y), and makes a file accessible by this package for this polynomial (for use in ssl_graph or modpol). To get the raw data file, download one of the files from Drew Sutherland's website (https://math.mit.edu/~drew/ClassicalModPolys.html), and save it as n_raw.dat in the folder modpolys. Then go through and make every line enclosed by square brackets, and replace every space by a comma. Once this is done, you can call this method to build the file. ");

default(parisize, "4096M");\\Must come at the end