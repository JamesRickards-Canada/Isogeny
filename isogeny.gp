parigp_version=version();
isogeny_library=strprintf("./libisogeny-%d-%d.so", parigp_version[1], parigp_version[2]);

install(getss,"G",,isogeny_library);
addhelp(getss,"getss(p): returns a finite field element of the j-value of a possible supersingular elliptic curve over F_p.");
install(getmodpol,"G");
addhelp(getmodpol,"getmodpol(l): returns the modular polynomial for l. The format is [v1, v2, ..., vn], which represents pol1*y^0+...+poln+y^n-1 and vi=[a1,...,am] means poli=a1x^0+...+amx^(m-1).");
install(ss_count,"lG");
addhelp(ss_count,"ss_count(p): returns the number of supersingular elliptic curves over F_p.");
install(ss_graph,"GG");
addhelp(ss_graph,"ss_graph(p, l): Returns the supersingular isogeny graph G(p, l). The output is [v, G], where v is the vector of possible j-invariants (as finite field elements), and G is the vector of Vecsmall of indices of where the ith element of v has directed arrows towards.");
install(ss_graphadjmat,"GDG");
addhelp(ss_graphadjmat,"ss_graphadjmat(p, {l}): Returns the adjacency matrix of the supersingular isogeny graph G(p, l); can also pass in p=ss_graph(p, l) only.");
install(ss_nbrs,"GGDG");
addhelp(ss_nbrs,"ss_nbrs(jval, l, {pol}): returns the j-values of the neighbours of jval in the l-isogeny graph.");


default(parisize, "4096M");\\Must come at the end