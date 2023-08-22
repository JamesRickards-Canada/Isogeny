parigp_version=version();
isogeny_library=strprintf("./libisogeny-%d-%d.so", parigp_version[1], parigp_version[2]);

install(getss,"G",,isogeny_library);
addhelp(getss,"getss(p): returns a finite field element of the j-value of a possible supersingular elliptic curve over F_p.");
install(getmodpol,"G");
addhelp(getmodpol,"getmodpol(l): returns the modular polynomial for l. The format is [v1, v2, ..., vn], which represents pol1*y^0+...+poln+y^n-1 and vi=[a1,...,am] means poli=a1x^0+...+amx^(m-1).");
install(ss_nbrs,"GGDG");
addhelp(ss_nbrs,"ss_nbrs(jval, l, {pol}): returns the j-values of the neighbours of jval in the l-isogeny graph.");

default(parisize, "4096M");\\Must come at the end