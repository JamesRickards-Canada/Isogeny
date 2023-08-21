parigp_version=version();
isogeny_library=strprintf("./libisogeny-%d-%d.so", parigp_version[1], parigp_version[2]);

install(getss,"G",,isogeny_library);
addhelp(getss,"getss(p): returns a finite field element of the j-value of a possible supersingular elliptic curve over F_p.");


default(parisize, "4096M");\\Must come at the end