gp.read("isogeny.gp")
from sage.libs.pari.convert_sage import gen_to_sage

#Makes the (sparse) adjacency matrix for the supersingular graph at (p, l)
#Faster than calling ssl_graph(p, l).adjacency_matrix(sparse = True), as we don't make a DiGraph element
def ssl_adjmat(p, l):
    gpgraph = gp.ssl_graph(p, l)
    edges = pari(gpgraph[2])
    nver = len(edges)
    M = matrix(ZZ, nver, sparse = True)
    for i in range(nver):
        vtx = edges[i]
        for j in range(len(vtx)):
            M[i, vtx[j] - 1] += 1
    return M

#Make the supersingular isogeny graph in PARI/GP, and import into Sage.
#We are careful to try to limit calls to gp. and gen_to_sage, as they are EXTREMELY SLOW.
#Seriously, they made an initial version over 150 times as slow.
def ssl_graph(p, l):
    gpgraph = gp.ssl_graph(p, l)
    jinv = pari(gpgraph[1])
    edges = pari(gpgraph[2])
    #First, we sort out the labels
    labels = []
    Fp2 = GF(p^2)
    #It's annoying to convert finite fields to Sage.
    F.<y> = ZZ[] #Make a polynomial ring in the integers
    #Relies on ffgen always outputting the same generator. We also use y, the same variable as ssl_graph.
    gpgen = pari(p^2).ffgen(y)
    #Let's find the min poly for gpgen=y, which is the variable used in jinv.
    mpol = y^2-eval(gen_to_sage((gpgen^2).Str()))
    mpol_Fp2 = mpol.change_ring(Fp2) #The polynomial
    gpgen_inFp2 = mpol_Fp2.any_root() #We have written it in Fp2=GF(p^2)
    sage_j = eval(gen_to_sage(jinv.Str()))#The j-invariants, as polynomials. The integral ones are OK, the others we must convert.
    notintclass = type(mpol)
    for j in sage_j:#Let's go through and convert them!
        if type(j) != notintclass:#Integer, i.e. F_p j-invariant.
            labels.append(Fp2(j))
            continue
        #Not an integer, so must do more work.
        coefs = j.list()
        labels.append(coefs[0] + coefs[1] * gpgen_inFp2)
    #Initialize empty graph, then loop over the edges, adding them in
    G = DiGraph(loops = True, multiedges = True)
    for i in range(len(edges)):
        vtx = edges[i]
        for j in range(len(vtx)):
            edge = vtx[j]
            G.add_edge(labels[i], labels[edge - 1])
    return G

