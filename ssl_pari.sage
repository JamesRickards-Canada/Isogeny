gp.read("isogeny.gp")
from sage.libs.pari.convert_sage import gen_to_sage

def ssl_graph(p, l):
	#Make the graph in PARI/GP
	gpgraph = gp.ssl_graph(p, l)
	jinv = gpgraph[1]
	jinv_conv = pari(jinv)
	edges = pari(gpgraph[2])
	#First, we sort out the labels
	labels = []
	Fp2 = GF(p^2)
	#It's annoying to convert finite fields to Sage.
	F.<x> = ZZ[] #Make a polynomial ring in the integers
	gpgen = gp.ffgen(jinv[1])
	#Let's find the min poly for gpgen, which is the variable used in jinv.
	coefs = [-gen_to_sage(coef) for coef in pari(gpgen^2).list()]
	coefs.append(1) #The min poly coefs for the generator
	mpol = F(coefs).change_ring(Fp2) #The polynomial
	gpgen_inFp2 = mpol.any_root() #We have written it in Fp2=GF(p^2)
	pari_itype = gp.type(0)
	for i in range(len(jinv)):#Let's go through and convert them!
			if gp.type(jinv_conv[i]) == pari_itype:#Integer, i.e. F_p j-invariant.
					labels.append(Fp2(gen_to_sage(jinv_conv[i])))
					continue
			#Not an integer, so must do more work.
			coefs = [gen_to_sage(coef) for coef in jinv_conv[i].list()]
			labels.append(coefs[0]+coefs[1]*gpgen_inFp2)
	#Initialize empty graph, then loop over the edges, adding them in
	G = DiGraph(loops = True, multiedges = True)
	for i in range(len(edges)):
			vtx = edges[i]
			for j in range(len(vtx)):
					edge = vtx[j]
					G.add_edge(labels[i], labels[edge - 1])
	return G	

