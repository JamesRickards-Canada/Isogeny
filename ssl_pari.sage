gp.read("isogeny.gp")

def ssl_graph(p, l):
  #Make the graph in PARI/GP
  x = pari(gp.ssl_graph(p, l))
  labels = [str(v) for v in x[0]]
  #Initialize empty graph, then loop over the edges, adding them in
  G = DiGraph(loops = True, multiedges = True)
  for i in range(len(x[1])):
    vtx = x[1][i]
    for j in range(len(vtx)):
      edge = vtx[j]
      G.add_edge(labels[i], labels[edge - 1])
  return G