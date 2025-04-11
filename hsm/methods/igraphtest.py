import igraph as ig

g = ig.Graph(directed=False)

mode = "seq"

if mode == "seq":
    # Import edges from sequence identity matrix
    with open("hsm/outs/Clustal/seq_edge_list.csv") as f:
        next(f) # Skip header
        edges = [line.strip().split(",") for line in f]
    g.add_vertices(len(edges))
    g.add_edges(edges)
elif mode == "struct":
    pass
