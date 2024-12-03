
import rustworkx as rx


def Connected(model, state):
    # Get the connected components of the graph given the state, and ask whether
    # the vertices of the edge are in the same connected component.
    G = model.lattice.graph
    subgraph = G.copy()
    
    excluded = []
    for edge in model.lattice.cubes:
        if edge not in model.occupied:
            u, v = edge.faces
            excluded.append((u.index, v.index))

    subgraph.remove_edges_from(excluded)

    edge = model.lattice.cubes[13]
    u, v = edge.faces

    return int(rx.has_path(subgraph, u.index, v.index))


def GraphConnected(model, state):
    edge = model.lattice.cubes[13]
    u, v = edge.at

    # Ignore all the edges that aren't occupied.
    G = model.lattice.graph
    subgraph = G.copy()

    excluded = []
    for e in G.edges():
        x, y = e.at
        if e.spin == 0: excluded.append((x.index, y.index))

    subgraph.remove_edges_from(excluded)

    return int(rx.has_path(subgraph, u.index, v.index))
