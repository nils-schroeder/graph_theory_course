import networkx as nx
from numpy.linalg import matrix_rank


def vertex_count(G: nx.Graph) -> int:
    return G.number_of_nodes()


def edge_count(G: nx.Graph) -> int:
    return G.number_of_edges()


def vertex_degrees(G: nx.Graph) -> tuple[int, ...]:
    degree_list = [v for _, v in G.degree()]
    degree_list.sort()
    return tuple(degree_list)


def algebraic_connectivity(G: nx.Graph) -> float:
    return nx.algebraic_connectivity(G, seed=1337)


def rank(G: nx.Graph) -> int:
    adjacency_matrix = nx.to_numpy_array(G)
    return int(matrix_rank(adjacency_matrix))


def weisfeiler_lehman_graph_hash(G: nx.Graph, iterations: int = 3) -> str:
    for node, data in G.nodes(data=True):
        charge = data.get('charge', '')
        element = data.get('element', '')
        G.nodes[node]['combined'] = f"{charge}{element}"

    return nx.weisfeiler_lehman_graph_hash(
        G,
        edge_attr='order',
        node_attr='combined',
        iterations=iterations
    )


def node_and_edge_label_histogram(G: nx.Graph) -> tuple[tuple[str, str, tuple[str, ...]], ...]:
    histogram = []
    for node, node_data in G.nodes(data=True):
        charge = node_data.get('charge', '')
        element = node_data.get('element', '')
        edge_labels = tuple(sorted(
            [edge_data.get('order', '') for edge, edge_data in G.edges(node, data=True)]))
        histogram.append((charge, element, edge_labels))
    return tuple(sorted(histogram))
