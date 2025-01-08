from typing import Dict, List, Tuple

import networkx as nx


def reaction_center_subgraph(its_graph, l_depth=0):
    rc_edges = set()
    for (u, v, d) in its_graph.edges(data=True):
        if d["standard_order"] != 0:
            rc_edges = rc_edges | {(u, v)}

    rc_subgraph = nx.edge_subgraph(its_graph, rc_edges)

    for _ in range(l_depth):
        for n in rc_subgraph.nodes:
            rc_edges = rc_edges | set(its_graph.edges(n))
        rc_subgraph = nx.edge_subgraph(its_graph, rc_edges)

    return rc_subgraph.copy()


def create_reaction_centers(
        data: list,
        r_id_key: str,
        l_depth: int = 0
) -> Dict[str, nx.Graph]:
    reaction_centers = dict()
    for reaction in data:
        r_id = reaction[r_id_key]
        its_graph = reaction['ITS']
        reaction_center = reaction_center_subgraph(its_graph, l_depth)
        reaction_centers[r_id] = reaction_center
    return reaction_centers


def isomorphism_clustering(
        reaction_centers: Dict[str, nx.Graph],
        prefiltered_clusters: Dict[..., List[str]]
) -> Tuple[Dict[str, List[str]], int]:
    final_clusters = dict()
    isomorphism_check_count = 0

    for prefiltered_cluster in prefiltered_clusters.values():
        new_clusters = dict()

        for r_id in prefiltered_cluster:
            rc = reaction_centers[r_id]
            for cluster_key in new_clusters.keys():
                isomorphism_check_count += 1
                if nx.is_isomorphic(
                        rc,
                        reaction_centers[cluster_key],
                        lambda n1, n2: n1['charge'] == n2['charge'] and n1['element'] == n2['element'],
                        lambda e1, e2: e1['order'] == e2['order']
                ):
                    new_clusters[cluster_key].append(r_id)
                    break
            else:
                new_clusters[r_id] = [r_id]

        final_clusters.update(new_clusters)

    return final_clusters, isomorphism_check_count
