import time
from collections import defaultdict
from typing import Callable, Any

import networkx as nx
import pandas as pd


def extract_reaction_centers(
        data: list,
        r_id_key: str,
        l_depth: int = 0
) -> dict[str, nx.Graph]:
    reaction_centers = dict()
    for reaction in data:
        r_id = reaction[r_id_key]
        its_graph = reaction['ITS']
        reaction_center = extract_single_reaction_center(its_graph, l_depth)
        reaction_centers[r_id] = reaction_center
    return reaction_centers


def extract_single_reaction_center(its_graph, l_depth=0):
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


def calculate_invariant_values(
        reaction_centers: dict[str, nx.Graph],
        invariant_functions: dict[str, Callable[[nx.Graph], Any]],
) -> dict[str, float | dict[str, Any]]:
    invariant_values = dict()

    for invariant_name, invariant_func in invariant_functions.items():
        start_time = time.time()
        invariant_values[invariant_name] = {
            "reactions": {}
        }

        for r_id, rc in reaction_centers.items():
            invariant_values[invariant_name]['reactions'][r_id] = invariant_func(rc)

        end_time = time.time()
        invariant_values[invariant_name]['execution_time_ms'] = (end_time - start_time) * 1000

    return invariant_values


def cluster_by_invariant_combinations(
        invariant_values: dict[str, float | dict[str, Any]],
        invariant_combinations: list[list[str]],
) -> dict[str, dict[str, float | dict[Any, list[str]]]]:
    invariant_clusters = dict()

    for invariant_combination in invariant_combinations:
        combined_invariant_name = '+'.join(invariant_combination)

        invariant_clusters[combined_invariant_name] = cluster_by_single_invariant_combination(invariant_values,
                                                                                              invariant_combination)

    return invariant_clusters


def cluster_by_single_invariant_combination(
        invariant_values: dict[str, float | dict[str, Any]],
        invariant_combination: list[str],
) -> dict[str, float | dict[Any, list[str]]]:
    start_time = time.time()
    total_execution_time_ms = 0
    invariant_aggregation = defaultdict(list)
    clusters = defaultdict(list)

    for invariant_name in invariant_combination:
        total_execution_time_ms += invariant_values[invariant_name]['execution_time_ms']

        for r_id, invariant_value in invariant_values[invariant_name]['reactions'].items():
            invariant_aggregation[r_id].append(invariant_values[invariant_name]['reactions'][r_id])

    for r_id, invariant_value in invariant_aggregation.items():
        invariant_tuple = tuple(invariant_value)
        clusters[invariant_tuple].append(r_id)

    end_time = time.time()
    total_execution_time_ms += (end_time - start_time) * 1000
    cluster_data = {
        'execution_time_ms': total_execution_time_ms,
        'clusters': clusters
    }
    return cluster_data


def benchmark_invariant_clusters(
        invariant_clusters: dict[str, dict[str, float | dict[Any, list[str]]]],
) -> pd.DataFrame:
    cluster_benchmarks = dict()

    for invariant_name, invariant_cluster_data in invariant_clusters.items():
        clusters = invariant_cluster_data['clusters']
        execution_time_ms = invariant_cluster_data['execution_time_ms']

        cluster_variance = 0
        for cluster in clusters.values():
            cluster_variance += len(cluster) ** 2

        cluster_benchmarks[invariant_name] = {
            "execution_time_ms": execution_time_ms,
            "num_clusters": len(clusters),
            "num_clusters/ms": len(clusters) / execution_time_ms,
            "variance": cluster_variance,
            "variance*ms": cluster_variance * execution_time_ms
        }

    return pd.DataFrame.from_dict(cluster_benchmarks).T


def cluster_invariant_clusters_by_isomorphism(
        reaction_centers: dict[str, nx.Graph],
        invariant_clusters: dict[str, dict[str, float | dict[Any, list[str]]]],
        include_no_invariant: bool = True,
) -> dict[str, dict[str, int | float | dict[str, list[str]]]]:
    isomorphism_clusters = dict()

    for invariant_name, invariant_cluster_data in invariant_clusters.items():
        prefiltered_clusters = invariant_cluster_data["clusters"]
        isomorphism_clusters[invariant_name] = cluster_by_isomorphism(reaction_centers, prefiltered_clusters)

    if include_no_invariant:
        trivial_cluster = {
            "trivial": list(reaction_centers.keys()),
        }
        isomorphism_clusters["no_invariant"] = cluster_by_isomorphism(reaction_centers, trivial_cluster)

    return isomorphism_clusters


def cluster_by_isomorphism(
        reaction_centers: dict[str, nx.Graph],
        prefiltered_clusters: dict[Any, list[str]]
) -> dict[str, float | int | dict[str, list[str]]]:
    start_time = time.time()
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

    end_time = time.time()
    cluster_data = {
        "execution_time_ms": (end_time - start_time) * 1000,
        "isomorphism_check_count": isomorphism_check_count,
        "clusters": final_clusters,
    }
    return cluster_data


def benchmark_isomorphism_clusters(
        isomorphism_clusters: dict[str, dict[str, int | float | dict[str, list[str]]]],
) -> pd.DataFrame:
    cluster_benchmarks = dict()

    for invariant_name, isomorphism_cluster_data in isomorphism_clusters.items():
        clusters = isomorphism_cluster_data['clusters']
        execution_time_ms = isomorphism_cluster_data['execution_time_ms']
        isomorphism_check_count = isomorphism_cluster_data['isomorphism_check_count']

        cluster_benchmarks[invariant_name] = {
            "execution_time_ms": execution_time_ms,
            "isomorphism_check_count": isomorphism_check_count,
            "num_clusters": len(clusters),
        }

    return pd.DataFrame.from_dict(cluster_benchmarks).T
