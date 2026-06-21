"""
Cascade — core algorithms.

Implements, with plain Python + networkx, the classic bioinformatics
algorithms the project is meant to showcase:
  - Graph construction from PPI edges
  - BFS / DFS traversal
  - Shortest path (Dijkstra, unweighted here so it's = BFS shortest path)
  - Clustering (connected components + greedy modularity communities)
  - Enrichment analysis (hypergeometric test, like a simplified GO/KEGG enrichment)
  - Drug target suggestion (lookup over enriched pathways)
"""

import itertools
from collections import deque
from scipy.stats import hypergeom
import networkx as nx

from data import GENE_PATHWAYS, PPI_EDGES, PATHWAY_DRUG_TARGETS, BACKGROUND_GENE_COUNT


def build_graph(genes):
    """Build the subgraph of the known PPI network induced by the uploaded genes
    plus their direct neighbors (so the network isn't just isolated dots)."""
    full_graph = nx.Graph()
    full_graph.add_edges_from(PPI_EDGES)

    genes_in_kb = [g for g in genes if g in full_graph]
    neighbors = set()
    for g in genes_in_kb:
        neighbors.update(full_graph.neighbors(g))

    node_set = set(genes_in_kb) | neighbors
    subgraph = full_graph.subgraph(node_set).copy()
    return subgraph, genes_in_kb


def bfs_traversal(graph, start):
    """Breadth-first traversal from `start`, returns visit order."""
    if start not in graph:
        return []
    visited = {start}
    order = [start]
    q = deque([start])
    while q:
        node = q.popleft()
        for nbr in sorted(graph.neighbors(node)):
            if nbr not in visited:
                visited.add(nbr)
                order.append(nbr)
                q.append(nbr)
    return order


def dfs_traversal(graph, start):
    """Depth-first traversal from `start`, returns visit order."""
    if start not in graph:
        return []
    visited = set()
    order = []

    def _dfs(node):
        visited.add(node)
        order.append(node)
        for nbr in sorted(graph.neighbors(node)):
            if nbr not in visited:
                _dfs(nbr)

    _dfs(start)
    return order


def shortest_paths_between_inputs(graph, genes_in_kb):
    """Shortest path (BFS-based, since the network is unweighted) between
    every pair of uploaded genes that both exist in the network."""
    results = []
    for a, b in itertools.combinations(sorted(set(genes_in_kb)), 2):
        if a in graph and b in graph and nx.has_path(graph, a, b):
            path = nx.shortest_path(graph, a, b)
            results.append({"from": a, "to": b, "length": len(path) - 1, "path": path})
    results.sort(key=lambda r: r["length"])
    return results


def cluster_network(graph):
    """Cluster the network into communities using greedy modularity
    maximization (falls back to connected components for tiny graphs)."""
    if graph.number_of_nodes() < 3 or graph.number_of_edges() == 0:
        components = list(nx.connected_components(graph))
        return [sorted(c) for c in components]

    communities = nx.algorithms.community.greedy_modularity_communities(graph)
    return [sorted(c) for c in communities]


def enrichment_analysis(genes_in_kb):
    """Hypergeometric enrichment test of uploaded genes against each pathway.

    p-value = P(X >= k) where:
      N = background gene pool size
      K = number of genes in that pathway (within our knowledge base)
      n = number of uploaded genes (found in knowledge base)
      k = number of uploaded genes found in that pathway
    """
    pathway_to_genes = {}
    for gene, pathways in GENE_PATHWAYS.items():
        for p in pathways:
            pathway_to_genes.setdefault(p, set()).add(gene)

    n = len(set(genes_in_kb))
    N = BACKGROUND_GENE_COUNT
    results = []
    for pathway, members in pathway_to_genes.items():
        k = len(set(genes_in_kb) & members)
        if k == 0:
            continue
        K = len(members)
        # P(X >= k)
        pval = hypergeom.sf(k - 1, N, K, n)
        results.append({
            "pathway": pathway,
            "genes_hit": sorted(set(genes_in_kb) & members),
            "hit_count": k,
            "pathway_size": K,
            "p_value": round(float(pval), 6),
        })
    results.sort(key=lambda r: r["p_value"])
    return results


def predict_affected_processes(enriched):
    """Very simple rule-based mapping from enriched pathways to higher-level
    biological processes/themes, just to demonstrate the 'predict affected
    biological processes' feature end to end."""
    theme_map = {
        "Apoptosis": "Programmed cell death may be dysregulated",
        "Cell Cycle": "Cell cycle checkpoint control may be disrupted",
        "DNA Damage Response": "Genome stability / DNA repair may be impaired",
        "PI3K-AKT Signaling": "Cell growth & survival signaling may be hyperactive",
        "MAPK Signaling": "Proliferation signaling may be hyperactive",
        "Angiogenesis": "Blood vessel formation may be altered",
        "Inflammation": "Immune/inflammatory response may be altered",
        "Metabolism": "Cellular energy metabolism may be reprogrammed",
    }
    return [
        {"pathway": e["pathway"], "process": theme_map.get(e["pathway"], "Effect unclear from current knowledge base")}
        for e in enriched
    ]


def suggest_drug_targets(enriched):
    suggestions = []
    for e in enriched:
        targets = PATHWAY_DRUG_TARGETS.get(e["pathway"], [])
        for t in targets:
            suggestions.append({"pathway": e["pathway"], "target": t})
    return suggestions


def analyze(genes):
    """Run the full pipeline and return everything the frontend needs."""
    genes = [g.strip().upper() for g in genes if g.strip()]
    graph, genes_in_kb = build_graph(genes)

    not_found = sorted(set(genes) - set(genes_in_kb))

    enriched = enrichment_analysis(genes_in_kb)
    clusters = cluster_network(graph)
    processes = predict_affected_processes(enriched)
    drug_targets = suggest_drug_targets(enriched)
    shortest = shortest_paths_between_inputs(graph, genes_in_kb)

    bfs_order = bfs_traversal(graph, genes_in_kb[0]) if genes_in_kb else []
    dfs_order = dfs_traversal(graph, genes_in_kb[0]) if genes_in_kb else []

    nodes = [{"id": n, "is_input": n in genes_in_kb} for n in graph.nodes()]
    edges = [{"source": u, "target": v} for u, v in graph.edges()]

    return {
        "input_genes": genes,
        "genes_found": genes_in_kb,
        "genes_not_found": not_found,
        "network": {"nodes": nodes, "edges": edges},
        "enrichment": enriched,
        "clusters": clusters,
        "predicted_processes": processes,
        "drug_targets": drug_targets,
        "shortest_paths": shortest,
        "traversal": {
            "start_gene": genes_in_kb[0] if genes_in_kb else None,
            "bfs_order": bfs_order,
            "dfs_order": dfs_order,
        },
    }
