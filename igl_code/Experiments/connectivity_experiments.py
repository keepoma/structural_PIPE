import networkx as nx
from networkx.algorithms.community import girvan_newman
from networkx.algorithms.community.quality import modularity


def compute_metrics(G):
    """
    Compute graph metrics for the given graph G,
    returns a dictionary with all metrics
    """

    metrics = {}

    # Weighted degree
    """
    Sum of the weights of all edges connected to a node
    """
    metrics['node_strength'] = dict(G.degree(weight='weight'))

    """
    Clustering coefficient = how many connections exist among a node’s neighbors compared to 
    the maximum number of possible connections.
    Average Clustering Coefficient is the mean of all individual nodes’ clustering coefficients in the network.
    Higher value = more local specialization
    """
    # Clustering coefficients per node and average clustering
    metrics['clustering'] = nx.clustering(G, weight='weight')
    metrics['avg_clustering'] = nx.average_clustering(G, weight='weight')

    # Global efficiency
    """
    It's the average inverse shortest path length in the network.
    A higher global efficiency indicates that, on average any node can 
    reach any other node through relatively few edges
    This allows us to estimate INTEGRATION
    "integration refers to the ability to quickly combine information from distant brain regions, and segregation refers
     to the ability for specialized processing to occur within dense regions" 
    """
    metrics['global_efficiency'] = nx.global_efficiency(G)

    # Average shortest path length (if graph is disconnected, use the largest connected component)
    """
    The mean of the shortest path lengths between all pairs of nodes in the network.
    A smaller value means that, on average, any pair of nodes can be reached in fewer steps
    """
    if nx.is_connected(G):
        metrics['avg_path_length'] = nx.average_shortest_path_length(G, weight='weight')
    else:
        largest_cc = max(nx.connected_components(G), key=len)
        G_cc = G.subgraph(largest_cc)
        metrics['avg_path_length'] = nx.average_shortest_path_length(G_cc, weight='weight')

    # Degree centrality
    metrics['degree_centrality'] = nx.degree_centrality(G)

    # Eigenvector centrality
    """
    Eigenvector centrality measures a node’s importance by also considering the importance of its neighbors.
    Whereas degree/strength looks only at how many or how strong the direct connections of a node are, 
    eigenvector centrality also weights who those neighbors are
    """
    metrics['eigenvector'] = nx.eigenvector_centrality_numpy(G, weight='weight')

    # Betweenness centrality
    """
    Quantifies how often the node lies on the shortest paths between all other node pairs in the network
    COMPUTATIONALLY INTENSIVE
    """
    metrics['betweenness'] = nx.betweenness_centrality(G, weight='weight')

    # Degree assortativity coefficient
    """
    Measures the correlation between the degrees of nodes that are connected by an edge.
    Ranges from -1 to 1
    A positive value means high-degree nodes tend to connect to other high-degree nodes (assortative).
    A negative value means high-degree nodes tend to connect to low-degree nodes (disassortative)
    """
    metrics['assortativity'] = nx.degree_assortativity_coefficient(G, weight='weight')

    # Modularity: Community detection using Girvan-Newman (take first level partition)
    """
    A measure of how well a network is divided into “communities” or modules, 
    where nodes in the same module are more densely interconnected with each other than with nodes in other modules
    A value near 1 means strong, well-defined communities; near 0 means little to no community structure.
    Negative or extremely low values suggest that the detected partitioning does not yield clear modules (very interconnected)
    """
    communities_gen = girvan_newman(G)
    first_partition = next(communities_gen)
    # Convert communities (tuple of sets) to list of lists
    partition = [list(c) for c in first_partition]
    metrics['communities'] = partition
    # Compute modularity of this partition
    metrics['modularity'] = modularity(G, first_partition, weight='weight')

    return metrics



def print_top_nodes(metrics, lookup, top_n=5):
    """
    Print the top nodes sorted by weighted degree, betweenness, and eigenvector centrality
    """

    centrality = metrics.get('degree_centrality', {})
    node_strength = metrics.get('node_strength', {})
    betweenness = metrics.get('betweenness', {})
    eigenvector = metrics.get('eigenvector', {})

    sorted_centrality = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
    sorted_strength = sorted(node_strength.items(), key=lambda x: x[1], reverse=True)
    sorted_betweenness = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)
    sorted_eigenvector = sorted(eigenvector.items(), key=lambda x: x[1], reverse=True)

    print("\nTop {} nodes by degree centrality:".format(top_n))
    for node, value in sorted_centrality[:top_n]:
        label = lookup.get(node, f"Node_{node}")
        print(f"{label} (ID: {node + 1}): {value}")

    print("\nTop {} nodes by weighted degree (node strength):".format(top_n))
    for node, value in sorted_strength[:top_n]:
        label = lookup.get(node, f"Node_{node}")
        print(f"{label} (ID: {node + 1}): {value}")

    print("\nTop {} nodes by betweenness centrality:".format(top_n))
    for node, value in sorted_betweenness[:top_n]:
        label = lookup.get(node, f"Node_{node}")
        print(f"{label} (ID: {node + 1}): {value:.4f}")

    print("\nTop {} nodes by eigenvector centrality:".format(top_n))
    for node, value in sorted_eigenvector[:top_n]:
        label = lookup.get(node, f"Node_{node}")
        print(f"{label} (ID: {node + 1}): {value:.4f}")


def compute_all():
    # Laptop files
    #sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/Processed/atlas/hcpmmp1.csv"
    #lookup_txt = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"

    # Desktop files
    sc_path = "/home/nikita/Nikita_MRI/me/atlas/hcpmmp1.csv"
    lookup_txt = "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt"

    # Build lookup dictionary for node labels
    lookup = lookup_dictionary(lookup_txt)

    # Load connectivity matrix and create graph.
    sc_matrix = threshold_matrix_by_clipping(sc_path)
    G = create_graph(sc_matrix)

    # Compute all graph metrics.
    metrics = compute_metrics(G)

    # Print computed global metrics.
    print("Average clustering coefficient:", metrics['avg_clustering'])
    print("Global efficiency:", metrics['global_efficiency'])
    print("Average shortest path length:", metrics['avg_path_length'])
    print("Degree assortativity coefficient:", metrics['assortativity'])
    print("Modularity (Girvan-Newman partition):", metrics['modularity'])
    print("Number of communities detected:", len(metrics['communities']), "\n")

    # Print top nodes for several metrics.
    print_top_nodes(metrics, lookup, top_n=5)

