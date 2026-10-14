import pandas as pd, networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

nodes = pd.read_csv("cytoscape/nodes.csv")
nodes = nodes[nodes["species"] != "AMR gene"].set_index("id")
e = pd.read_csv("cytoscape/edges_plasmids.csv")

COUNTRY = {"Malawi": "#1b9e77", "SouthAfrica": "#d95f02", "Zambia": "#7570b3"}
TIER = {"local co-occurrence": "#bdbdbd", "persistent reservoir": "#fdae61", "broad regional dispersal": "#d7191c"}
SOURCE = {"clinical": "#000000", "animal_food": "#8c510a", "environmental": "#2166ac", "community_unspecified": "#999999"}
MARK = {"E_coli": "o", "K_pneumoniae": "s", "Salmonella": "^"}

def collapse(df):
    return (df.groupby(["source", "target"], as_index=False)
              .agg(tier=("tier", "first"), cross=("cross_species", "max"), n_plasmids=("plasmid_pair", "nunique")))

def draw(df, name, title):
    c = collapse(df)
    G = nx.Graph()
    for r in c.itertuples():
        G.add_edge(r.source, r.target, tier=r.tier, cross=bool(r.cross))
    pos = nx.spring_layout(G, seed=42, k=0.5, iterations=300)
    fig, ax = plt.subplots(figsize=(13, 10))
    for t, col in TIER.items():
        for cross in (False, True):
            el = [(u, v) for u, v, d in G.edges(data=True) if d["tier"] == t and d["cross"] == cross]
            if el:
                nx.draw_networkx_edges(G, pos, edgelist=el, edge_color=col, width=3.5 if cross else 1, alpha=0.8, ax=ax)
    for sp, mk in MARK.items():
        ns = [n for n in G if nodes.loc[n, "species"] == sp]
        if ns:
            nx.draw_networkx_nodes(G, pos, nodelist=ns, node_shape=mk, node_size=140,
                node_color=[COUNTRY[nodes.loc[n, "country"]] for n in ns],
                edgecolors=[SOURCE[nodes.loc[n, "source_category"]] for n in ns], linewidths=2, ax=ax)
    labelled = {n for u, v, d in G.edges(data=True) if d["cross"] for n in (u, v)}
    short = lambda n: n.replace("K_pneumoniae_", "K_").replace("Salmonella_", "S_").replace("E_coli_", "E_")
    nx.draw_networkx_labels(G, pos, labels={n: short(n) for n in labelled}, font_size=6, ax=ax)
    h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=9, label=k) for k, c in COUNTRY.items()]
    h += [Line2D([0], [0], marker=m, color="w", markerfacecolor="#cccccc", markeredgecolor="k", markersize=9, label=k) for k, m in MARK.items()]
    h += [Line2D([0], [0], marker="o", color="w", markerfacecolor="w", markeredgecolor=c, markeredgewidth=2, markersize=9, label="source: " + k) for k, c in SOURCE.items()]
    h += [Line2D([0], [0], color=c, lw=2, label=k) for k, c in TIER.items()]
    h.append(Line2D([0], [0], color="k", lw=3.5, label="cross-species (thick edge)"))
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1, 1), fontsize=8, frameon=False)
    ax.set_title(title); ax.axis("off")
    plt.savefig(f"cytoscape/{name}.png", dpi=200, bbox_inches="tight")
    plt.savefig(f"cytoscape/{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(name, G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges")

big = e[e["min_plasmid_size"] >= 5000]
draw(big, "network_5kb", "Isolates sharing near-identical plasmids (plasmids >= 5 kb)")
draw(big[big["cross_species"]], "network_cross_species", "Cross-species plasmid sharing (plasmids >= 5 kb)")
