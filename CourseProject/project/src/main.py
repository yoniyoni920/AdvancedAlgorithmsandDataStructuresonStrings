import time, os
from OPM import order_preserving_match
from CTM import discover_ctm_motifs
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import networkx as nx


def plot_top_10_motifs_shapes(data, motifs, length, frequency, save_path, title_prefix):
    """
    Generates a line plot showing the normalized shapes of the top 10 most frequent motifs.

    Parameters:
    -----------
    data : list
        The raw temperature anomaly data.
    motifs : list
        The list of found motifs (OPM or CTM) sorted by frequency.
    length : int
        The length (L) of the motif.
    frequency : str
        'monthly' or 'annual' indicating the time scale.
    save_path : str
        The file path to save the resulting plot.
    title_prefix : str
        Prefix for the chart title (e.g., 'OPM' or 'CTM').

    What the Axes Represent:
    ------------------------
    X-Axis (Time Steps): Represents the relative time progression within the motif.
                         For L=3, steps are 0, 1, 2.
    Y-Axis (Normalized Value): Represents the temperature anomaly scaled between 0 and 1.
                               0 is the lowest point in that specific sequence, 1 is the highest.
                               This allows for structural comparison regardless of absolute temperature.
    """
    if not motifs:
        return

    top_10 = motifs[:10]
    plt.figure(figsize=(12, 7))

    for i, (prefix_ranks, (count, positions, _)) in enumerate(top_10):
        # We take the first occurrence of the motif to represent its shape
        pos = positions[0]
        seq = data[pos:pos + length]

        # Min-Max Normalization to bring all shapes to the same 0-1 scale
        seq_norm = (seq - np.min(seq)) / (np.ptp(seq) + 1e-8)

        plt.plot(range(length), seq_norm, marker='o', linewidth=2,
                 label=f"Rank {i + 1}: {prefix_ranks} ({count}x)")

    plt.title(f"{title_prefix} - Top 10 Most Frequent Motifs \n(L={length}, {frequency})", fontsize=15)
    plt.xlabel("Time Steps within Motif (Relative)", fontsize=12)
    plt.ylabel("Normalized Anomaly Value (0-1)", fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_ctm_tree(motif_tuple, length, count, save_path):
    """
        Visualizes a Cartesian Tree for a given CTM (Cartesian Tree Motif) motif.

        Parameters
        ----------
        motif_tuple : tuple or list
            The motif represented as a sequence or parent-distance vector.
            Each element corresponds to a value in the motif.
        length : int
            Length of the motif (number of elements, L).
        count : int
            The frequency of this motif in the dataset.
        save_path : str
            File path to save the generated tree visualization as a PNG image.

        Description
        -----------
        This function constructs a Cartesian Tree from the input motif, computes
        hierarchical positions for the nodes, and generates a visual representation
        using NetworkX and Matplotlib.

        Workflow
        --------
        1. build_cartesian_tree(arr)
           - Constructs the parent array for the motif using a stack-based approach.
           - Each node points to its nearest smaller value to the left (parent).
           - Time complexity: O(L) where L is the motif length.
           - Output: parent array where parent[i] is the index of the parent of node i.

        2. get_tree_pos(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5)
           - Recursively calculates the (x, y) coordinates for each node.
           - Ensures horizontal spacing between siblings and vertical spacing between levels.
           - Returns a dictionary mapping nodes to positions for drawing.

        3. Build DiGraph (G)
           - Nodes correspond to elements of the motif.
           - Edges represent parent-child relationships derived from the Cartesian Tree.
           - The root node is automatically identified as the node with no parent.

        4. Draw the tree using nx.draw()
           - Nodes are colored light coral (#F08080).
           - Labels display the actual motif values.
           - Edges are gray arrows showing direction from parent to child.
           - The figure is saved to `save_path` as a PNG.

        Example
        -------
        >>> motif = [3, 1, 4, 2]
        >>> plot_ctm_tree(motif, length=4, count=7, save_path='CTM_Tree_L4.png')
        (Generates and saves a tree diagram showing the parent-child structure.)

        Notes
        -----
        - Node size, font size, and arrow size are fixed for clarity.
        - If the tree has no valid root, a spring layout is used as a fallback.
        - Useful for visualizing motif patterns in time series or other sequence data.
        """
    def get_tree_pos(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
        pos = {root: (xcenter, vert_loc)}
        children = list(G.neighbors(root))
        if children:
            dx = width / len(children)
            nextx = xcenter - width / 2 - dx / 2
            for child in children:
                nextx += dx
                pos.update(get_tree_pos(G, child, width=dx, vert_gap=vert_gap,
                                        vert_loc=vert_loc - vert_gap, xcenter=nextx))
        return pos

    def build_cartesian_tree(arr):
        stack = []
        parent = [-1] * len(arr)
        for i in range(len(arr)):
            last = -1
            while stack and arr[stack[-1]] > arr[i]:
                last = stack.pop()
            if stack:
                parent[i] = stack[-1]
            if last != -1:
                parent[last] = i
            stack.append(i)
        return parent

    parent = build_cartesian_tree(motif_tuple)
    G = nx.DiGraph()
    root = -1
    for i, p in enumerate(parent):
        G.add_node(i)
        if p != -1:
            G.add_edge(p, i)
        else:
            root = i

    pos = get_tree_pos(G, root) if root != -1 else nx.spring_layout(G)

    plt.figure(figsize=(10, 6))

    labels = {i: str(motif_tuple[i]) for i in G.nodes()}

    nx.draw(
        G, pos, with_labels=True, labels=labels,
        node_size=1500, node_color='#F08080',
        font_size=11, font_weight='bold',
        arrows=True, arrowsize=20, edge_color='gray'
    )

    motif_str = ', '.join(map(str, motif_tuple))
    title_text = f"CTM Motif Tree\n" \
                 f"Pattern: ({motif_str})\n" \
                 f"Length: {length} | Frequency: {count}"

    plt.title(title_text, fontsize=14, fontweight='bold', pad=30)

    plt.tight_layout()

    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.close()

# ---------------------------------------------
# Load the temperature anomaly dataset
# ---------------------------------------------
def load_data(filename):
    """
        Parses the Berkeley Earth global temperature dataset file.

        Parameters:
        -----------
        filename : str
            The path to the .txt file containing the dataset.

        Returns:
        --------
        tuple
            Returns two tuples:
            1. (monthly_anom, start_year_monthly, start_month_monthly)
            2. (annual_anom, start_year_annual)

        How it works:
        -------------
        The function iterates through the file line by line. It skips comment lines (starting with '%').
        It splits each line to extract the Year, Month, Monthly Anomaly (Column 3), and Annual Anomaly (Column 5).
        It filters out 'NaN' values and stores the valid float data into lists.

        What it does:
        -------------
        Converts the raw text file into usable Python lists for analysis.
        """
    monthly_anom = []
    annual_anom = []
    start_year_monthly = None
    start_month_monthly = None
    start_year_annual = None

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip():
                continue
            parts = line.split()
            try:
                year = int(parts[0])
                month = int(parts[1])

                ma = parts[2]
                if ma != 'NaN':
                    if start_year_monthly is None:
                        start_year_monthly = year
                        start_month_monthly = month
                    monthly_anom.append(float(ma))

                aa = parts[4]
                if aa != 'NaN':
                    if start_year_annual is None:
                        start_year_annual = year
                    annual_anom.append(float(aa))
            except (IndexError, ValueError):
                continue

    return (monthly_anom, start_year_monthly, start_month_monthly), \
           (annual_anom, start_year_annual)


# ---------------------------------------------
# Analyze one motif group to describe its meaning
# ---------------------------------------------
def analyze_motif_group(data, motifs, years, length):
    """
    Calculates statistical properties of a specific motif group.
    """
    if not motifs:
        return None

    prefix_ranks, (count, positions, _) = motifs[0]
    slopes, means, motif_years = [], [], []

    for pos in positions:
        seq = data[pos:pos + length]
        if len(seq) == length:
            slopes.append(seq[-1] - seq[0])
            means.append(np.mean(seq))
            motif_years.append(years[pos])

    avg_slope = np.mean(slopes)
    avg_mean = np.mean(means)
    time_span = f"{int(min(motif_years))}–{int(max(motif_years))}"
    trend_type = "warming trend" if avg_slope > 0 else "cooling trend" if avg_slope < 0 else "stable pattern"

    return {
        "name": prefix_ranks,
        "count": count,
        "slope": avg_slope,
        "mean": avg_mean,
        "time_span": time_span,
        "trend_type": trend_type
    }


# ---------------------------------------------
# Plot Grand Summary (Consolidated Single Graph)
# ---------------------------------------------
def plot_grand_summary(all_opm_data, all_ctm_data, save_path):
    """
    Plots a single Consolidated chart showing Top Motifs for ALL lengths mixed together.
    """

    colors_opm = {
        3: '#87CEFA',
        6: '#4169E1',
        9: '#000080',
        12: '#8A2BE2'
    }

    colors_ctm = {
        3: '#FFD700',
        6: '#FF4500',
        9: '#B22222',
        12: '#4B0082'
    }
    def flatten_and_sort(data_dict, limit=20):
        combined_list = []
        for length, motifs in data_dict.items():
            for m in motifs:
                m['origin_length'] = length
                combined_list.append(m)

        # Sort by Count (Frequency) Descending
        combined_list.sort(key=lambda x: x['count'], reverse=True)
        return combined_list[:limit]


    top_opm_overall = flatten_and_sort(all_opm_data, limit=100)
    top_ctm_overall = flatten_and_sort(all_ctm_data, limit=100)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(26, 12))
    fig.suptitle("Comparison of Top Motifs Across All Lengths", fontsize=24, y=0.95)

    # --- Plot OPM (Left) ---
    if top_opm_overall:
        counts = [m['count'] for m in top_opm_overall]
        labels = [f"{m['tuple']}" for m in top_opm_overall]
        bar_colors = [colors_opm.get(m['origin_length'], 'blue') for m in top_opm_overall]

        bars = ax1.bar(range(len(counts)), counts, color=bar_colors, alpha=0.9, edgecolor='black')
        ax1.set_ylabel("Frequency", fontsize=14)
        ax1.set_title("Top OPM Motifs (All Lengths Mixed)", fontsize=18, fontweight='bold')
        ax1.set_xticks(range(len(counts)))
        ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=10, fontweight='bold')
        ax1.grid(axis='y', linestyle='--', alpha=0.5)

        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom',
                     fontsize=11, fontweight='bold')

        legend_patches = [mpatches.Patch(color=color, label=f'Length {L}') for L, color in colors_opm.items()]
        ax1.legend(handles=legend_patches, title="Motif Length", fontsize=12)

    else:
        ax1.text(0.5, 0.5, "No Motifs Found", ha='center')

    if top_ctm_overall:
        counts = [m['count'] for m in top_ctm_overall]
        labels = [f"{m['tuple']}" for m in top_ctm_overall]
        bar_colors = [colors_ctm.get(m['origin_length'], 'red') for m in top_ctm_overall]

        bars = ax2.bar(range(len(counts)), counts, color=bar_colors, alpha=0.9, edgecolor='black')
        ax2.set_ylabel("Frequency", fontsize=14)
        ax2.set_title("Top CTM Motifs (All Lengths Mixed)", fontsize=18, fontweight='bold')
        ax2.set_xticks(range(len(counts)))
        ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=10, fontweight='bold')
        ax2.grid(axis='y', linestyle='--', alpha=0.5)

        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom',
                     fontsize=11, fontweight='bold')

        legend_patches = [mpatches.Patch(color=color, label=f'Length {L}') for L, color in colors_ctm.items()]
        ax2.legend(handles=legend_patches, title="Motif Length", fontsize=12)

    else:
        ax2.text(0.5, 0.5, "No Motifs Found", ha='center')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(save_path, dpi=300)
    plt.close(fig)

def plot_summary_across_lengths(summary_opm, summary_ctm, save_path):
    """
    Generates a summary plot comparing the top motif across different lengths (L).

    Parameters:
    -----------
    summary_opm : list
        List of tuples containing OPM results for different lengths.
    summary_ctm : list
        List of tuples containing CTM results for different lengths.
    save_path : str
        File path where the image will be saved.

    Returns:
    --------
    None (Saves an image to disk).

    How it works:
    -------------
    It extracts the frequency counts and motif labels from the input summaries.
    It creates a side-by-side bar chart using Matplotlib: one subplot for OPM and one for CTM.
    It annotates the bars with the exact frequency count.

    What it does:
    -------------
    Visualizes how the most frequent pattern changes as we change the time scale (L).
    """

    lengths = [x[0] for x in summary_opm]

    # Prepare data for OPM
    opm_counts = [x[1] for x in summary_opm]
    opm_labels = [f"L={x[0]}\n{x[2]}" for x in summary_opm]  # Label: L=3 \n (1,2,3)

    # Prepare data for CTM
    ctm_counts = [x[1] for x in summary_ctm]
    ctm_labels = [f"L={x[0]}\n{x[2]}" for x in summary_ctm]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9))
    fig.suptitle("Comparison of the Most Frequent Motif for Each Length (L)", fontsize=20)

    # Plot OPM Summary
    bars1 = ax1.bar(range(len(lengths)), opm_counts, color='royalblue', alpha=0.9, edgecolor='black')
    ax1.set_title("Most Frequent OPM Motifs", fontsize=16)
    ax1.set_ylabel("Frequency", fontsize=14)
    ax1.set_xticks(range(len(lengths)))
    ax1.set_xticklabels(opm_labels, fontsize=11, rotation=0)  # Shows motif at bottom
    ax1.grid(axis='y', linestyle='--', alpha=0.5)

    # Add count numbers
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom', fontsize=12,
                 fontweight='bold')

    # Plot CTM Summary
    bars2 = ax2.bar(range(len(lengths)), ctm_counts, color='crimson', alpha=0.9, edgecolor='black')
    ax2.set_title("Most Frequent CTM Motifs", fontsize=16)
    ax2.set_ylabel("Frequency", fontsize=14)
    ax2.set_xticks(range(len(lengths)))
    ax2.set_xticklabels(ctm_labels, fontsize=11, rotation=0)  # Shows motif at bottom
    ax2.grid(axis='y', linestyle='--', alpha=0.5)

    # Add count numbers
    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom', fontsize=12,
                 fontweight='bold')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(save_path, dpi=300)
    plt.close(fig)


def plot_top_motifs_barchart(opm_motifs, ctm_motifs, length, frequency, save_path, data, years):
    """
    Generates a detailed bar chart of the top 10 most frequent motifs.

    Parameters:
    -----------
    opm_motifs : list
        List of OPM motifs found by the algorithm.
    ctm_motifs : list
        List of CTM motifs found by the algorithm.
    length : int
        The motif length (L) being analyzed.
    frequency : str
        'monthly' or 'annual'.
    save_path : str
        Path to save the .png file.
    data : list
        The raw temperature data (used for analysis logic if needed).
    years : array
        The years array corresponding to the data.

    Returns:
    --------
    None (Saves an image to disk).

    How it works:
    -------------
    1. Selects the top 10 motifs from OPM and CTM results.
    2. Extracts the tuple representation (e.g., '(1, 2, 3)') to use as the X-axis label.
    3. Plots two bar charts side-by-side: OPM (Blue) and CTM (Red).
    4. Adds the exact frequency count on top of every bar.

    What it does:
    -------------
    Provides a clear visual comparison of which specific patterns appear most often in the data.
    """
    # Get top 10
    top_opm = opm_motifs[:10]
    top_ctm = ctm_motifs[:10]

    # Helper to generate detailed labels
    def get_labels_and_data(motifs):
        labels = []
        counts = []
        structured_data = []
        for m in motifs:
            info = analyze_motif_group(data, [m], years, length)
            trend_name = info['trend_type'].replace(" trend", "").capitalize()
            motif_string = str(m[0])
            label = f"{trend_name}\n{motif_string}"
            labels.append(label)
            counts.append(m[1][0])

            structured_data.append({
                'trend': trend_name,
                'tuple': motif_string,
                'count': m[1][0]
            })
        return labels, counts, structured_data

    opm_labels, opm_counts, opm_struct = get_labels_and_data(top_opm)
    ctm_labels, ctm_counts, ctm_struct = get_labels_and_data(top_ctm)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))
    fig.suptitle(f"Top 10 Most Frequent Motifs (L={length}, {frequency})", fontsize=20)

    # Plot OPM
    if top_opm:
        bars1 = ax1.bar(range(len(top_opm)), opm_counts, color='royalblue', alpha=0.8, edgecolor='black')
        ax1.set_title("Order Preserving Motifs (OPM)", fontsize=16)
        ax1.set_ylabel("Frequency", fontsize=14)
        ax1.set_xticks(range(len(top_opm)))
        ax1.set_xticklabels(opm_labels, rotation=45, ha='right', fontsize=10, fontweight='bold')
        ax1.grid(axis='y', linestyle='--', alpha=0.5)
        for bar in bars1:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom',
                     fontsize=12, fontweight='bold')
    else:
        ax1.text(0.5, 0.5, "No OPM Motifs found", ha='center')

    # Plot CTM
    if top_ctm:
        bars2 = ax2.bar(range(len(top_ctm)), ctm_counts, color='crimson', alpha=0.8, edgecolor='black')
        ax2.set_title("Cartesian Tree Motifs (CTM)", fontsize=16)
        ax2.set_ylabel("Frequency", fontsize=12)
        ax2.set_xticks(range(len(top_ctm)))
        ax2.set_xticklabels(ctm_labels, rotation=45, ha='right', fontsize=10, fontweight='bold')
        ax2.grid(axis='y', linestyle='--', alpha=0.5)
        for bar in bars2:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom',
                     fontsize=12, fontweight='bold')
    else:
        ax2.text(0.5, 0.5, "No CTM Motifs found", ha='center')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(save_path, dpi=300)
    plt.close(fig)

    return opm_struct, ctm_struct
# ---------------------------------------------
# Analyze one motif group to describe its meaning
# ---------------------------------------------
def analyze_motif_group(data, motifs, years, length):
    """
    Calculates statistical properties of a specific motif group (Trend, Slope, Time Span).

    Parameters:
    -----------
    data : list
        The raw temperature data.
    motifs : list
        The specific motif group to analyze.
    years : array
        The years array matching the data.
    length : int
        The length (L) of the motif.

    Returns:
    --------
    dict
        A dictionary containing: Name, Count, Slope, Mean, Time Span, and Trend Type.

    How it works:
    -------------
    It iterates over all starting 'positions' where this motif occurs.
    For each occurrence, it slices the raw data to get the sequence.
    It calculates the slope (End Value - Start Value) and the Mean.
    It determines the trend (Warming if slope > 0, Cooling if slope < 0).

    What it does:
    -------------
    Converts a raw motif (e.g., '(1,2,3)') into a meaningful description like "Warming Trend".
    """
    if not motifs:
        return None

    prefix_ranks, (count, positions, _) = motifs[0]
    slopes, means, motif_years = [], [], []

    for pos in positions:
        seq = data[pos:pos + length]
        if len(seq) == length:
            slopes.append(seq[-1] - seq[0])
            means.append(np.mean(seq))
            motif_years.append(years[pos])

    avg_slope = np.mean(slopes)
    avg_mean = np.mean(means)
    time_span = f"{int(min(motif_years))}–{int(max(motif_years))}"
    trend_type = "warming trend ↑" if avg_slope > 0 else "cooling trend ↓" if avg_slope < 0 else "stable pattern →"

    return {
        "name": prefix_ranks,
        "count": count,
        "slope": avg_slope,
        "mean": avg_mean,
        "time_span": time_span,
        "trend_type": trend_type
    }


# ---------------------------------------------
# Summarize motifs statistically (for text output)
# ---------------------------------------------
def summarize_motif(data, motifs, length, frequency, start_info):
    """
        Generates a text summary for the top 3 motifs to be written to the results file.

        Parameters:
        -----------
        data : list
            The raw data.
        motifs : list
            The list of found motifs.
        length : int
            Motif length L.
        frequency : str
            'monthly' or 'annual'.
        start_info : tuple/int
            Information about the start year/month.

        Returns:
        --------
        list
            A list of strings, where each string is a sentence describing a motif.

        How it works:
        -------------
        It processes the top 3 motifs.
        It calculates the average temperature change (trend) for each.
        It constructs a formatted string describing the motif rank, count, time span, and physical meaning.

        What it does:
        -------------
        Provides a human-readable summary for the 'analysis_results.txt' file.
        """
    summaries = []
    if not motifs:
        return ["No motifs found."]

    year = start_info[0] if isinstance(start_info, tuple) else start_info
    month = start_info[1] if isinstance(start_info, tuple) else 1
    start_val = year + (month - 1) / 12.0
    years = start_val + (np.arange(len(data)) / (12.0 if frequency == 'monthly' else 1.0))

    for prefix_ranks, (count, positions, windows) in motifs[:3]:  # Analyze top 3 motifs
        means, trends = [], []
        motif_years = []
        for pos in positions:
            seq = data[pos:pos + length]
            means.append(np.mean(seq))
            trends.append(seq[-1] - seq[0])
            motif_years.append(years[pos])

        avg_mean = np.mean(means)
        avg_trend = np.mean(trends)
        time_span = f"{int(min(motif_years))} And {int(max(motif_years))}"
        trend_desc = "warming trend" if avg_trend > 0 else "cooling trend" if avg_trend < 0 else "stable pattern"

        summaries.append(
            f" Motif {prefix_ranks} appeared {count} times ({frequency}), "
            f"between {time_span}. "
            f"It represents an average {trend_desc} of {avg_trend:.3f}C over {length} units."
        )

    return summaries


# ---------------------------------------------
# Plot main time series, motif highlights, and interpretations
# ---------------------------------------------
def plot_motifs(data, motifs_opm, motifs_ctm, length, frequency, save_global, save_zoom, start_info):
    """
        Creates the main visualization: Global Time Series and Average Motif Shape.

        Parameters:
        -----------
        data : list
            Raw temperature data.
        motifs_opm : list
            List of OPM motifs.
        motifs_ctm : list
            List of CTM motifs.
        length : int
            Motif length L.
        frequency : str
            'monthly' or 'annual'.
        save_global : str
            Path to save global plot.
        save_zoom : str
            Path to save zoomed plot.
        start_info : tuple/int
            Start year/month info.

        Returns:
        --------
        None (Saves images to disk).

        How it works:
        -------------
        1. Plot 1 (Top): Plots the full temperature timeline (Gray).
           It overlays the occurrences of the #1 OPM motif (Blue) and #1 CTM motif (Red).
           It adds text boxes with statistical interpretation (Slope, Count).
        2. Plot 2 (Bottom): Calculates the 'Average Shape' of the motif.
           It takes all occurrences, normalizes them to 0-1 range, and computes the mean vector.
           It visualizes this mean shape to show the abstract pattern.
        3. Generates a second image (Zoomed View) focusing on the first occurrence of the OPM motif.

        What it does:
        -------------
        Visualizes where the motifs appear in history and what their average shape looks like.
        """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 9), gridspec_kw={'height_ratios': [2, 1]})
    fig.subplots_adjust(hspace=0.4)

    year = start_info[0] if isinstance(start_info, tuple) else start_info
    month = start_info[1] if isinstance(start_info, tuple) else 1
    start_val = year + (month - 1) / 12.0
    years = start_val + (np.arange(len(data)) / 12.0)

    # ----------------------
    # (1) Time-Series Plot
    # ----------------------
    ax1.plot(years, data, color='gray', alpha=0.4, label='Temperature Anomaly', zorder=1)

    def add_motif_to_plot(motifs, color, style, width, name, z):
        if motifs:
            _, (count, positions, _) = motifs[0]
            for i, pos in enumerate(positions):
                ax1.plot(years[pos:pos + length], data[pos:pos + length],
                         color=color, linestyle=style, linewidth=width,
                         label=f'{name} ({count}x)' if i == 0 else "", zorder=z)

    add_motif_to_plot(motifs_opm, 'royalblue', '-', 2.5, 'OPM Motif', 3)
    add_motif_to_plot(motifs_ctm, 'crimson', '--', 2, 'CTM Motif', 2)

    ax1.set_title(f"Recurring Temperature Patterns (L={length}, {frequency})", fontsize=15)
    ax1.set_xlabel("Year", fontsize=12)
    ax1.set_ylabel("Temperature Anomaly (°C)", fontsize=12)
    ax1.legend(loc="upper left", frameon=True, shadow=True)
    ax1.grid(True, linestyle=':', alpha=0.6)

    # ----------------------
    # (2) Average Normalized Motif Shapes
    # ----------------------
    def plot_avg_shape(motifs, color, label):
        if motifs:
            _, (_, positions, _) = motifs[0]
            motifs_arrays = [data[pos:pos + length] for pos in positions if pos + length <= len(data)]
            motifs_arrays = np.array(motifs_arrays)
            motifs_norm = (motifs_arrays - np.min(motifs_arrays, axis=1, keepdims=True)) / \
                          (np.ptp(motifs_arrays, axis=1, keepdims=True) + 1e-8)

            avg_shape = np.mean(motifs_norm, axis=0)
            ax2.plot(np.arange(length), avg_shape, color=color, linewidth=2.5, label=label)
            ax2.fill_between(np.arange(length),
                             avg_shape - np.std(motifs_norm, axis=0),
                             avg_shape + np.std(motifs_norm, axis=0),
                             color=color, alpha=0.2)

    plot_avg_shape(motifs_opm, 'royalblue', 'OPM Average Shape')
    plot_avg_shape(motifs_ctm, 'crimson', 'CTM Average Shape')

    ax2.set_title("Average Motif Shapes (Normalized 0–1)", fontsize=13)
    ax2.set_xlabel("Relative Time Within Motif", fontsize=12)
    ax2.set_ylabel("Normalized Anomaly", fontsize=12)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(frameon=True, shadow=True)

    # ----------------------
    # (3) Add Interpretation Boxes on the Graph
    # ----------------------
    opm_info = analyze_motif_group(data, motifs_opm, years, length)
    ctm_info = analyze_motif_group(data, motifs_ctm, years, length)

    if opm_info:
        text = (f"OPM Motif:\n"
                f"• Type: {opm_info['trend_type']}\n"
                f"• ΔT ≈ {opm_info['slope']:.3f}°C\n"
                f"• Avg anomaly: {opm_info['mean']:.3f}°C\n"
                f"• {opm_info['count']} occurrences\n"
                f"• Years: {opm_info['time_span']}")
        ax1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=3, frameon=True, shadow=True)
        ax1.text(0.02, 0.95, text, transform=ax1.transAxes, fontsize=10,
                 color='royalblue', va='top',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='royalblue'))
    if ctm_info:
        text = (f"CTM Motif:\n"
                f"• Type: {ctm_info['trend_type']}\n"
                f"• ΔT ≈ {ctm_info['slope']:.3f}°C\n"
                f"• Avg anomaly: {ctm_info['mean']:.3f}°C\n"
                f"• {ctm_info['count']} occurrences\n"
                f"• Years: {ctm_info['time_span']}")
        ax1.text(0.68, 0.95, text, transform=ax1.transAxes, fontsize=10,
                 color='crimson', va='top',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='crimson'))

    plt.savefig(save_global, dpi=300, bbox_inches='tight')

    if motifs_opm:
        first_pos = motifs_opm[0][1][1][0]
        ax1.set_xlim(years[max(0, first_pos - 24)], years[min(len(data) - 1, first_pos + length + 24)])
        ax1.set_title(f"Zoomed View – OPM Motif (Year {int(years[first_pos])})", fontsize=14)
        plt.savefig(save_zoom, dpi=300, bbox_inches='tight')

    plt.close(fig)


# ---------------------------------------------
# Main Execution
# ---------------------------------------------
def main():
    """
        Main execution function for the Climate Motif Analysis.

        Input:
        ------
        None (Reads files from '../dataset/' and calls other functions).

        How it works:
        -------------
        1. Sets up output directories.
        2. Loads the Berkeley Earth Dataset.
        3. Scales float data to integers (x100) for algorithmic efficiency.
        4. Loops through defined motif lengths L = [3, 6, 9, 12].
        5. For each length:
           - Runs 'order_preserving_match' (OPM) to find rank-based patterns.
           - Runs 'discover_ctm_motifs' (CTM) to find tree-based patterns.
           - Calls plotting functions to generate Time Series and Bar Charts.
           - Writes a textual summary of the results to 'analysis_results.txt'.

        What it does:
        -------------
        Orchestrates the entire project pipeline from data loading to result generation.
        """
    os.makedirs('../output', exist_ok=True)
    (m_data, m_year, m_month), (a_data, a_year) = load_data('../dataset/Berkeley Earth global temperature.txt')

    SCALE = 100
    monthly_scaled = [int(round(x * SCALE)) for x in m_data]
    annual_scaled = [int(round(x * SCALE)) for x in a_data]
    grand_summary_opm = {}
    grand_summary_ctm = {}
    K = [2, 3, 5, 7, 10]  # Minimum motif frequencies
    lengths = [3,6,9,12]     # Motif lengths
    top_n = 10

    with open('../output/analysis_results.txt', 'w') as f:
        for length, k in zip(lengths, K):
            if length == 12:
                data_to_match = annual_scaled
                original_data = a_data
                frequency = 'annual'
                time_info = a_year
            else:
                data_to_match = monthly_scaled
                original_data = m_data
                frequency = 'monthly'
                time_info = (m_year, m_month)
            y_start = time_info[0] if isinstance(time_info, tuple) else time_info
            m_start = time_info[1] if isinstance(time_info, tuple) else 1
            step = 1.0 if frequency == 'annual' else (1.0 / 12.0)
            start_val = y_start + (m_start - 1) / 12.0
            years = start_val + (np.arange(len(original_data)) * step)
            start = time.time()
            opm_recurring_motifs = sorted(order_preserving_match(data_to_match, length, k),
                                          key=lambda x: x[1][0], reverse=True)
            ctm_recurring_motifs = sorted(discover_ctm_motifs(data_to_match, length, k),
                                          key=lambda x: x[1][0], reverse=True)
            end = time.time()
            run_time = end - start

            save_shapes_opm = f"../output/Top10_Shapes_OPM_L{length}.png"
            save_shapes_ctm = f"../output/Top10_Shapes_CTM_L{length}.png"
            #Generate Top10_shapes_OPM
            plot_top_10_motifs_shapes(original_data, opm_recurring_motifs, length,
                                      frequency, save_shapes_opm, "OPM")
            #Generate Top10_shapes_CTM
            plot_top_10_motifs_shapes(original_data, ctm_recurring_motifs, length,
                                      frequency, save_shapes_ctm, "CTM")
            # --- Draw CTM motif trees ---
            os.makedirs('../output/CTM_Trees', exist_ok=True)
            for i, (prefix_ranks, (count, positions, windows)) in enumerate(ctm_recurring_motifs[:5]):  # top 5
                motif_tuple = prefix_ranks
                tree_path = f"../output/CTM_Trees/CTM_Tree_L{length}_#{i + 1}.png"
                plot_ctm_tree(motif_tuple, length, count, tree_path)

            f.write(f"\n=== Analysis for Motif Length {length} ({frequency}) ===\n")
            f.write(f"Minimum occurrences: {k}\n")
            f.write(f"Top {top_n} recurring motifs:\n")

            if not opm_recurring_motifs and not ctm_recurring_motifs:
                f.write("No motifs found.\n")
                continue

            save_global = f"../output/Global_{frequency}_L{length}.png"
            save_zoom = f"../output/Zoom_{frequency}_L{length}.png"
            plot_motifs(original_data, opm_recurring_motifs, ctm_recurring_motifs,
                        length, frequency, save_global, save_zoom, time_info)
            save_bar = f"../output/BarChart_{frequency}_L{length}.png"

            opm_data, ctm_data = plot_top_motifs_barchart(
                opm_recurring_motifs, ctm_recurring_motifs,
                length, frequency, save_bar,
                original_data, years)
            grand_summary_opm[length] = opm_data
            grand_summary_ctm[length] = ctm_data
            # Write detailed motif data
            f.write("\n--- OPM Results ---\n")
            for prefix_ranks, (count, positions, windows) in opm_recurring_motifs[:top_n]:
                f.write("\n" + "-" * 40)
                f.write(f"\nMotif order representation: {prefix_ranks}")
                f.write(f"\nFrequency: {count}")
                f.write(f"\nPositions: {positions}")
                f.write(f"\nSubstrings: {windows}")

            f.write("\n\n--- CTM Results ---\n")
            for prefix_ranks, (count, positions, windows) in ctm_recurring_motifs[:top_n]:
                f.write("\n" + "-" * 40)
                f.write(f"\nMotif order representation: {prefix_ranks}")
                f.write(f"\nFrequency: {count}")
                f.write(f"\nPositions: {positions}")
                f.write(f"\nSubstrings: {windows}")

            # Add interpretation summary
            f.write("\n\n--- Interpretation Summary ---\n")
            opm_summary = summarize_motif(original_data, opm_recurring_motifs, length, frequency, time_info)
            ctm_summary = summarize_motif(original_data, ctm_recurring_motifs, length, frequency, time_info)
        print("Generating Grand Summary Plot...")
        plot_grand_summary(grand_summary_opm, grand_summary_ctm, "../output/Grand_Summary_All_Motifs.png")
    print(f"Run time: {run_time}")
    print("Analysis complete.")
    print("Results saved in '../output/analysis_results.txt'")
    print("Plots generated in '../output/' folder.")


if __name__ == '__main__':
    main()

