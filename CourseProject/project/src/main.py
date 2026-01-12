import time, os
from OPM import order_preserving_match
from CTM import discover_ctm_motifs
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------
# Load the temperature anomaly dataset
# ---------------------------------------------
def load_data(filename):
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
    """Analyze trend, mean, and time span for a given motif group."""
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
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 9), gridspec_kw={'height_ratios': [2, 1]})
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
    os.makedirs('../output', exist_ok=True)
    (m_data, m_year, m_month), (a_data, a_year) = load_data('../dataset/Berkeley Earth global temperature.txt')

    SCALE = 100
    monthly_scaled = [int(round(x * SCALE)) for x in m_data]
    annual_scaled = [int(round(x * SCALE)) for x in a_data]

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

            start = time.time()
            opm_recurring_motifs = sorted(order_preserving_match(data_to_match, length, k),
                                          key=lambda x: x[1][0], reverse=True)
            ctm_recurring_motifs = sorted(discover_ctm_motifs(data_to_match, length, k),
                                          key=lambda x: x[1][0], reverse=True)
            end = time.time()
            run_time = end - start

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

            f.write("\n[OPM Summary]\n" + "\n".join(opm_summary))
            f.write("\n\n[CTM Summary]\n" + "\n".join(ctm_summary))
            f.write(f"\n\nExecution time: {run_time:.2f} seconds\n")

    print("Analysis complete.")
    print("Results saved in '../output/analysis_results.txt'")
    print("Plots generated in '../output/' folder.")


if __name__ == '__main__':
    main()
