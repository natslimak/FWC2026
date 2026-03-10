import numpy as np


def report_line_geometry(chain_length, rope_length, r_anchor, r_fair, depth, z_fair):
    """Print a quick straight-length vs total-length slack check and return values."""
    L_total = chain_length + rope_length
    d_straight = np.sqrt((r_anchor - r_fair) ** 2 + (depth + z_fair) ** 2)
    slack_pct = 100.0 * (L_total - d_straight) / d_straight

    print(f"Line geometry check: straight={d_straight:.2f} m, length={L_total:.2f} m, slack={slack_pct:.1f}%")
    if slack_pct < 5:
        print("Warning: line is very tight. Consider increasing line lengths or reducing anchor radius.")
    elif slack_pct > 50:
        print("Warning: line is very slack. Consider reducing line lengths or increasing anchor radius.")
    else:
        print("Line geometry looks reasonable.")

    return L_total, d_straight, slack_pct


def print_line_tensions_by_leg(system, angles, title="\nLine tensions by mooring leg:"):
    """Print chain and rope end tensions grouped by mooring leg."""
    print(title)
    for i in range(len(angles)):
        chain_line = system.lineList[2 * i]
        rope_line = system.lineList[2 * i + 1]
        print(
            f"Leg {i+1}: "
            f"chain(L{2*i+1}) TA={chain_line.TA:.2f} N, TB={chain_line.TB:.2f} N | "
            f"rope(L{2*i+2}) TA={rope_line.TA:.2f} N, TB={rope_line.TB:.2f} N"
        )


def report_post_load_line_status(system, max_utilization=0.5):
    """Report only seabed-contact length (grounded line)."""
    _ = max_utilization  # kept for backward compatibility with existing calls

    print("\nPost-load seabed-contact length:")
    print("Line | L_total(m) L_on_seabed(m) Seabed(%)")

    for i, line in enumerate(system.lineList):
        L_total = float(line.L)
        L_on_seabed = float(getattr(line, "LBot", 0.0))
        seabed_pct = 100.0 * L_on_seabed / L_total if L_total > 0 else np.nan

        print(f"{i+1:>4} | {L_total:10.3f} {L_on_seabed:13.3f} {seabed_pct:8.2f}")
