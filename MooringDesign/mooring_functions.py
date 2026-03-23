import numpy as np
from matplotlib.lines import Line2D


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


def report_load_line_status(system, max_utilization=0.5):
    """Report only seabed-contact length (grounded line)."""
    _ = max_utilization  # kept for backward compatibility with existing calls

    print("\nSeabed-contact length:")
    print("Line | L_total(m) L_on_seabed(m) Seabed(%)")

    for i, line in enumerate(system.lineList):
        L_total = float(line.L)
        L_on_seabed = float(getattr(line, "LBot", 0.0))
        seabed_pct = 100.0 * L_on_seabed / L_total if L_total > 0 else np.nan

        print(f"{i+1:>4} | {L_total:10.3f} {L_on_seabed:13.3f} {seabed_pct:8.2f}")


def annotate_mooring_plot(system, ax, show_small_body_axes=True, body_axes_len=0.25, show_line_numbers=True):
    """Add body axes triad and line-number labels to an existing MoorPy plot axis."""
    if show_small_body_axes:
        x0, y0, z0, *_ = system.bodyList[0].r6
        ax.plot([x0, x0 + body_axes_len], [y0, y0], [z0, z0], color='r')
        ax.plot([x0, x0], [y0, y0 + body_axes_len], [z0, z0], color='g')
        ax.plot([x0, x0], [y0, y0], [z0, z0 + body_axes_len], color='b')

    if show_line_numbers:
        for i, line in enumerate(system.lineList):
            r_mid = 0.5 * (line.rA + line.rB)
            ax.text(r_mid[0], r_mid[1], r_mid[2], f"{i+1}", color='k', fontsize=9)


def plot_body_axes(system, ax, body_axes_len=0.25, color='r', enabled=True):
    """Draw a small body-fixed axes triad at the body position using a single color."""
    if not enabled:
        return

    x0, y0, z0, *_ = system.bodyList[0].r6
    ax.plot([x0, x0 + body_axes_len], [y0, y0], [z0, z0], color=color)
    ax.plot([x0, x0], [y0, y0 + body_axes_len], [z0, z0], color=color)
    ax.plot([x0, x0], [y0, y0], [z0, z0 + body_axes_len], color=color)


def style_chain_rope_plot(
    ax,
    n_lines,
    chain_color='black',
    rope_color='teal',
    title='Initial Mooring Configuration (no thrust)',
    legend_loc='upper right',
    legend_anchor=(1.10, 0.98),
):
    """Color chain/rope segments and attach a matching legend/title to an existing axis."""
    for i, artist in enumerate(ax.lines[-n_lines:]):
        artist.set_color(chain_color if i % 2 == 0 else rope_color)

    ax.set_title(title)
    legend_handles = [
        Line2D([0], [0], color=chain_color, lw=2, label='Chain'),
        Line2D([0], [0], color=rope_color, lw=2, label='Rope'),
    ]
    ax.legend(handles=legend_handles, loc=legend_loc, bbox_to_anchor=legend_anchor, frameon=True)


def style_thrust_comparison_plot(
    ax,
    title='Deflection under Turbine Thrust',
    baseline_color='black',
    thrust_color='red',
    legend_loc='upper right',
):
    """Set title and legend for no-thrust vs thrust equilibrium comparison."""
    ax.set_title(title)
    legend_handles = [
        Line2D([0], [0], color=baseline_color, lw=2, label='No thrust (baseline)'),
        Line2D([0], [0], color=thrust_color, lw=2, label='With turbine thrust'),
    ]
    ax.legend(handles=legend_handles, loc=legend_loc, frameon=True)


def print_body_response(r6, label='Body response'):
    """Print body displacement and rotation from a 6-DOF state vector [x,y,z,rx,ry,rz]."""
    x, y, z, rx, ry, rz = r6
    print(f"{label}: surge={x:.4f} m, sway={y:.4f} m, heave={z:.4f} m")
    print(
        f"{label}: roll={np.degrees(rx):.3f} deg, "
        f"pitch={np.degrees(ry):.3f} deg, yaw={np.degrees(rz):.3f} deg"
    )


def check_anchor_stability(system, anchor_weight_N, angles):
    """Check if anchors can withstand line loads.

    For each mooring leg, reports tension components and safety margin.
    Vertical component of tension (from line angle) tries to pull anchor up.
    Anchor weight resists uplift.
    """
    print("\nAnchor Stability Report")
    print("Leg | Chain_TA(N) | Angle(°) | T_vertical(N) | T_horizontal(N) | Safety_Margin(N) | Status")
    
    all_safe = True
    for i in range(len(angles)):
        chain_line = system.lineList[2 * i]  # chain is the first line of each leg
        
        # Get anchor and connector positions
        r_anchor = chain_line.rA  # anchor point at seabed
        r_connector = chain_line.rB  # connector point (transition to rope)
        
        # Calculate vertical and horizontal distances
        dz = r_connector[2] - r_anchor[2]  # vertical rise
        dx_dy = np.sqrt((r_connector[0] - r_anchor[0])**2 + (r_connector[1] - r_anchor[1])**2)  # horizontal distance
        
        # Line angle from horizontal
        theta_rad = np.arctan2(dz, dx_dy)
        theta_deg = np.degrees(theta_rad)
        
        # Decompose tension at anchor
        T_total = chain_line.TA
        T_vertical = T_total * np.sin(theta_rad)
        T_horizontal = T_total * np.cos(theta_rad)
        
        # Safety check: anchor weight must exceed vertical pull
        safety_margin = anchor_weight_N - T_vertical
        is_safe = safety_margin > 0
        all_safe = all_safe and is_safe
        
        status = "OK" if is_safe else "FAIL"
        print(f"{i+1:>3} | {T_total:10.2f} | {theta_deg:7.2f} deg | {T_vertical:13.2f} | {T_horizontal:15.2f} | {safety_margin:16.2f} | {status}")
    
    if not all_safe:
        print("\nWARNING: Anchor uplift risk! Increase anchor weight or reduce applied load.")
    else:
        print("\nAll anchors stable against uplift.")
    
    return all_safe
