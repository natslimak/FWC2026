import numpy as np
import moorpy as mp
from line_type_specs import register_chain_and_rope
from mooring_functions import (
    print_line_tensions_by_leg,
    report_line_geometry,
    report_post_load_line_status,
)

# ----- System Geometry Parameters -----
depth     = 2                               # water depth [m]
angles    = np.radians([60, 180, 300])      # line headings list [rad]
rAnchor   = 3.0                             # anchor radius/spacing [m]
zFair     = 0                               # fairlead z elevation [m] (below SWL)
rFair     = 0.1                             # fairlead radius [m]

# Mixed line definition (anchor side + fairlead side)
chainLength = 1.5                           # anchor-side chain length [m]
ropeLength  = 2.5                           # fairlead-side rope length [m]

# Plot options for small-scale models
show_line_numbers = True
show_small_body_axes = True
body_axes_len = 0.25                        # [m] custom triad length if enabled

# Floater Geometry 
m_floater = 80.0                        # kg
rho_water = 1025.0                      # kg/m^3 (seawater)
v_disp = m_floater / rho_water          # ~0.078 m^3 to float 80 kg
a = 2.0                                 # m, side length of an equilateral triangle (for estimating waterplane area)
AWP_est = (np.sqrt(3) / 4.0) * a**2     # m^2 (if full triangle is at waterplane)
rM_est = np.array([0.0, 0.0, 0.15])     # m, simple initial guess of mass center above the waterline (for estimating hydrostatic stiffness)


# ----- Mooring System and Floating Body -----

# Create new MoorPy System and set its depth
ms = mp.System(depth=depth)

# Check that the line geometry is reasonable (not too slack or too tight)
L_total, d_straight, slack_pct = report_line_geometry(
    chainLength, ropeLength, rAnchor, rFair, depth, zFair)

# Add line types from separate specification file
chainType, ropeType = register_chain_and_rope(ms)  # <- add to line_type_specs.py detials from manufacturer
print(f"Using line types: chain={chainType}, rope={ropeType}")

# Add a floating body 
# Type -1: Use when you want to get the Cmoor matrix for the lines without a floating body (e.g. for a fixed-bottom test case or to isolate line stiffness). The body will be massless and volume-less, so it won't affect the equilibrium or stiffness results, but it allows you to use the line geometry and tensioning features of MoorPy without needing a full floating body definition.
# Type 0: Use when you want to full C matrix including the hydrostatic stiffness contributions from the displaced volume and waterplane area. The body will be initialized with a default mass and volume, but these can be adjusted as needed. This is a good option for small-scale models where you want to capture the basic hydrostatic effects without needing a detailed floating body definition.
ms.addBody(
    mytype=-1, r6=np.zeros(6),
    m=m_floater,
    v=v_disp,
    AWP=AWP_est,
    rM=rM_est
)

# For each line heading, set the anchor point, the fairlead point, and the line itself
for i, angle in enumerate(angles):

    # create end points for the line
    ms.addPoint(1, [rAnchor*np.cos(angle), rAnchor*np.sin(angle), -depth])   # create anchor point (type 0, fixed)
    ms.addPoint(1, [  rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])   # create fairlead point (type 0, fixed)

    # Create a free connector point where chain transitions to rope
    x_anchor = rAnchor*np.cos(angle)
    y_anchor = rAnchor*np.sin(angle)
    x_fair = rFair*np.cos(angle)
    y_fair = rFair*np.sin(angle)
    frac_chain = chainLength/(chainLength + ropeLength)
    x_mid = x_anchor + frac_chain*(x_fair - x_anchor)
    y_mid = y_anchor + frac_chain*(y_fair - y_anchor)
    z_mid = -depth + frac_chain*(zFair + depth)
    ms.addPoint(0, [x_mid, y_mid, z_mid])

    # Attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
    ms.bodyList[0].attachPoint(3*i+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair])

    # Add two line segments: chain (anchor->connector) then rope (connector->fairlead)
    anchor_id = 3*i + 1
    fairlead_id = 3*i + 2
    connector_id = 3*i + 3

    ms.addLine(chainLength, chainType, pointA=anchor_id, pointB=connector_id)
    ms.lineList[-1].nNodes = 30
    ms.addLine(ropeLength, ropeType, pointA=connector_id, pointB=fairlead_id)
    ms.lineList[-1].nNodes = 30


# ----- Run the Model -----
ms.initialize()                                             # make sure everything's connected

ms.solveEquilibrium()                                       # equilibrate
fig, ax = ms.plot(draw_body=False)                          # hide default 5 m body axes

if show_small_body_axes:
    x0, y0, z0, *_ = ms.bodyList[0].r6
    ax.plot([x0, x0 + body_axes_len], [y0, y0], [z0, z0], color='r')
    ax.plot([x0, x0], [y0, y0 + body_axes_len], [z0, z0], color='g')
    ax.plot([x0, x0], [y0, y0], [z0, z0 + body_axes_len], color='b')

if show_line_numbers:
    for i, line in enumerate(ms.lineList):
        r_mid = 0.5 * (line.rA + line.rB)
        ax.text(r_mid[0], r_mid[1], r_mid[2], f"{i+1}", color='k', fontsize=9)

ms.unload("mooring_design.txt")                                     # export to MD input file

ms.bodyList[0].f6Ext = np.array([0, 0, 0, 0, 0, 0])       # apply an external force on the body
#                              Fx, Fy, Fz, Mx, My, Mz
ms.solveEquilibrium(tol=1e-6)                                       # equilibrate
#fig, ax = ms.plot(ax=ax, color='red', draw_body=False)      # plot displaced config without default body axes


# ----- Linearized restoring (stiffness) matrix about current equilibrium -----
# DOF order is typically: [surge, sway, heave, roll, pitch, yaw]
# For the whole structure including the lines 
# C = ms.getSystemStiffness(DOFtype='free', dx=0.01, dth=0.01)

# Compute mooring stiffness matrix (floater)
Cmoor = ms.getCoupledStiffness(lines_only=True)
Ctotal = ms.getCoupledStiffness(lines_only=False)
Chydro = Ctotal - Cmoor

print("\nC_Moor Matrix (6x6):")
print(Cmoor)

# Print line tensions grouped by mooring leg (chain + rope)
print_line_tensions_by_leg(ms, angles)

# Baseline (pre-thrust) slack and holding margin.
report_post_load_line_status(ms, max_utilization=0.5)

# Other important outputs to check:
# Total chain weight per leg (dry and submerged).
chain_props = ms.lineTypes[chainType]
chain_weight_dry = chainLength * chain_props.get("m", chain_props.get("massden")) * 9.81
chain_weight_wet = chainLength * chain_props.get("w", np.nan)  # 'w' is submerged weight per length [N/m]
print(f"\nChain dry weight per leg: {chain_weight_dry:.2f} N")
print(f"Chain submerged weight per leg: {chain_weight_wet:.2f} N")



# ----- Test of the system when turbines are operating -----

# Switch to free-body mode for physical deflection under applied thrust.
ms.bodyList[0].type = 0
ms.bodyList[0].f6Ext = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
ms.solveEquilibrium(tol=1e-6)

# Comparison plot: black = no thrust, red = with thrust.
fig_cmp, ax_cmp = ms.plot(draw_body=False, color='black')
if show_small_body_axes:
    x0, y0, z0, *_ = ms.bodyList[0].r6
    ax_cmp.plot([x0, x0 + body_axes_len], [y0, y0], [z0, z0], color='k', alpha=0.6)
    ax_cmp.plot([x0, x0], [y0, y0 + body_axes_len], [z0, z0], color='k', alpha=0.6)
    ax_cmp.plot([x0, x0], [y0, y0], [z0, z0 + body_axes_len], color='k', alpha=0.6)

# Turbine thrust test (2 turbines, 10 kgf each)
n_turbines = 2
thrust_kgf_each = 5 # <- adjust this value with the actual thrust per turbine in kgf (1 kgf = 9.81 N)
thrust_total_N = n_turbines * thrust_kgf_each * 9.81


# Apply thrust in +x direction. Format: [Fx, Fy, Fz, Mx, My, Mz]
ms.bodyList[0].f6Ext = np.array([thrust_total_N, 0.0, 0.0, 0.0, 0.0, 0.0])
ms.solveEquilibrium(tol=1e-6)
fig_cmp, ax_cmp = ms.plot(ax=ax_cmp, color='red', draw_body=False)
if show_small_body_axes:
    x1, y1, z1, *_ = ms.bodyList[0].r6
    ax_cmp.plot([x1, x1 + body_axes_len], [y1, y1], [z1, z1], color='r')
    ax_cmp.plot([x1, x1], [y1, y1 + body_axes_len], [z1, z1], color='r')
    ax_cmp.plot([x1, x1], [y1, y1], [z1, z1 + body_axes_len], color='r')

x, y, z, rx, ry, rz = ms.bodyList[0].r6
print("\n=== RESULTS WHEN TURBINE IS OPERATING ===")
print(f"\nApplied turbine thrust: {n_turbines} x {thrust_kgf_each:.1f} kgf = {thrust_total_N:.2f} N")
print(f"Body response: surge={x:.4f} m, sway={y:.4f} m, heave={z:.4f} m")
print(
    f"Body response: roll={np.degrees(rx):.3f} deg, "
    f"pitch={np.degrees(ry):.3f} deg, yaw={np.degrees(rz):.3f} deg"
)

if ms.bodyList[0].type != 0:
    print("Note: body type is not free (type != 0). To see physical deflection, set body type to 0.")


# Print line tensions grouped by mooring leg (chain + rope)
print_line_tensions_by_leg(ms, angles)

# Post-thrust slack and holding margin.
report_post_load_line_status(ms, max_utilization=0.5)

