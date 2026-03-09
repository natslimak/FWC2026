import numpy as np
import moorpy as mp
from moorpy.helpers import getLineProps

# ----- System Geometry Parameters -----
depth     = 2                               # water depth [m]
angles    = np.radians([60, 180, 300])      # line headings list [rad]
rAnchor   = 3.0                            # anchor radius/spacing [m]
zFair     = 0                            # fairlead z elevation [m] (below SWL)
rFair     = 0.1                            # fairlead radius [m]

# Mixed line definition (anchor side + fairlead side)
chainLength = 1.5                           # anchor-side chain length [m]
ropeLength  = 2.5                           # fairlead-side rope length [m]
chainType   = "chain120"                    # identifier string for chain segment
ropeType    = "polyester180"                # identifier string for rope segment

# Plot options for small-scale models
show_small_body_axes = True
body_axes_len = 0.25                        # [m] custom triad length if enabled



# ----- Mooring System and Floating Body -----

# Create new MoorPy System and set its depth
ms = mp.System(depth=depth)

# Check that the line geometry is reasonable (not too slack or too tight)
L_total = chainLength + ropeLength
d_straight = np.sqrt((rAnchor - rFair)**2 + (depth + zFair)**2)
slack_pct = 100.0*(L_total - d_straight)/d_straight
print(f"Line geometry check: straight={d_straight:.2f} m, length={L_total:.2f} m, slack={slack_pct:.1f}%")
if slack_pct < 5:
    print("Warning: line is very tight. Consider increasing line lengths or reducing anchor radius.")
elif slack_pct > 50:
    print("Warning: line is very slack. Consider reducing line lengths or increasing anchor radius.")
else:
    print("Line geometry looks reasonable.")

# Add line types used in the mixed line
ms.lineTypes[chainType] = getLineProps(20, "chain", source="default", name=chainType)
ms.lineTypes[ropeType] = getLineProps(60, "polyester", source="default", name=ropeType)

# Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
ms.addBody(0, np.zeros(6), m=1e6, v=1e3, rM=100, AWP=1e3)

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

ms.unload("mooring_design.txt")                                     # export to MD input file

ms.bodyList[0].f6Ext = np.array([0, 0, 0, 0, 0, 0])       # apply an external force on the body
#                              Fx, Fy, Fz, Mx, My, Mz
ms.solveEquilibrium()                                       # equilibrate
fig, ax = ms.plot(ax=ax, color='red', draw_body=False)      # plot displaced config without default body axes


# ----- Linearized restoring (stiffness) matrix about current equilibrium -----
# DOF order is typically: [surge, sway, heave, roll, pitch, yaw]
# For the whole structure including the lines 
#C = ms.getSystemStiffness(DOFtype='free', dx=0.01, dth=0.01)

# Compute mooring stiffness matrix (floater)
Cmoor = ms.getCoupledStiffness()

print("\nMooring Restoring Matrix (6x6):")
print(Cmoor)

# Print line tensions
print("\nLine tensions:")
for i, line in enumerate(ms.lineList):
    print(f"Line {i+1}: {line.TA:.2f} N")