import numpy as np
import moorpy as mp
from moorpy.helpers import getLineProps

# ----- choose some system geometry parameters -----

depth     = 600                             # water depth [m]
angles    = np.radians([60, 180, 300])      # line headings list [rad]
rAnchor   = 1600                            # anchor radius/spacing [m]
zFair     = -21                             # fairlead z elevation [m]
rFair     = 20                              # fairlead radius [m]
lineLength= 1800                            # line unstretched length [m]
typeName  = "chain"                         # identifier string for the line type
# typeName options: "chain", "chain_studlink", "hmpe", "nylon", "polyester", "wire"

# ----- set up the mooring system and floating body -----

# Create new MoorPy System and set its depth
ms = mp.System(depth=depth)

# add a line type
ms.lineTypes[typeName] = getLineProps(120, "chain", source="default", name=typeName)  # this would be 120 mm chain

# Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
ms.addBody(0, np.zeros(6), m=1e6, v=1e3, rM=100, AWP=1e3)

# For each line heading, set the anchor point, the fairlead point, and the line itself
for i, angle in enumerate(angles):

    # create end Points for the line
    ms.addPoint(1, [rAnchor*np.cos(angle), rAnchor*np.sin(angle), -depth])   # create anchor point (type 0, fixed)
    ms.addPoint(1, [  rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])   # create fairlead point (type 0, fixed)

    # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
    ms.bodyList[0].attachPoint(2*i+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair])

    # add a Line going between the anchor and fairlead Points
    ms.addLine(lineLength, typeName, pointA=2*i+1, pointB=2*i+2)


ms.initialize()                                             # make sure everything's connected

ms.solveEquilibrium()                                       # equilibrate
fig, ax = ms.plot()                                         # plot the system in original configuration
ms.unload("sample.txt")                                     # export to MD input file

ms.bodyList[0].f6Ext = np.array([3e6, 0, 0, 0, 0, 0])       # apply an external force on the body
ms.solveEquilibrium()                                       # equilibrate
fig, ax = ms.plot(ax=ax, color='red')                       # plot the system in displaced configuration (on the same plot, in red)



""" Getting the line properties database:"""

import moorpy
import moorpy.helpers as mh
from pathlib import Path
import inspect

print("MoorPy version:", getattr(moorpy, "__version__", "unknown"))

# 1) Preferred: list types from the default DB
if hasattr(mh, "loadLineProps"):
    db = mh.loadLineProps(source="default")
    print("\nTypes/materials in source='default':")
    for k in sorted(db.keys()):
        print(" -", k)
else:
    print("\nloadLineProps() not found in this MoorPy version.")

# 2) Show where DB files live in your installed package
pkg_dir = Path(moorpy.__file__).resolve().parent
print("\nPossible line-property files in package:")
for p in sorted(pkg_dir.rglob("*")):
    if p.suffix.lower() in {".txt", ".csv", ".yaml", ".yml", ".json"}:
        name = p.name.lower()
        if "line" in name or "prop" in name or "rope" in name or "chain" in name:
            print(" -", p)

# 3) If needed, inspect getLineProps implementation directly
print("\ngetLineProps signature:", inspect.signature(mh.getLineProps))