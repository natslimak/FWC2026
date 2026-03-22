from moorpy.helpers import getLineProps

# Central place to define chain/rope line types used by mooring_design1.py.
# Set mode='catalog' to use MoorPy defaults, or mode='custom' for exact values.

CHAIN_SPEC = {
    "name": "chain120",
    "mode": "custom",  # 'catalog' or 'custom'

    "catalog": {
        "size_mm": 20,
        "material": "chain",
        "source": "default",
    },
    "custom": {
        "d_vol": 0.018,         # [m] volume-equivalent diameter
        "massden": 1.99,        # [kg/m] linear mass density
        "EA": 10100,            # [N] axial stiffness
        "d_nom": 0.01,          # [m] nominal diameter (optional)
        "MBL": 4.6649822e5,     # [N] minimum breaking load (optional)
        "Cd": 1.333,
        "Ca": 1.0,
        "CdAx": 0.639,
        "CaAx": 0.5,
        "cost": 20.68,
        "notes": "custom chain",
        "material": "chain",
    },
}

ROPE_SPEC = {
    "name": "polyester180",
    "mode": "catalog",  # 'catalog' or 'custom'
    
    "catalog": {
        "size_mm": 60,
        "material": "polyester",
        "source": "default",
    },
    "custom": {
        "d_vol": 0.0828,
        "massden": 3.228,
        "EA": 1.2315186e7,
        "d_nom": 0.06,
        "MBL": 6.549729e5,
        "Cd": 1.2,
        "Ca": 1.0,
        "CdAx": 0.2,
        "CaAx": 0.0,
        "cost": 34.33,
        "notes": "custom rope",
        "material": "polyester",
    },
}


def _add_line_type_from_spec(system, spec):
    name = spec["name"]
    mode = spec["mode"].lower()

    if mode == "catalog":
        cat = spec["catalog"]
        system.lineTypes[name] = getLineProps(
            cat["size_mm"],
            cat["material"],
            source=cat.get("source", "default"),
            name=name,
        )
    elif mode == "custom":
        custom = spec["custom"]
        required = ["d_vol", "massden", "EA"]
        missing = [k for k in required if k not in custom]
        if missing:
            raise ValueError(f"Missing custom line fields for '{name}': {missing}")

        system.addLineType(name, custom["d_vol"], custom["massden"], custom["EA"], name=name)
        line_type = system.lineTypes[name]
        line_type["d_nom"] = custom.get("d_nom", custom["d_vol"])
        line_type["MBL"] = custom.get("MBL", 0.0)
        line_type["Cd"] = custom.get("Cd", 0.0)
        line_type["Ca"] = custom.get("Ca", 0.0)
        line_type["CdAx"] = custom.get("CdAx", 0.0)
        line_type["CaAx"] = custom.get("CaAx", 0.0)
        line_type["cost"] = custom.get("cost", 0.0)
        line_type["notes"] = custom.get("notes", "custom line type")
        line_type["material"] = custom.get("material", name)
    else:
        raise ValueError(f"Invalid mode '{mode}' for line type '{name}'. Use 'catalog' or 'custom'.")

    return name


def register_chain_and_rope(system):
    chain_type = _add_line_type_from_spec(system, CHAIN_SPEC)
    rope_type = _add_line_type_from_spec(system, ROPE_SPEC)
    return chain_type, rope_type
