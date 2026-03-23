#!/usr/bin/env python3
"""Fix Unicode encoding issues in mooring_design1.py"""

# Read the file
with open('mooring_design1.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Replace problematic characters
replacements = [
    ('N²·s/rad', 'N^2 s/rad'),
    ('N² s/rad', 'N^2 s/rad'),
    ('Angle(░)', 'Angle (deg)'),
]

for old, new in replacements:
    content = content.replace(old, new)

# Write back
with open('mooring_design1.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("Encoding issues fixed!")
