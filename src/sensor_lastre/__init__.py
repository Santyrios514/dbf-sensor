"""Volumen utilizable y lastre del sensor remolcado (DBF 2026-27, UPB).

Internamente todo está en SI (m, kg, s, rad). La E/S (YAML, CSV) usa mm, g, cm³ y grados.
"""

__version__ = "0.1.0"

G0 = 9.80665  # m/s^2, gravedad estándar (única constante física permitida fuera de ISA)
