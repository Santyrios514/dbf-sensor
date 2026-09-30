"""Optimización en malla de la geometría del sensor remolcado (spec v3).

Módulo nuevo sobre `sensor_lastre`: importa y envuelve sus funciones (volúmenes, masas, llenado
de lastre, CG/SM/trim, Barrowman, puente con OpenRocket) sin modificarlas. Ver
docs/mapa_reuso.md. Internamente todo está en SI; la E/S usa mm, g y grados.
"""

__version__ = "0.1.0"
