"""Lectura del `.ork` (ZIP con un `rocket.ork` en XML) sin JVM, para verificar la configuración.

Solo lee lo que usa este proyecto: NoseCone, BodyTube, Transition y FreeformFinSet.
Valores en SI, tal como los guarda OpenRocket.
"""

from __future__ import annotations

import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

ETIQUETAS = {"nosecone": "NoseCone", "bodytube": "BodyTube", "transition": "Transition",
             "freeformfinset": "FreeformFinSet"}


@dataclass
class ComponenteXML:
    tipo: str
    nombre: str
    campos: dict[str, str] = field(default_factory=dict)
    material: str = ""
    densidad: float = float("nan")
    puntos: list[tuple[float, float]] = field(default_factory=list)
    padre: str | None = None

    def num(self, clave: str) -> float:
        """Valor numérico de un campo; admite 'auto 0.0325' (radio automático)."""
        v = self.campos[clave].split()
        return float(v[-1])

    def es_auto(self, clave: str) -> bool:
        return self.campos.get(clave, "").startswith("auto")


@dataclass
class OrkXML:
    version_formato: str
    creador: str
    referencia: str
    componentes: list[ComponenteXML]

    def por_nombre(self, nombre: str) -> ComponenteXML:
        for c in self.componentes:
            if c.nombre == nombre:
                return c
        raise KeyError(f"no hay componente '{nombre}' en el .ork "
                       f"(hay: {[c.nombre for c in self.componentes]})")


def _xml(ruta: str | Path) -> ET.Element:
    ruta = Path(ruta)
    if zipfile.is_zipfile(ruta):
        with zipfile.ZipFile(ruta) as z:
            nombres = [n for n in z.namelist() if n.endswith(".ork")]
            if len(nombres) != 1:
                raise ValueError(f"{ruta}: se esperaba un único .ork dentro del ZIP, hay {nombres}")
            return ET.fromstring(z.read(nombres[0]))
    return ET.parse(ruta).getroot()


def leer(ruta: str | Path) -> OrkXML:
    raiz = _xml(ruta)
    rocket = raiz.find("rocket")
    comps: list[ComponenteXML] = []

    def visitar(elem: ET.Element, padre: str | None):
        for hijo in elem:
            if hijo.tag == "subcomponents":
                visitar(hijo, padre)
                continue
            nombre_el = hijo.find("name")
            nombre = nombre_el.text if nombre_el is not None else None
            if hijo.tag in ETIQUETAS:
                c = ComponenteXML(tipo=ETIQUETAS[hijo.tag], nombre=nombre, padre=padre)
                for campo in hijo:
                    if campo.tag in ("subcomponents", "finpoints"):
                        continue
                    if campo.tag == "material":
                        c.material = campo.text or ""
                        c.densidad = float(campo.get("density"))
                    elif campo.text is not None:
                        c.campos[campo.tag] = campo.text.strip()
                    if campo.tag in ("axialoffset", "position"):
                        c.campos[campo.tag + "_method"] = campo.get("method") or campo.get("type") or ""
                fp = hijo.find("finpoints")
                if fp is not None:
                    c.puntos = [(float(p.get("x")), float(p.get("y"))) for p in fp.findall("point")]
                comps.append(c)
            sub = hijo.find("subcomponents")
            if sub is not None:
                visitar(sub, nombre)

    visitar(rocket, None)
    return OrkXML(version_formato=raiz.get("version", ""), creador=raiz.get("creator", ""),
                  referencia=(rocket.findtext("referencetype") or "").strip(), componentes=comps)
