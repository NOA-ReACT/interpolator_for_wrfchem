"""Handles reading of species maps from toml files"""

from dataclasses import dataclass
from pathlib import Path
import tomli

UNITS = {
    "kg": 1,
    "g": 1e-3,
    "mg": 1e-6,
    "ug": 1e-9,
}


@dataclass
class Units:
    source: str
    target: str


class SpeciesMap:
    name: str
    units: Units
    aliases_source: dict[str, str]
    aliases_target: dict[str, str]
    map: dict[str, dict[str, float]]

    required_source_species: list[str]
    target_species: list[str]

    def __init__(self, path: Path | str):
        self.path = path

        with open(path, "rb") as fin:
            data = tomli.load(fin)

        # Check metadata
        if "meta" not in data:
            raise ValueError(f"Missing 'meta' section in {path}")
        for field in ["name", "version"]:
            if field not in data["meta"]:
                raise ValueError(f"Missing '{field}' field in {path}")
        if data["meta"]["version"] != 1:
            raise ValueError(f"Only version 1 species files are supported")
        self.name = data["meta"]["name"]

        # Check units
        if "units" not in data:
            raise ValueError(f"Missing 'units' section in {path}")
        source_unit_name = data["units"]["source"]
        target_unit_name = data["units"]["target"]

        if source_unit_name not in UNITS:
            raise ValueError(f"Unknown source unit {source_unit_name}")
        if target_unit_name not in UNITS:
            raise ValueError(f"Unknown target unit {target_unit_name}")

        self.units = Units(UNITS[source_unit_name], UNITS[target_unit_name])

        # Store alias maps
        self.aliases_source = data.get("aliases_source", {})
        self.aliases_target = data.get("aliases_target", {})

        # Store species map, verify data types
        self.map = data["species_map"]
        for spec in self.map:
            if not isinstance(self.map[spec], dict):
                raise ValueError(f"Species map for {spec} is not a dict")
            for component in self.map[spec]:
                if not isinstance(self.map[spec][component], float):
                    raise ValueError(
                        f"Species coefficient for {spec} -> {component} is not a float"
                    )

        # Compute required source keys
        self.required_source_species = []
        for spec in self.map:
            for component in self.map[spec].keys():
                alias = self.aliases_source.get(component, component)
                if alias not in self.required_source_species:
                    self.required_source_species.append(alias)

        # Compute target species
        self.target_species = []
        for spec in self.map:
            alias = self.aliases_target.get(spec, spec)
            if alias not in self.target_species:
                self.target_species.append(spec)

    def __str__(self) -> str:
        return f"SpeciesMap({self.path}, {self.name}): {', '.join(self.target_species)}"
