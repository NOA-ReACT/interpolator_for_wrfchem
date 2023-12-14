"""Handles reading of species maps from toml files"""

from pathlib import Path
import tomli


class SpeciesMap:
    name: str
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
