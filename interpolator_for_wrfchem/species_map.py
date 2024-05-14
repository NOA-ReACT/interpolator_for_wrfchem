"""Handles reading of species maps from toml files"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import tomli

SI_PREFIX = {"M": 6, "k": 3, "": 0, "m": -3, "u": -6, "μ": -6, "n": -9}


def check_prefix(prefix: str) -> None:
    if prefix not in SI_PREFIX:
        raise ValueError(f"Unknown SI prefix {prefix}")


def convert_si(value: Any, from_prefix: str, to_prefix: str):
    if from_prefix == to_prefix:
        return value
    if from_prefix not in SI_PREFIX:
        raise ValueError(f"Unknown SI prefix {from_prefix}")
    if to_prefix not in SI_PREFIX:
        raise ValueError(f"Unknown SI prefix {to_prefix}")
    return value * 10 ** (SI_PREFIX[from_prefix] - SI_PREFIX[to_prefix])


@dataclass
class Units:
    global_model: str = field(default="")
    regional_model: str = field(default="")


@dataclass
class Species:
    coeffs: dict[str, float]
    units: Units
    weight: float = field(default=1.0)
    offset: float = field(default=0.0)


class SpeciesMap:
    name: str
    units: Units
    aliases_source: dict[str, str]
    aliases_target: dict[str, str]
    map: dict[str, Species]

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
        if data["meta"]["version"] != 2:
            raise ValueError(f"Only version 2 species files are supported")
        self.name = data["meta"]["name"]

        # Store alias maps
        self.aliases_source = data.get("aliases_source", {})
        self.aliases_target = data.get("aliases_target", {})

        # Parse species map
        self.map = {}
        for target, raw_spec in data["species_map"].items():
            if not isinstance(raw_spec, dict):
                raise ValueError(f"Species map for {target} is not a dict")

            coeffs = {}
            if "coeffs" not in raw_spec:
                raise ValueError(f"Missing 'coeffs' field for {target}")
            for source, coeff in raw_spec["coeffs"].items():
                if not isinstance(coeff, float):
                    raise ValueError(
                        f"Coefficient for {target} -> {source} is not a float"
                    )
                coeffs[source] = coeff

            global_model_unit = None
            regional_model_unit = None
            units = raw_spec.get("units", {"global": "", "regional": ""})
            if "global" in units:
                global_model_unit = units["global"]
                check_prefix(global_model_unit)
            if "regional" in units:
                regional_model_unit = units["regional"]
                check_prefix(regional_model_unit)
            units = Units(
                global_model=global_model_unit, regional_model=regional_model_unit
            )

            self.map[target] = Species(
                coeffs=coeffs,
                units=units,
                weight=raw_spec.get("weight", 1.0),
                offset=raw_spec.get("offset", 0.0),
            )

        # Compute required source keys
        self.required_source_species = []
        for spec in self.map:
            for component in self.map[spec].coeffs.keys():
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
