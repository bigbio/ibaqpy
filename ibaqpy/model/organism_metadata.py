import json

from dataclasses import dataclass, field
from importlib.resources import open_text
from typing import ClassVar, Optional


@dataclass
class OrganismDescription:
    registry: ClassVar[dict[str, "OrganismDescription"]] = {}

    name: str
    genome_size: int
    histone_proteins: list[str] = field(default_factory=list, repr=False)
    histone_entries: list[str] = field(default_factory=list, repr=False)

    @classmethod
    def get(cls, key, default=None) -> "Optional[OrganismDescription]":
        return cls.registry.get(key.upper(), default)

    @classmethod
    def registered_organisms(cls):
        return cls.registry.keys()

    def __post_init__(self):
        self.registry[self.name.upper()] = self


for v in json.load(open_text("ibaqpy.data", "organisms.json")).values():
    OrganismDescription(**v)