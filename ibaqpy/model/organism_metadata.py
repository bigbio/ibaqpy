import json

from dataclasses import dataclass, field
from importlib.resources import open_text
from typing import ClassVar, Optional


@dataclass
class OrganismDescription:
    """
    Represents an organism's metadata, including its name, genome size, and histone proteins.

    Attributes
    ----------
        registry (ClassVar[dict[str, OrganismDescription]]): A class-level dictionary to store
            registered organisms by their name in uppercase.
        name (str): The name of the organism.
        genome_size (int): The size of the organism's genome.
        histone_proteins (list[str]): A list of histone proteins associated with the organism.
        histone_entries (list[str]): A list of histone entries associated with the organism.

    Methods
    -------
        get(key, default=None) -> Optional[OrganismDescription]: Retrieves an organism description
            from the registry using the given key, returning the default if not found.
        registered_organisms(): Returns a list of all registered organism names.
    """

    registry: ClassVar[dict[str, "OrganismDescription"]] = {}

    name: str
    genome_size: int
    histone_proteins: list[str] = field(default_factory=list, repr=False)
    histone_entries: list[str] = field(default_factory=list, repr=False)

    @classmethod
    def get(cls, key, default=None) -> "Optional[OrganismDescription]":
        """
        Retrieve an organism description from the registry using the given key.

        Parameters:
            key (str): The name of the organism to retrieve, case-insensitive.
            default (Optional[OrganismDescription]): The value to return if the organism is not found.

        Returns:
            Optional[OrganismDescription]: The organism description if found, otherwise the default value.
        """
        return cls.registry.get(key.upper(), default)

    @classmethod
    def registered_organisms(cls):
        """
        Return a list of all registered organism names.

        Returns:
            dict_keys: A view of the registry's keys representing organism names.
        """
        return cls.registry.keys()

    def __post_init__(self):
        """
        Automatically register the organism in the class-level registry
        using its name in uppercase as the key.
        """
        self.registry[self.name.upper()] = self


for v in json.load(open_text("ibaqpy.data", "organisms.json")).values():
    OrganismDescription(**v)
