from enum import Enum, auto

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import ClassVar, Iterator, Union, Optional


class QuantificationCategory(Enum):
    TMT = auto()
    ITRAQ = auto()
    LFQ = auto()

    @classmethod
    def from_str(cls, name: str) -> "QuantificationCategory":
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    @classmethod
    def classify(
        cls, labels: Union[list[str], set[str]]
    ) -> tuple["Optional[QuantificationCategory]", "Optional[IsobaricLabel]"]:
        label_scheme = None
        label_category = None

        if len(labels) == 1 and any("label free" in s.lower() for s in labels):
            label_category = cls.LFQ

        elif any("tmt" in s.lower() for s in labels):
            label_category = cls.TMT
            if (
                len(labels) > 11
                or "TMT134N" in labels
                or "TMT133C" in labels
                or "TMT133N" in labels
                or "TMT132C" in labels
                or "TMT132N" in labels
            ):
                label_scheme = IsobaricLabel.TMT16plex
            elif len(labels) == 11 or "TMT131C" in labels:
                label_scheme = IsobaricLabel.TMT11plex
            elif len(labels) > 6:
                label_scheme = IsobaricLabel.TMT10plex
            else:
                label_scheme = IsobaricLabel.TMT6plex

        elif any("itraq" in s.lower() for s in labels):
            label_category = cls.ITRAQ
            if len(labels) > 4:
                label_scheme = IsobaricLabel.ITRAQ8plex
            else:
                label_scheme = IsobaricLabel.ITRAQ4plex

        else:
            raise ValueError(
                f"Cannot infer labeling scheme from {labels}, only support label free, TMT and ITRAQ experiment!"
            )
        return label_category, label_scheme


class IsobaricLabel(Enum):
    TMT6plex = auto()
    TMT10plex = auto()
    TMT11plex = auto()
    TMT16plex = auto()

    ITRAQ4plex = auto()
    ITRAQ8plex = auto()

    @classmethod
    def from_str(cls, name: str) -> "IsobaricLabel":
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def channels(self) -> 'IsobaricLabelSpec':
        return IsobaricLabelSpec.registry[self.name]


@dataclass
class IsobaricLabelSpec(Mapping[str, int]):
    registry: ClassVar[dict[str, "IsobaricLabelSpec"]] = {}

    name: str
    channels: dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        self.registry[self.name] = self

    @property
    def id(self):
        try:
            return IsobaricLabel[self.name]
        except ValueError:
            return None

    def __getitem__(self, key: str) -> int:
        return self.channels[key]

    def __iter__(self) -> Iterator[str]:
        yield from self.channels

    def __len__(self) -> int:
        return len(self.channels)

    def __contains__(self, key) -> bool:
        return key in self.channels

TMT16plex = IsobaricLabelSpec("TMT16plex", {
    "TMT126": 1,
    "TMT127N": 2,
    "TMT127C": 3,
    "TMT128N": 4,
    "TMT128C": 5,
    "TMT129N": 6,
    "TMT129C": 7,
    "TMT130N": 8,
    "TMT130C": 9,
    "TMT131N": 10,
    "TMT131C": 11,
    "TMT132N": 12,
    "TMT132C": 13,
    "TMT133N": 14,
    "TMT133C": 15,
    "TMT134N": 16,
})

TMT11plex = IsobaricLabelSpec("TMT11plex", {
    "TMT126": 1,
    "TMT127N": 2,
    "TMT127C": 3,
    "TMT128N": 4,
    "TMT128C": 5,
    "TMT129N": 6,
    "TMT129C": 7,
    "TMT130N": 8,
    "TMT130C": 9,
    "TMT131N": 10,
    "TMT131C": 11,
})

TMT10plex = IsobaricLabelSpec("TMT10plex", {
    "TMT126": 1,
    "TMT127N": 2,
    "TMT127C": 3,
    "TMT128N": 4,
    "TMT128C": 5,
    "TMT129N": 6,
    "TMT129C": 7,
    "TMT130N": 8,
    "TMT130C": 9,
    "TMT131": 10,
})

TMT6plex = IsobaricLabelSpec("TMT6plex", {
    "TMT126": 1,
    "TMT127": 2,
    "TMT128": 3,
    "TMT129": 4,
    "TMT130": 5,
    "TMT131": 6,
})

ITRAQ4plex = IsobaricLabelSpec("ITRAQ4plex", {"ITRAQ114": 1, "ITRAQ115": 2, "ITRAQ116": 3, "ITRAQ117": 4})

ITRAQ8plex = IsobaricLabelSpec("ITRAQ8plex", {
    "ITRAQ113": 1,
    "ITRAQ114": 2,
    "ITRAQ115": 3,
    "ITRAQ116": 4,
    "ITRAQ117": 5,
    "ITRAQ118": 6,
    "ITRAQ119": 7,
    "ITRAQ121": 8,
})
