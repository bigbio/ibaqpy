from enum import Enum, auto

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import ClassVar, Iterator, Union, Optional


class QuantificationCategory(Enum):
    """
    An enumeration representing different quantification categories used in proteomics.

    Attributes:
        TMT: Represents Tandem Mass Tag quantification.
        ITRAQ: Represents Isobaric Tags for Relative and Absolute Quantitation.
        LFQ: Represents Label-Free Quantification.

    Methods:
        from_str(name: str) -> QuantificationCategory:
            Converts a string to a QuantificationCategory enum member.

        classify(labels: Union[list[str], set[str]]) -> tuple[Optional[QuantificationCategory], Optional[IsobaricLabel]]:
            Classifies a set of labels into a quantification category and determines the isobaric label scheme.
    """

    TMT = auto()
    ITRAQ = auto()
    LFQ = auto()

    @classmethod
    def from_str(cls, name: str) -> "QuantificationCategory":
        """
        Converts a string representation of a quantification category to its corresponding
        QuantificationCategory enum member.

        Parameters
        ----------
            name (str): The name of the quantification category.

        Returns
        -------
            QuantificationCategory: The corresponding enum member.

        Raises
        ------
            KeyError: If the provided name does not match any quantification category.
        """
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    @classmethod
    def classify(
        cls, labels: Union[list[str], set[str]]
    ) -> tuple["Optional[QuantificationCategory]", "Optional[IsobaricLabel]"]:
        """
        Classifies a set of labels into a quantification category and determines the isobaric label scheme.

        Parameters
        ----------
            labels (Union[list[str], set[str]]): A collection of label strings to classify.

        Returns
        -------
            tuple[Optional[QuantificationCategory], Optional[IsobaricLabel]]:
            A tuple containing the quantification category and the isobaric label scheme, if applicable.

        Raises
        ------
            ValueError: If the labels do not correspond to a known quantification category.
        """
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
    """
    An enumeration for different isobaric labeling schemes used in proteomics.

    Attributes
    ----------
        TMT6plex: Represents the TMT 6-plex labeling scheme.
        TMT10plex: Represents the TMT 10-plex labeling scheme.
        TMT11plex: Represents the TMT 11-plex labeling scheme.
        TMT16plex: Represents the TMT 16-plex labeling scheme.
        ITRAQ4plex: Represents the ITRAQ 4-plex labeling scheme.
        ITRAQ8plex: Represents the ITRAQ 8-plex labeling scheme.

    Methods
    -------
        from_str(name: str) -> IsobaricLabel:
            Converts a string to an IsobaricLabel enum member.
        channels() -> IsobaricLabelSpec:
            Retrieves the channel specifications for the isobaric label.
    """

    TMT6plex = auto()
    TMT10plex = auto()
    TMT11plex = auto()
    TMT16plex = auto()

    ITRAQ4plex = auto()
    ITRAQ8plex = auto()

    @classmethod
    def from_str(cls, name: str) -> "IsobaricLabel":
        """
        Converts a string representation of a quantification category to its corresponding
        QuantificationCategory enum member.

        Parameters
        ----------
            name (str): The name of the quantification category.

        Returns
        -------
            QuantificationCategory: The corresponding enum member.

        Raises
        ------
            KeyError: If the provided name does not match any quantification category.
        """
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def channels(self) -> "IsobaricLabelSpec":
        """
        Retrieves the channel specifications associated with the isobaric label.

        Returns:
            IsobaricLabelSpec: The channel specifications for the current isobaric label.
        """
        return IsobaricLabelSpec.registry[self.name]


@dataclass
class IsobaricLabelSpec(Mapping[str, int]):
    """
    A data class representing the specifications of isobaric labels, including their
    name and channel mappings. This class supports dictionary-like access to channel
    information and maintains a registry of all instances.

    Attributes:
        registry (ClassVar[dict[str, IsobaricLabelSpec]]): A class-level registry of all
            isobaric label specifications.
        name (str): The name of the isobaric label.
        channels (dict[str, int]): A mapping of channel names to their respective indices.

    Methods:
        __post_init__(): Registers the instance in the class-level registry.
        id: Retrieves the corresponding IsobaricLabel enum member for the label name.
        __getitem__(key: str) -> int: Returns the index of the specified channel.
        __iter__() -> Iterator[str]: Iterates over the channel names.
        __len__() -> int: Returns the number of channels.
        __contains__(key) -> bool: Checks if a channel name exists in the channels.
    """

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


TMT16plex = IsobaricLabelSpec(
    "TMT16plex",
    {
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
    },
)

TMT11plex = IsobaricLabelSpec(
    "TMT11plex",
    {
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
    },
)

TMT10plex = IsobaricLabelSpec(
    "TMT10plex",
    {
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
    },
)

TMT6plex = IsobaricLabelSpec(
    "TMT6plex",
    {
        "TMT126": 1,
        "TMT127": 2,
        "TMT128": 3,
        "TMT129": 4,
        "TMT130": 5,
        "TMT131": 6,
    },
)

ITRAQ4plex = IsobaricLabelSpec(
    "ITRAQ4plex", {"ITRAQ114": 1, "ITRAQ115": 2, "ITRAQ116": 3, "ITRAQ117": 4}
)

ITRAQ8plex = IsobaricLabelSpec(
    "ITRAQ8plex",
    {
        "ITRAQ113": 1,
        "ITRAQ114": 2,
        "ITRAQ115": 3,
        "ITRAQ116": 4,
        "ITRAQ117": 5,
        "ITRAQ118": 6,
        "ITRAQ119": 7,
        "ITRAQ121": 8,
    },
)
