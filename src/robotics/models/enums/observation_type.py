from enum import StrEnum


class ObservationTypeEnum(StrEnum):
    """Enumeration for the type of an observation."""

    HOT_CALIBRATION = "HOT_CALIBRATION"
    COLD_CALIBRATION = "COLD_CALIBRATION"
    TARGET_OBSERVATION = "TARGET_OBSERVATION"
