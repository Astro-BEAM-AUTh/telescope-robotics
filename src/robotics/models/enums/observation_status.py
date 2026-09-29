from enum import StrEnum


class ObservationStatusEnum(StrEnum):
    """Enumeration for the status of an observation."""

    PENDING = "PENDING"
    IN_PROGRESS = "IN_PROGRESS"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
