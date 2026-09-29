from enum import Enum


class CentralFrequencyEnum(Enum):
    """Enumeration for central frequencies."""

    FREQ_1420_MHZ = 1420.0
    FREQ_1670_MHZ = 1670.0
    FREQ_22000_MHZ = 22000.0


class BandwidthEnum(Enum):
    """Enumeration for observation bandwidths around the central frequencies."""

    BW_1_5_MHZ = 1.5
    BW_1_75_MHZ = 1.75
    BW_2_5_MHZ = 2.5
    BW_2_75_MHZ = 2.75
    BW_3_MHZ = 3.0
    BW_3_84_MHZ = 3.84
    BW_5_5_MHZ = 5.5
    BW_6_MHZ = 6.0
    BW_7_MHZ = 7.0
    BW_8_75_MHZ = 8.75
    BW_10_MHZ = 10.0
    BW_12_MHZ = 12.0
    BW_14_MHZ = 14.0
    BW_20_MHZ = 20.0
    BW_28_MHZ = 28.0
