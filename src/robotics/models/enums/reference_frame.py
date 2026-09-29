from enum import StrEnum


class ReferenceFrameEnum(StrEnum):
    """
    Enumeration for the velocity reference frame.

    - TOPO: Topocentric frame, the raw, uncorrected frame of our antenna sitting exactly where it is on Earth. Used for calibration, debugging, and observing satellites.
    - LSRK: Local Standard of Rest Kinematic frame, a theoretical frame moving in a perfect, circular orbit around the center of the Milky Way. It removes the Earth's spin, the Earth's orbit, and the Sun's messy localized wobbling. This should be our software's default setting, it is the standard for galactic radio astronomy. Anytime a user observes 1.42 GHz Milky Way Hydrogen or local 22 GHz galactic water masers, they must use LSRK.
    - BARY: Barycentric frame, the center of mass of our Solar System. It removes the Earth's orbit and spin, but keeps the Sun's movement. Used for extragalactic observations. If the user points the telescope at a quasar billions of light-years away, or a pulsar in another galaxy, the local spinning of our Milky Way doesn't matter. Barycentric is the standard for anything outside our galaxy.
    - HELIO: Heliocentric frame, the exact center of the Sun. It is incredibly close mathematically to Barycentric (since the Sun holds most of the Solar System's mass). However, many older optical astronomy catalogs and specific pulsar databases publish their coordinates and velocities in the Heliocentric frame. We include this for compatibility with those specific databases.
    - GEO: Geocentric frame, the center of the Earth. It removes the Earth's daily rotation, but keeps the Earth's orbit around the Sun. Used when observing the Moon, near-Earth asteroids, or geostationary satellites.
    """

    TOPO = "TOPO"  # Topocentric: The raw, uncorrected frame of our antenna sitting exactly where it is on Earth.
    # When should the user select it: When calibrating the telescope, debugging software, or looking for man-made interference
    #                         (RFI) and artificial satellites. If they are looking at a satellite, they want TOPO
    #                         because the satellite orbits the Earth along with them.
    LSRK = "LSRK"  # Local Standard of Rest Kinematic: A theoretical frame moving in a perfect, circular orbit around
    #                                                  the center of the Milky Way. It removes the Earth's spin, the
    #                                                  Earth's orbit, and the Sun's messy localized wobbling.
    # When should the user select it: This should be our software's **default** setting, it is the standard for galactic radio astronomy.
    #                         Anytime a user observes 1.42 GHz Milky Way Hydrogen or local 22 GHz galactic water masers,
    #                         they must use LSRK.
    BARY = "BARY"  # Barycentric: The center of mass of our Solar System. It removes the Earth's orbit and spin, but keeps the Sun's movement.
    # When should the user select it: For extragalactic observations. If the user points the telescope at a quasar billions of
    #                         light-years away, or a pulsar in another galaxy, the local spinning of our Milky Way doesn't
    #                         matter. Barycentric is the standard for anything outside our galaxy.
    HELIO = "HELIO"  # Heliocentric: The exact center of the Sun.
    # When should the user select it: It is incredibly close mathematically to Barycentric (since the Sun holds most of the
    #                                 Solar System's mass). However, many older optical astronomy catalogs and specific
    #                                 pulsar databases publish their coordinates and velocities in the Heliocentric frame.
    #                                 We include this for compatibility with those specific databases.
    GEO = "GEO"  # Geocentric: The center of the Earth. It removes the Earth's daily rotation, but keeps the Earth's orbit around the Sun.
    # When should the user select it: When observing the Moon, near-Earth asteroids, or geostationary satellites.
