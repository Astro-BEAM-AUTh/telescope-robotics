"""Telescope observation database model."""

from datetime import datetime

import sqlalchemy as sa
from pydantic import AnyUrl
from sqlalchemy import event
from sqlmodel import CheckConstraint, Field, Relationship, SQLModel, String, Text
from sqlmodel._compat import SQLModelConfig

from robotics.models.enums.frequencies import BandwidthEnum, CentralFrequencyEnum
from robotics.models.enums.observation_status import ObservationStatusEnum
from robotics.models.enums.observation_type import ObservationTypeEnum
from robotics.models.enums.reference_frame import ReferenceFrameEnum
from robotics.models.user import User
from robotics.utils.time_utils import utc_now


class ObservationBase(SQLModel):
    """Base model for telescope observations shared between API and database."""

    # Observation target information
    target_name: str | None = Field(description="Name of the observation target", sa_type=String(255))

    # Telescope pointing information
    ra: float = Field(description="Right Ascension in degrees")
    dec: float = Field(description="Declination in degrees")

    # Telescope configuration
    bandwidth: BandwidthEnum = Field(
        description="Bandwidth in MHz",
        default=BandwidthEnum.BW_1_5_MHZ,
        sa_column_kwargs={"server_default": sa.text("'BW_1_5_MHZ'")},
    )
    center_frequency: CentralFrequencyEnum = Field(
        description="Center frequency in MHz",
        default=CentralFrequencyEnum.FREQ_1420_MHZ,
        sa_column_kwargs={"server_default": sa.text("'FREQ_1420_MHZ'")},
    )
    velocity_frame: ReferenceFrameEnum = Field(
        description="Velocity reference frame",
        default=ReferenceFrameEnum.LSRK,
        sa_column_kwargs={"server_default": sa.text("'LSRK'")},
    )
    observation_type: ObservationTypeEnum = Field(
        description="Type of observation",
        default=ObservationTypeEnum.TARGET_OBSERVATION,
        sa_column_kwargs={"server_default": sa.text("'TARGET_OBSERVATION'")},
    )
    fft_size: int = Field(default=1024, description="FFT size for spectral observations", sa_column_kwargs={"server_default": sa.text("1024")})
    integration_time: float = Field(description="Integration time in seconds")

    # Planning information
    planned_start: datetime | None = Field(default=None, description="Planned start time for the observation")

    # Output information
    output_filename: str = Field(description="Filename for the observation data", sa_type=String(1000))
    receive_csv: bool = Field(
        default=False,
        description="Whether the user wants to receive a CSV containing frequency-power pairs",
        sa_column_kwargs={"server_default": sa.false()},
    )
    perform_data_analysis: bool = Field(
        default=True,
        description="Whether to perform data analysis and generate plots",
        sa_column_kwargs={"server_default": sa.true()},
    )

    model_config: SQLModelConfig = {
        "json_schema_extra": {
            "example": {
                "target_name": "M31",
                "ra": 10.68470833,
                "dec": 41.26875,
                "bandwidth": 20.0,
                "center_frequency": 1420.0,
                "velocity_frame": "LSRK",
                "observation_type": "TARGET_OBSERVATION",
                "fft_size": 1024,
                "integration_time": 600.0,
                "planned_start": "2026-11-15T20:00:00",
                "output_filename": "m31_observation",
                "receive_csv": False,
                "perform_data_analysis": True,
            },
        },
    }


class Observation(ObservationBase, table=True):
    """Database model for telescope observations."""

    __tablename__ = "observations"
    __table_args__ = (
        CheckConstraint("ra >= 0 AND ra < 360", name="ck_observations_ra_range"),
        CheckConstraint("dec >= -90 AND dec <= 90", name="ck_observations_dec_range"),
        CheckConstraint("integration_time > 0", name="ck_observations_integration_time_positive"),
    )

    id: int | None = Field(default=None, primary_key=True)

    # User relationship
    user_id: int = Field(foreign_key="users.id", description="ID of the user who submitted the observation")
    user: User = Relationship(back_populates="observations")

    # Status tracking
    status: ObservationStatusEnum = Field(
        default=ObservationStatusEnum.PENDING,
        description="Current observation status",
        sa_column_kwargs={"server_default": sa.text("'PENDING'")},
    )
    completed_on: datetime | None = Field(default=None, description="Timestamp of completion")

    # Additional metadata
    created_on: datetime = Field(
        default_factory=utc_now,
        description="Record creation timestamp",
        sa_column_kwargs={"server_default": sa.func.now()},
    )
    updated_on: datetime = Field(
        default_factory=utc_now,
        description="Record update timestamp",
        sa_column_kwargs={"server_default": sa.func.now()},
    )

    csv_download_url: AnyUrl | None = Field(
        default=None,
        description="Pre-signed URL for downloading the CSV file (if receive_csv is True)",
        sa_type=Text,
    )
    analysis_results_url: AnyUrl | None = Field(
        default=None,
        description="Pre-signed URL for downloading the analysis results (if perform_data_analysis is True)",
        sa_type=Text,
    )
    data_download_url: AnyUrl | None = Field(
        default=None,
        description="Pre-signed URL for downloading the observation data",
        sa_type=Text,
    )


@event.listens_for(Observation, "before_update")
def update_updated_on(_: object, __: object, target: Observation) -> None:
    """Automatically update the 'updated_on' timestamp on record update."""
    target.updated_on = utc_now()
