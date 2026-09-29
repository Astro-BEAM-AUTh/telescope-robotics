import logging

import inject
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from robotics.database import get_db
from robotics.models.enums.observation_status import ObservationStatusEnum
from robotics.models.observation import Observation
from robotics.utils.time_utils import utc_now

logger = logging.getLogger("astro_robotics")


@inject.params(session=get_db)
async def claim_next_pending_observation(session: AsyncSession) -> Observation | None:
    """Fetch the earliest pending observation for processing if any."""
    result = await session.exec(
        select(Observation).where(Observation.status == ObservationStatusEnum.PENDING).order_by(Observation.created_on.asc()).limit(1),
    )
    observation = result.first()

    if observation is None:
        return None

    await mark_observation(observation, ObservationStatusEnum.IN_PROGRESS, session)

    logger.info(f"Claimed observation {observation.id} for processing")

    return observation


@inject.params(session=get_db)
async def mark_observation(observation: Observation, status: ObservationStatusEnum, session: AsyncSession) -> None:
    """
    Mark an observation with the given status.

    Args:
        observation (Observation): The observation to mark.
        status (ObservationStatusEnum): The status to set for the observation.
        session (AsyncSession): The database session to use for the update.
    """
    observation.status = status
    observation.updated_on = utc_now()
    logger.info(f"Marked observation {observation.id} as {status.name}")
    await session.flush([observation])
