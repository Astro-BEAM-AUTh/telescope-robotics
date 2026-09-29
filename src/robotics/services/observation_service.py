import logging

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from robotics.database import inject_db_session
from robotics.models.enums.observation_status import ObservationStatusEnum
from robotics.models.observation import Observation
from robotics.utils.time_utils import utc_now

logger = logging.getLogger("astro_robotics")


@inject_db_session
async def claim_next_pending_observation(session: AsyncSession) -> Observation | None:
    """Fetch the earliest pending observation for processing if any."""
    result = await session.exec(
        select(Observation).where(Observation.status == ObservationStatusEnum.PENDING).order_by(Observation.created_on.asc()).limit(1),
    )
    observation = result.first()

    if observation is None:
        logger.info("No pending observations found")
        return None

    logger.info(f"Found pending observation {observation.id}")

    await mark_observation(observation, ObservationStatusEnum.IN_PROGRESS)

    logger.info(f"Claimed observation {observation.id} for processing")

    return observation


@inject_db_session
async def mark_observation(observation: Observation, status: ObservationStatusEnum, session: AsyncSession) -> None:
    """
    Mark an observation with the given status.

    Args:
        observation (Observation): The observation to mark.
        status (ObservationStatusEnum): The status to set for the observation.
        session (AsyncSession): The database session to use for the update.
    """
    observation.status = status
    now = utc_now()
    observation.updated_on = now

    if status == ObservationStatusEnum.COMPLETED:
        observation.completed_on = now

    logger.info(f"Marked observation {observation.id} as {status.name}")
    observation = await session.merge(observation)
    await session.flush([observation])
