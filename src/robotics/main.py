import asyncio

from robotics.configs.config import settings
from robotics.configs.custom_logging import setup_logger
from robotics.database import close_database_connection, initialize_database_connection
from robotics.models.enums.observation_status import ObservationStatusEnum
from robotics.services.observation_service import claim_next_pending_observation, mark_observation

logger = setup_logger("astro_robotics")


def main() -> None:
    logger.info("Starting Astro Robotics application")
    logger.info(f"Application Name: {settings.app_name}")
    logger.info(f"Application Version: {settings.app_version}")
    logger.info(f"Environment: {settings.environment}")
    logger.info(f"Debug Mode: {settings.debug}")

    asyncio.run(async_main())

    logger.info("Astro Robotics application has stopped")
    logger.info("Exiting main function")


async def async_main() -> None:
    await asyncio.gather(
        normal_processing_cycle(),
        # TODO: Add here the cycles and background tasks for emergency shutdown and other background operations  # noqa: FIX002, TD002
    )


async def normal_processing_cycle() -> None:
    logger.info("Initializing database connection")
    initialize_database_connection()

    try:
        while True:
            try:
                logger.info("Looking for the next pending observation")
                observation = await claim_next_pending_observation()
                if observation is None:
                    logger.info(f"No pending observations found; polling again in {settings.polling_interval} seconds")
                    await asyncio.sleep(settings.polling_interval)
                    continue

                # Process the observation here
                logger.info(f"Processing observation: {observation}")
                try:
                    await process_observation(
                        observation,
                    )  # TODO @dyka3773: Implement the actual observation processing logic  # noqa: FIX002

                    await mark_observation(observation, ObservationStatusEnum.COMPLETED)
                except Exception:
                    logger.exception(f"Error processing observation {observation}")
                    await mark_observation(observation, ObservationStatusEnum.FAILED)
            except Exception:
                logger.exception("Unexpected error in normal processing cycle")
                await asyncio.sleep(settings.polling_interval)
    finally:
        logger.info("Closing database connection")
        await close_database_connection()
        logger.info("Exiting processor function")


async def process_observation(observation) -> None:
    """Process a single observation."""
    await asyncio.sleep(10)


if __name__ == "__main__":
    main()
