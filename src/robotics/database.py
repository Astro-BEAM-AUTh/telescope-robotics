import functools
from collections.abc import AsyncGenerator, Awaitable, Callable
from contextlib import asynccontextmanager

from sqlalchemy.ext.asyncio import AsyncEngine, async_sessionmaker, create_async_engine
from sqlmodel.ext.asyncio.session import AsyncSession

from robotics.configs.config import settings

# Database engine
engine: AsyncEngine | None = None

# Session factory
async_session_maker: async_sessionmaker[AsyncSession] | None = None


def initialize_database_connection() -> None:
    """Initialize database connection and session factory."""
    global engine, async_session_maker  # noqa: PLW0603

    engine = create_async_engine(
        str(settings.database_url),
        echo=settings.db_echo,
        pool_pre_ping=True,
    )

    async_session_maker = async_sessionmaker(
        engine,
        class_=AsyncSession,
        expire_on_commit=False,
        autocommit=False,
        autoflush=False,
    )

    if engine is None or async_session_maker is None:
        msg = "Failed to initialize database connection"
        raise RuntimeError(msg)


async def close_database_connection() -> None:
    """Close database connection."""
    if engine is not None:
        await engine.dispose()


@asynccontextmanager
async def _get_db_session() -> AsyncGenerator[AsyncSession]:
    """
    Get a database session.

    Yields:
        AsyncSession: Database session

    Usage:
        async with get_db_session() as session:
            result = await session.execute(query)
    """
    if async_session_maker is None:
        msg = "Database not initialized. Call initialize_database_connection() first."
        raise RuntimeError(msg)

    async with async_session_maker() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise


async def get_db() -> AsyncGenerator[AsyncSession]:
    """
    Dependency to inject database session.

    Yields:
        AsyncSession: Database session
    """
    async with _get_db_session() as session:
        yield session


def inject_db_session[**P, T](func: Callable[P, Awaitable[T]]) -> Callable[P, Awaitable[T]]:
    """
    Decorate an async function so it receives a ``session`` keyword argument automatically.

    A new session (and its own commit/rollback transaction) is opened only when the caller
    hasn't already supplied one, so nested calls can share a single transaction by passing
    their own ``session`` through explicitly.
    """

    @functools.wraps(func)
    async def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        if "session" in kwargs:
            return await func(*args, **kwargs)

        async with _get_db_session() as session:
            return await func(*args, session=session, **kwargs)

    return wrapper
