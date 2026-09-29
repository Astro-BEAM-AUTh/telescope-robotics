"""Application configuration using Pydantic Settings."""

from importlib import metadata

from pydantic import Field, PostgresDsn
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=False,
        extra="ignore",
    )

    app_name: str = Field(default="Astro BEAM Robotics", description="Application name")
    app_version: str = Field(default=metadata.version("robotics"), description="Application version")
    environment: str = Field(default="DEV", description="Application environment")
    debug: bool = Field(default=True, description="Debug mode")  # Only for DEV

    polling_interval: int = Field(default=10, description="Polling interval in seconds")

    # Database settings
    database_url: PostgresDsn = Field(
        default="postgresql+asyncpg://postgres:postgres@localhost:5432/astro_beam",
        description="PostgreSQL database URL",
    )
    db_echo: bool = Field(default=True, description="Echo SQL queries")  # Only for DEV


# Global settings instance
settings = Settings()
