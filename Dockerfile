FROM ghcr.io/astral-sh/uv:bookworm-slim

# Install the project into `/app`
WORKDIR /app

COPY . /app

# Ensure TLS trust store is present for outbound HTTPS requests.
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []

CMD ["uv", "run", "robotics"]