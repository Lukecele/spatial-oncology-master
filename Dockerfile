# Multi-stage reproducible Docker container for spatial-sanctuary (Osimertinib Master Study)
# Author: Luca Celebrano (University of Naples Federico II)
# License: GNU General Public License v3.0 (GPL-3.0)

FROM python:3.11-slim as base

# Prevent Python from writing .pyc and enable unbuffered output
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# Install minimal build and runtime libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy dependency manifests
COPY pyproject.toml setup.py requirements.txt README.md /app/

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code and tests
COPY src/ /app/src/
COPY data/ /app/data/
COPY tests/ /app/tests/

# Install spatial-sanctuary package in production mode
RUN pip install --no-cache-dir .

# Create unprivileged user
RUN useradd -m -u 1000 appuser && \
    mkdir -p /app/results && \
    chown -R appuser:appuser /app
USER appuser

# Healthcheck / validation test execution
RUN python3 -m unittest discover -s tests -v

ENTRYPOINT ["spatial-sanctuary"]
CMD ["simulate"]
