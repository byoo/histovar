# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies (needed for compiling some C extensions like pysam if pre-built wheels aren't available)
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy only the dependency files first to leverage Docker cache
COPY requirements.txt pyproject.toml ./

# Install Python dependencies
# We install torch explicitly to ensure compatibility, though it is in requirements
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code
COPY histovar/ ./histovar/
COPY data/ ./data/

# Install the package in editable mode (or standard mode)
RUN pip install --no-cache-dir -e .

# Define the entrypoint
ENTRYPOINT ["histovar"]
