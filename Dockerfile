# Use Python 3.13 slim as base image
FROM python:3.13-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV DEBIAN_FRONTEND=noninteractive

# Create a non-root user
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Install system dependencies and bioinformatics tools
RUN apt-get update && apt-get install -y \
    # Basic utilities
    curl \
    wget \
    git \
    vim \
    less \
    gzip \
    # Build tools for Python packages
    build-essential \
    gcc \
    g++ \
    make \
    # Libraries needed for pysam and other bioinformatics packages
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    # Additional dependencies
    ca-certificates \
    gnupg \
    # Bioinformatics tools
    samtools \
    bcftools \
    tabix \
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set the working directory
WORKDIR /workspace

# Configure Poetry: install packages globally, don't create virtual env
ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_CREATE=false \
    POETRY_CACHE_DIR='/tmp/poetry_cache' \
    POETRY_VERSION=1.8.5 \
    PYTHONPATH=/workspace

# Install Poetry using pip (we'll use Poetry for all other package management)
RUN pip install poetry==1.8.5

# Copy poetry files
COPY pyproject.toml poetry.lock* ./

# Install dependencies globally (as root)
RUN poetry install --with dev --no-root

# Install Jupyter kernel system-wide (as root)
# RUN python -m ipykernel install --name=python3 --display-name="Python 3" --prefix=/usr/local

# Change ownership of workspace to vscode user
RUN chown -R $USERNAME:$USERNAME /workspace

# Switch to non-root user
USER $USERNAME

# Set up shell for better experience
RUN echo 'alias ll="ls -alF"' >> ~/.bashrc \
    && echo 'alias la="ls -A"' >> ~/.bashrc \
    && echo 'alias l="ls -CF"' >> ~/.bashrc

# Verify bioinformatics tools are installed
RUN samtools --version && bcftools --version && tabix --version

# Verify Python packages are accessible
RUN python -c "import pandas, pysam, matplotlib, seaborn, jupyter; print('✅ All packages accessible globally')"

ENV PYTHONPATH=/workspace:/workspace/src

# Set default command
CMD ["/bin/bash"]
