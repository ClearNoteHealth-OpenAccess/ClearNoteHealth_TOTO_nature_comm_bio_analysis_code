FROM rocker/r-ver:4.4.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    zlib1g-dev \
    cmake \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install R libraries
WORKDIR /tmp
COPY ./R/r-requirements.txt /tmp/

RUN Rscript -e 'packages <- readLines("/tmp/r-requirements.txt"); \
    install.packages(packages, repos="https://cloud.r-project.org/")'

# Install python uv
RUN wget -qO- https://astral.sh/uv/install.sh | sh

# Install python libraries
COPY ./uv.lock /tmp/
COPY ./pyproject.toml /tmp/
COPY ./.python-version /tmp/
RUN /root/.local/bin/uv sync
ENV PATH="/tmp/.venv/bin:$PATH"

# Set working directory
WORKDIR /cnh

# Set entrypoint
CMD ["R", "python"]
