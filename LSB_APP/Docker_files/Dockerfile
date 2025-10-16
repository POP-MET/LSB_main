# Start from Node.js Debian base (official n8n image is Node-based)
FROM node:20-bullseye

# Set working directory
WORKDIR /home/node

# Install bash, curl, wget, git, and Python
RUN apt-get update && apt-get install -y \
    bash \
    curl \
    wget \
    git \
    python3 \
    python3-pip \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Detect architecture and install appropriate Miniconda
ARG ARCH
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"; \
    elif [ "$ARCH" = "aarch64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    wget $MINICONDA_URL -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Install n8n globally
RUN npm install -g n8n

# Expose n8n port and basic auth
ENV N8N_PORT=5678
ENV N8N_BASIC_AUTH_ACTIVE=true
ENV N8N_BASIC_AUTH_USER=testuser
ENV N8N_BASIC_AUTH_PASSWORD=testpassword


RUN mkdir -p /home/node/.n8n && chown -R node:node /home/node/.n8n

VOLUME ["/home/node/.n8n"]

# Switch to non-root user
USER node

# Default command
CMD ["n8n", "start"]
