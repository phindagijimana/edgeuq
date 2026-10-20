FROM python:3.11-slim-bookworm

# MRtrix3 only provides mrconvert for .mif node images; NIfTI inputs do not need it.
ARG WITH_MRTRIX=1
RUN if [ "$WITH_MRTRIX" = "1" ]; then \
      apt-get update && apt-get install -y --no-install-recommends mrtrix3 && rm -rf /var/lib/apt/lists/*; \
    fi

WORKDIR /opt/edgeuq
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir . && edgeuq --version

ENV OMP_NUM_THREADS=1
WORKDIR /data
ENTRYPOINT ["edgeuq"]
CMD ["--help"]
