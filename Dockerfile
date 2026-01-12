FROM continuumio/miniconda3

WORKDIR /app

# Copy environment description
COPY env_docker.yml .

# Create the conda environment
RUN conda env create -f env_docker.yml
RUN pip install pyradiomics

# Activate env automatically
SHELL ["conda", "run", "-n", "pdff", "/bin/bash", "-c"]

# Copy your pipeline code
COPY . .

# Run your pipeline
CMD ["python", "src/pipeline.py"]
