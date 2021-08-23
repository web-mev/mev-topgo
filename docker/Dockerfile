FROM debian:stretch

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
  apt-get install -y build-essential \
    wget

# Install all software under /opt/software
RUN mkdir -p /opt/software

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  -O /opt/software/miniconda.sh \
  && bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# Install the bioC dependencies:
ADD spec_file.yaml /opt/software/
RUN conda env create --name topgo_env --file=/opt/software/spec_file.yaml

# Add the sctk dimensional reduction script
ADD topgo.R /opt/software
ADD create_gene_mappings.R /opt/software
ADD create_go_mapping_files.R /opt/software

# add to the PATH. Otherwise, conda can swallow return codes
ENV PATH="/opt/software/miniconda/envs/topgo_env/bin:${PATH}"

# Create the gene mapping and GO term files:
mkdir /data
Rscript /opt/software/create_gene_mappings.R human /data/human_gene_mappings.tsv
Rscript /opt/software/create_gene_mappings.R mouse /data/mouse_gene_mappings.tsv
Rscript /opt/software/create_go_mapping_files.R human /data/human_go_mappings.tsv
Rscript /opt/software/create_go_mapping_files.R mouse /data/mouse_go_mappings.tsv

ENTRYPOINT ["/bin/bash"]
