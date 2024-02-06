FROM mambaorg/micromamba

# for building the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# Build the environment based on the yml file
RUN micromamba install -y -n base -f /tmp/environment.yml \
  && micromamba clean --all --yes

# Add the topGO scripts
ADD topgo.R /usr/local/bin
ADD create_gene_mappings.R /usr/local/bin
ADD create_go_mapping_files.R /usr/local/bin
ADD run.sh /usr/local/bin

