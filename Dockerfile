FROM mambaorg/micromamba

ADD env.yml /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ADD ./scripts/executable_scripts/generate_annotations.sh /usr/local/bin/generate_annotations.sh
ADD ./scripts/executable_scripts/generate_annotations_ABC.sh /usr/local/bin/generate_annotations_ABC.sh