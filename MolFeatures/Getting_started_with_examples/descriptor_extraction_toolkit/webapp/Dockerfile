FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    Getting_started_with_examples/descriptor_extraction_toolkit/webapp/environment.yml \
    /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace/MolFeatures

ENV DESCRIPYTOR_ROOT=/workspace/MolFeatures
ENV PYTHONPATH=/workspace/MolFeatures:/workspace/MolFeatures/Getting_started_with_examples/descriptor_extraction_toolkit
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
ENV XTBPATH=/opt/conda/share/xtb
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /workspace/MolFeatures/Getting_started_with_examples/descriptor_extraction_toolkit/webapp
EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD /opt/conda/bin/python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=5).read()" || exit 1

CMD ["/opt/conda/bin/streamlit", "run", "app.py"]
