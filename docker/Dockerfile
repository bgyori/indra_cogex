FROM labsyspharm/indra

ARG INDRA_NEO4J_URL
ARG INDRA_NEO4J_USER
ARG INDRA_NEO4J_PASSWORD

RUN pip3 install git+https://github.com/bgyori/indra_cogex.git#egg=indra_cogex[web,gunicorn,gsea] && \
    python -m indra.ontology.bio build && \
    python -c "from indra_cogex.client.enrichment.utils import build_caches; build_caches();"

ENTRYPOINT python -m indra_cogex.apps.cli --port 5000 --host "0.0.0.0" --with-gunicorn --workers 4 --timeout 300
