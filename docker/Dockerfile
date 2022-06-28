FROM labsyspharm/indra

ARG INDRA_NEO4J_URL
ARG INDRA_NEO4J_USER
ARG INDRA_NEO4J_PASSWORD

RUN git clone https://github.com/bgyori/indra_cogex.git && \
    cd indra_cogex && \
    # pip install with [web,gunicorn,gsea] \
    pip3 install -e .[web,gunicorn,gsea] && \
    cd .. && \
    python -m indra.ontology.bio build && \
    python -c "from indra_cogex.client.enrichment.utils import build_caches; build_caches();" && \
    # Install nodejs so npm can be used: https://github.com/nodesource/distributions/tree/igsu/16.x-lts#installation-instructions \
    curl -fsSL https://deb.nodesource.com/setup_16.x | sudo -E bash - \
    sudo apt-get install -y nodejs && \
    # Get indralab-vue, install it and then pack it to a local tgz file \
    git clone --single-branch --branch vue3 https://github.com/indralab/indralab-vue.git && \
    cd indralab-vue && \
    npm install && \
    npm run build && \
    npm pack && \
    cd .. && \
    mv indralab-vue/indralab-vue-*.tgz indra_cogex/src/indra_cogex/apps/chat_page/app/ && \
    # Run deployment script
    cd indra_cogex/src/indra_cogex/apps/chat_page/app/ && \
    ilTgz=`ls -1 indralab-vue-*.tgz` && \
    ./deploy.sh -i ${ilTgz} -s /chat

ENTRYPOINT python -m indra_cogex.apps.cli --port 5000 --host "0.0.0.0" --with-gunicorn --workers 4 --timeout 300
