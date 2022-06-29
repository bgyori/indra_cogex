# -*- coding: utf-8 -*-

"""Full INDRA CoGEx web app suite."""

import logging
import os

from flask import Flask
from flask_bootstrap import Bootstrap4
from indralab_auth_tools.auth import auth, config_auth

from indra_cogex.apps.constants import INDRA_COGEX_EXTENSION, STATIC_DIR, TEMPLATES_DIR
from indra_cogex.apps.chat_page import chat_blueprint
from indra_cogex.apps.curator import curator_blueprint
from indra_cogex.apps.data_display import data_display_blueprint
from indra_cogex.apps.gla.gene_blueprint import gene_blueprint
from indra_cogex.apps.gla.metabolite_blueprint import metabolite_blueprint
from indra_cogex.apps.home import home_blueprint
from indra_cogex.apps.queries_web import api
from indra_cogex.client.neo4j_client import Neo4jClient
from indra_cogex.client.enrichment.utils import build_caches


logger = logging.getLogger(__name__)


app = Flask(__name__, template_folder=TEMPLATES_DIR, static_folder=STATIC_DIR)
app.register_blueprint(auth)
app.register_blueprint(home_blueprint)
app.register_blueprint(gene_blueprint)
app.register_blueprint(metabolite_blueprint)
app.register_blueprint(data_display_blueprint)
app.register_blueprint(curator_blueprint)
app.register_blueprint(chat_blueprint)
api.init_app(app)

app.extensions[INDRA_COGEX_EXTENSION] = Neo4jClient()

config_auth(app)

# Secret key must be set to use flask-wtf, but there's no *really*
# secure information in this app so it's okay to set randomly
app.config["WTF_CSRF_ENABLED"] = False
app.config["SECRET_KEY"] = os.urandom(32)
app.config["SWAGGER_UI_DOC_EXPANSION"] = "list"
app.config["EXPLAIN_TEMPLATE_LOADING"] = False

bootstrap = Bootstrap4(app)

build_caches()
