"""Curation app for INDRA CoGEx."""

import logging
import time
from typing import Collection

import flask
from flask import Response, redirect, render_template, url_for
from flask_jwt_extended import jwt_optional
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

from indra_cogex.apps import proxies
from indra_cogex.apps.proxies import client
from indra_cogex.client.curation import (
    get_go_curation_hashes,
    get_goa_hashes,
    get_mesh_curation_hashes,
    get_ppi_hashes,
)
from indra_cogex.client.queries import get_stmts_for_stmt_hashes

from .curator.utils import get_unfinished_hashes
from .utils import render_statements, resolve_email

__all__ = [
    "curator_blueprint",
]

logger = logging.getLogger(__name__)
curator_blueprint = flask.Blueprint("curator", __name__, url_prefix="/curate")


class GeneOntologyForm(FlaskForm):
    """A form for choosing a GO term."""

    term = StringField(
        "Gene Ontology Term",
        validators=[DataRequired()],
        description="Choose a gene ontology term to curate (e.g., GO:0003677)",
    )
    submit = SubmitField("Submit")


@curator_blueprint.route("/go/", methods=["GET", "POST"])
def gene_ontology():
    """A home page for GO curation."""
    form = GeneOntologyForm()
    if form.is_submitted():
        return redirect(url_for(f".{curate_go.__name__}", term=form.term.data))
    return render_template("curation/go_form.html", form=form)


@curator_blueprint.route("/go/<term>", methods=["GET"])
@jwt_optional
def curate_go(term: str):
    hashes = get_go_curation_hashes(go_term=("GO", term), client=client)
    return _render_hashes(hashes, title=f"GO Curator: {term}")


class MeshDiseaseForm(FlaskForm):
    """A form for choosing a MeSH disease term."""

    term = StringField(
        "MeSH Term",
        validators=[DataRequired()],
        description="Choose a MeSH disease to curate (e.g., D006009)",
    )
    submit = SubmitField("Submit")


@curator_blueprint.route("/mesh-disease/", methods=["GET", "POST"])
def mesh():
    """A home page for MeSH Disease curation."""
    form = MeshDiseaseForm()
    if form.is_submitted():
        return redirect(url_for(f".{curate_mesh.__name__}", term=form.term.data))
    return render_template("curation/mesh_form.html", form=form)


@curator_blueprint.route("/mesh-disease/<term>", methods=["GET"])
@jwt_optional
def curate_mesh(term: str):
    hashes = get_mesh_curation_hashes(mesh_term=("MESH", term), client=client)
    return _render_hashes(hashes, title=f"MeSH Curator: {term}")


def _render_hashes(
    hashes: Collection[int], title: str, filter_curated: bool = True
) -> Response:
    logger.info(f"Getting statements for {len(hashes)} hashes")
    start_time = time.time()
    stmts, evidence_counts = get_stmts_for_stmt_hashes(
        hashes[: proxies.limit], evidence_limit=10, return_evidence_counts=True
    )
    logger.info(f"Got statements in {time.time() - start_time:.2f} seconds")
    _, _, email = resolve_email()
    return render_statements(
        stmts,
        user_email=email,
        title=title,
        filter_curated=filter_curated,
        evidence_counts=evidence_counts,
    )


@curator_blueprint.route("/ppi", methods=["GET"])
@jwt_optional
def ppi():
    """The PPI curator looks for the highest evidences for PPIs that don't appear in a database."""
    hashes = get_ppi_hashes(client=client, limit=proxies.limit)
    return _render_hashes(hashes, title="PPI Curator")


@curator_blueprint.route("/goa", methods=["GET"])
@jwt_optional
def goa():
    """The GO Annotation curator looks for the highest evidence gene-GO term relations that don't appear in GOA."""
    hashes = get_goa_hashes(client=client, limit=proxies.limit)
    return _render_hashes(hashes, title="GO Annotation Curator")


@curator_blueprint.route("/conflicts", methods=["GET"])
@jwt_optional
def conflicts():
    """Curate statements with conflicting prior curations."""
    hashes = get_unfinished_hashes(client=client)
    return _render_hashes(hashes, title="Conflict Resolver", filter_curated=False)