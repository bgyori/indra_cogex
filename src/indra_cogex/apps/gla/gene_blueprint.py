"""Gene-centric analysis blueprint."""

from pathlib import Path
from typing import Dict, List, Mapping, Tuple

import flask
import pandas as pd
from flask import url_for
from flask_wtf import FlaskForm
from indra.databases import hgnc_client
from wtforms import BooleanField, SubmitField, TextAreaField
from wtforms.validators import DataRequired

from indra_cogex.apps.constants import INDRA_COGEX_WEB_LOCAL
from indra_cogex.apps.proxies import client
from indra_cogex.client.enrichment.continuous import (
    get_human_scores,
    get_mouse_scores,
    indra_downstream_gsea,
    indra_upstream_gsea,
    reactome_gsea,
    wikipathways_gsea,
)

from .fields import (
    alpha_field,
    correction_field,
    file_field,
    indra_path_analysis_field,
    keep_insignificant_field,
    minimum_belief_field,
    minimum_evidence_field,
    permutations_field,
    source_field,
    species_field,
)
from ...client.enrichment.continuous import get_rat_scores, go_gsea
from ...client.enrichment.discrete import (
    EXAMPLE_GENE_IDS,
    go_ora,
    indra_downstream_ora,
    indra_upstream_ora,
    phenotype_ora,
    reactome_ora,
    wikipathways_ora,
)
from ...client.enrichment.signed import (
    EXAMPLE_NEGATIVE_HGNC_IDS,
    EXAMPLE_POSITIVE_HGNC_IDS,
    reverse_causal_reasoning,
)

__all__ = ["gene_blueprint"]

gene_blueprint = flask.Blueprint("gla", __name__, url_prefix="/gene")

genes_field = TextAreaField(
    "Genes",
    description="Paste your list of gene symbols, HGNC gene identifiers, or"
    ' CURIEs here or click here to use <a href="#" onClick="exampleGenes()">an'
    " example list of human genes</a> related to COVID-19.",
    validators=[DataRequired()],
)
positive_genes_field = TextAreaField(
    "Positive Genes",
    description="Paste your list of gene symbols, HGNC gene identifiers, or CURIEs here",
    validators=[DataRequired()],
)
negative_genes_field = TextAreaField(
    "Negative Genes",
    description="Paste your list of gene symbols, HGNC gene identifiers, or"
    ' CURIEs here or click here to use <a href="#" onClick="exampleGenes()">an'
    " example list</a> related to prostate cancer.",
    validators=[DataRequired()],
)


def parse_genes_field(s: str) -> Tuple[Dict[str, str], List[str]]:
    """Parse a gene field string."""
    records = {
        record.strip().strip('"').strip("'").strip()
        for line in s.strip().lstrip("[").rstrip("]").split()
        if line
        for record in line.strip().split(",")
        if record.strip()
    }
    hgnc_ids = []
    errors = []
    for entry in records:
        if entry.lower().startswith("hgnc:"):
            hgnc_ids.append(entry.lower().replace("hgnc:", "", 1))
        elif entry.isnumeric():
            hgnc_ids.append(entry)
        else:  # probably a symbol
            hgnc_id = hgnc_client.get_current_hgnc_id(entry)
            if hgnc_id:
                hgnc_ids.append(hgnc_id)
            else:
                errors.append(entry)
    genes = {hgnc_id: hgnc_client.get_hgnc_name(hgnc_id) for hgnc_id in hgnc_ids}
    return genes, errors


class DiscreteForm(FlaskForm):
    """A form for discrete gene set enrichment analysis."""

    genes = genes_field
    indra_path_analysis = indra_path_analysis_field
    minimum_evidence = minimum_evidence_field
    minimum_belief = minimum_belief_field
    alpha = alpha_field
    correction = correction_field
    keep_insignificant = keep_insignificant_field
    if INDRA_COGEX_WEB_LOCAL:
        local_download = BooleanField("local_download")
    submit = SubmitField("Submit")

    def parse_genes(self) -> Tuple[Mapping[str, str], List[str]]:
        """Resolve the contents of the text field."""
        return parse_genes_field(self.genes.data)


class SignedForm(FlaskForm):
    """A form for signed gene set enrichment analysis."""

    positive_genes = positive_genes_field
    negative_genes = negative_genes_field
    minimum_evidence = minimum_evidence_field
    minimum_belief = minimum_belief_field
    alpha = alpha_field
    # correction = correction_field
    keep_insignificant = keep_insignificant_field
    submit = SubmitField("Submit")

    def parse_positive_genes(self) -> Tuple[Mapping[str, str], List[str]]:
        """Resolve the contents of the text field."""
        return parse_genes_field(self.positive_genes.data)

    def parse_negative_genes(self) -> Tuple[Mapping[str, str], List[str]]:
        """Resolve the contents of the text field."""
        return parse_genes_field(self.negative_genes.data)


class ContinuousForm(FlaskForm):
    """A form for continuous gene set enrichment analysis."""

    file = file_field
    species = species_field
    permutations = permutations_field
    alpha = alpha_field
    keep_insignificant = keep_insignificant_field
    source = source_field
    minimum_evidence = minimum_evidence_field
    minimum_belief = minimum_belief_field
    submit = SubmitField("Submit")

    def get_scores(self) -> Dict[str, float]:
        """Get scores dictionary."""
        name = self.file.data.filename
        sep = "," if name.endswith("csv") else "\t"
        df = pd.read_csv(self.file.data, sep=sep)
        if self.species.data == "rat":
            scores = get_rat_scores(df)
        elif self.species.data == "mouse":
            scores = get_mouse_scores(df)
        elif self.species.data == "human":
            scores = get_human_scores(df)
        else:
            raise ValueError
        return scores


@gene_blueprint.route("/discrete", methods=["GET", "POST"])
def discretize_analysis():
    """Render the home page."""
    form = DiscreteForm()
    if form.validate_on_submit():
        method = form.correction.data
        alpha = form.alpha.data
        keep_insignificant = form.keep_insignificant.data
        minimum_evidence_count = form.minimum_evidence.data
        minimum_belief = form.minimum_belief.data
        genes, errors = form.parse_genes()
        gene_set = set(genes)

        go_results = go_ora(
            client,
            gene_set,
            method=method,
            alpha=alpha,
            keep_insignificant=keep_insignificant,
        )
        wikipathways_results = wikipathways_ora(
            client,
            gene_set,
            method=method,
            alpha=alpha,
            keep_insignificant=keep_insignificant,
        )
        reactome_results = reactome_ora(
            client,
            gene_set,
            method=method,
            alpha=alpha,
            keep_insignificant=keep_insignificant,
        )
        phenotype_results = phenotype_ora(
            gene_set,
            client=client,
            method=method,
            alpha=alpha,
            keep_insignificant=keep_insignificant,
        )
        if form.indra_path_analysis.data:
            indra_upstream_results = indra_upstream_ora(
                client,
                gene_set,
                method=method,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
                minimum_evidence_count=minimum_evidence_count,
                minimum_belief=minimum_belief,
            )
            indra_downstream_results = indra_downstream_ora(
                client,
                gene_set,
                method=method,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
                minimum_evidence_count=minimum_evidence_count,
                minimum_belief=minimum_belief,
            )
        else:
            indra_upstream_results = None
            indra_downstream_results = None

        if INDRA_COGEX_WEB_LOCAL and form.local_download.data:
            downloads = Path.home().joinpath("Downloads")
            go_results.to_csv(
                downloads.joinpath("go_results.tsv"), sep="\t", index=False
            )
            wikipathways_results.to_csv(
                downloads.joinpath("wikipathways_results.tsv"), sep="\t", index=False
            )
            reactome_results.to_csv(
                downloads.joinpath("reactome_results.tsv"), sep="\t", index=False
            )
            phenotype_results.to_csv(
                downloads.joinpath("phenotype_results.tsv"), sep="\t", index=False
            )
            if form.indra_path_analysis.data:
                indra_downstream_results.to_csv(
                    downloads.joinpath("indra_downstream_results.tsv"),
                    sep="\t",
                    index=False,
                )
                indra_upstream_results.to_csv(
                    downloads.joinpath("indra_upstream_results.tsv"),
                    sep="\t",
                    index=False,
                )
            flask.flash(f"Downloaded files to {downloads}")
            return flask.redirect(url_for(f".{discretize_analysis.__name__}"))

        return flask.render_template(
            "gene_analysis/discrete_results.html",
            genes=genes,
            errors=errors,
            method=method,
            alpha=alpha,
            go_results=go_results,
            wikipathways_results=wikipathways_results,
            reactome_results=reactome_results,
            phenotype_results=phenotype_results,
            indra_downstream_results=indra_downstream_results,
            indra_upstream_results=indra_upstream_results,
        )

    return flask.render_template(
        "gene_analysis/discrete_form.html",
        form=form,
        example_hgnc_ids=", ".join(EXAMPLE_GENE_IDS),
    )


@gene_blueprint.route("/signed", methods=["GET", "POST"])
def signed_analysis():
    """Render the signed gene set enrichment analysis form."""
    form = SignedForm()
    if form.validate_on_submit():
        # method = form.correction.data
        # alpha = form.alpha.data
        positive_genes, positive_errors = form.parse_positive_genes()
        negative_genes, negative_errors = form.parse_negative_genes()
        results = reverse_causal_reasoning(
            client=client,
            positive_hgnc_ids=positive_genes,
            negative_hgnc_ids=negative_genes,
            alpha=form.alpha.data,
            keep_insignificant=form.keep_insignificant.data,
            minimum_evidence_count=form.minimum_evidence.data,
            minimum_belief=form.minimum_belief.data,
        )
        return flask.render_template(
            "gene_analysis/signed_results.html",
            positive_genes=positive_genes,
            positive_errors=positive_errors,
            negative_genes=negative_genes,
            negative_errors=negative_errors,
            results=results,
            # method=method,
            # alpha=alpha,
        )
    return flask.render_template(
        "gene_analysis/signed_form.html",
        form=form,
        example_positive_hgnc_ids=", ".join(EXAMPLE_POSITIVE_HGNC_IDS),
        example_negative_hgnc_ids=", ".join(EXAMPLE_NEGATIVE_HGNC_IDS),
    )


@gene_blueprint.route("/continuous", methods=["GET", "POST"])
def continuous_analysis():
    """Render the continuous analysis form."""
    form = ContinuousForm()
    if form.validate_on_submit():
        scores = form.get_scores()
        source = form.source.data
        alpha = form.alpha.data
        permutations = form.permutations.data
        keep_insignificant = form.keep_insignificant.data
        if source == "go":
            results = go_gsea(
                client=client,
                scores=scores,
                permutation_num=permutations,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
            )
        elif source == "wikipathways":
            results = wikipathways_gsea(
                client=client,
                scores=scores,
                permutation_num=permutations,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
            )
        elif source == "reactome":
            results = reactome_gsea(
                client=client,
                scores=scores,
                permutation_num=permutations,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
            )
        elif source == "indra-upstream":
            results = indra_upstream_gsea(
                client=client,
                scores=scores,
                permutation_num=permutations,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
                minimum_evidence_count=form.minimum_evidence.data,
                minimum_belief=form.minimum_belief.data,
            )
        elif source == "indra-downstream":
            results = indra_downstream_gsea(
                client=client,
                scores=scores,
                permutation_num=permutations,
                alpha=alpha,
                keep_insignificant=keep_insignificant,
                minimum_evidence_count=form.minimum_evidence.data,
                minimum_belief=form.minimum_belief.data,
            )
        else:
            raise ValueError

        return flask.render_template(
            "gene_analysis/continuous_results.html",
            source=source,
            results=results,
        )
    return flask.render_template(
        "gene_analysis/continuous_form.html",
        form=form,
    )
