from functools import lru_cache
from typing import Counter, Tuple

from flask import Blueprint, current_app, render_template

from indra_cogex.apps.constants import edge_labels
from indra_cogex.apps.proxies import client

from ...client.queries import get_edge_counter, get_node_counter

__all__ = [
    "home_blueprint",
]

home_blueprint = Blueprint("home", __name__)


@lru_cache(1)
def _get_counters() -> Tuple[Counter, Counter]:
    node_counter = get_node_counter(client=client)
    edge_counter = get_edge_counter(client=client)
    return node_counter, edge_counter


def _figure_number(n: int):
    if n > 1_000_000:
        lead = n / 1_000_000
        if lead < 10:
            return round(lead, 1), "M"
        else:
            return round(lead), "M"
    if n > 1_000:
        lead = n / 1_000
        if lead < 10:
            return round(lead, 1), "K"
        else:
            return round(lead), "K"
    else:
        return n, ""


@home_blueprint.route("/")
def home():
    """Render the home page."""
    node_counter, edge_counter = _get_counters()
    return render_template(
        "home.html",
        format_number=_figure_number,
        node_counter=node_counter,
        edge_counter=edge_counter,
        edge_labels=edge_labels,
        blueprints=current_app.blueprints,
    )
