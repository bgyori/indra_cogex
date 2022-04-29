# -*- coding: utf-8 -*-

"""Run the HPOA processor using ``python -m indra_cogex.sources.hpoa``."""

import click
from more_click import verbose_option

from . import HpDiseasePhenotypeProcessor, HpPhenotypeGeneProcessor


@click.command()
@verbose_option
@click.pass_context
def _main(ctx: click.Context):
    ctx.invoke(HpDiseasePhenotypeProcessor.get_cli())
    ctx.invoke(HpPhenotypeGeneProcessor.get_cli())


if __name__ == "__main__":
    _main()
