from more_click import make_web_command

from indra_cogex.apps.wsgi import app

cli = make_web_command(app=app)

if __name__ == "__main__":
    cli()
