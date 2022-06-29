from pathlib import Path
from indra_cogex.apps.wsgi import app


if __name__ == '__main__':
    base_template = app.jinja_env.get_template("chat/chat_page.html")
    index_path = Path("./public/index.html")

    with app.app_context(), app.test_request_context():
        with index_path.open("w") as f:
            f.write(base_template.render())

    print(f"Wrote index.html to {index_path.absolute().as_posix()}")
