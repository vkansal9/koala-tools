# Runs a tool's own web handler inside the browser (Pyodide), without a
# server: the handler writes its page into memory and the page is shown.
# The tool scripts themselves are unchanged, so the same file also runs on
# your own computer with:  python3 tools/<name>.py

import importlib
import io


def _browser_urlopen(url, timeout=None):
    """urlopen replacement for the browser (used for name look-ups)."""
    from pyodide.http import open_url
    try:
        text = open_url(url).read()
    except Exception as exc:  # network or cross-origin failure
        raise OSError(str(exc))

    class Reply:
        def read(self):
            return text.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    return Reply()


def run_tool(name, query):
    """The HTML page the tool would serve for '/?<query>'."""
    module = importlib.import_module(name)
    if hasattr(module, "urlopen"):
        module.urlopen = _browser_urlopen

    class Handler(module.Handler):
        def __init__(self, path):  # no socket: just the request path
            self.path = path
            self.wfile = io.BytesIO()

        def send_response(self, code, message=None):
            pass

        def send_header(self, keyword, value):
            pass

        def end_headers(self):
            pass

        def send_error(self, code, message=None, explain=None):
            self.wfile.write(("<p>Error %d</p>" % code).encode())

    handler = Handler("/?" + query if query else "/")
    handler.do_GET()
    return handler.wfile.getvalue().decode("utf-8")
