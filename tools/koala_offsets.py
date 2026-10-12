#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# KOALA: Offsets for a given position angle
#

import html
import math
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse, parse_qs

decimals = 2

SCALES = {"WIDE": 1.25, "NARROW": 0.7}  # arcsec/spaxel


def calculate(pa_text, s):
    """Return the (N, E) offsets for c->a, a->b, b->c, b->d."""
    pa = math.radians(float(pa_text))
    pa60 = math.radians(60 - float(pa_text))

    ab = (s * math.sin(pa), -s * math.cos(pa))
    bc = (-s * (math.sin(pa60) + math.sin(pa)),
          s * (math.cos(pa) - math.cos(pa60)))
    ca = (s * math.sin(pa60), s * math.cos(pa60))
    bd = (-s * math.sqrt(3) * math.cos(pa), -s * math.sqrt(3) * math.sin(pa))

    return [tuple(round(x, decimals) for x in pair) for pair in (ca, ab, bc, bd)]


STYLE = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif;
       margin: 24px; color: #111; background: #ececec; }
h1 { font-size: 16px; margin: 0 0 16px 0; }
table { border-collapse: separate; border-spacing: 10px 8px; }
td { white-space: nowrap; text-align: center; }
td.l { text-align: left; }
input { width: 7em; font-size: 14px; padding: 3px 4px; }
button { font-size: 14px; padding: 5px 14px; }
.y { background: yellow; padding: 4px 8px; }
.err { color: #b00020; font-weight: bold; }
"""

ROWS = ["Offset c ---&gt; a :", "Offset a ---&gt; b :",
        "Offset b ---&gt; c :", "Offset b ---&gt; d :"]


def render(pa_text, scale, result, error):
    if result:
        cells = ["<td>%s</td><td>%s</td>" % pair for pair in result]
    else:
        cells = ["<td>-</td><td>-</td>"] * 4
    if error:
        cells[1] = '<td colspan="2" class="err">%s</td>' % html.escape(error)

    rows = "\n".join('<tr><td class="l">%s</td>%s</tr>' % (label, c)
                     for label, c in zip(ROWS, cells))

    return """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>KOALA: Offsets for a given position angle</title>
<style>%s</style></head>
<body>
<h1>KOALA: Offsets for a given position angle</h1>
<form method="get" action="/">
<table>
<tr><td class="l"><span class="y">Position Angle (PA):</span></td>
    <td colspan="2"><input name="pa" value="%s" autocomplete="off" autofocus></td>
    <td><button type="submit" name="mode" value="WIDE">WIDE</button></td>
    <td><button type="submit" name="mode" value="NARROW">NARROW</button></td></tr>
<tr><td></td><td>N</td><td>E</td><td>%s</td><td>arcsec/spaxel</td></tr>
%s
</table>
</form>
</body></html>""" % (STYLE, html.escape(pa_text, quote=True),
                     scale if scale is not None else "-", rows)


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        url = urlparse(self.path)
        if url.path != "/":
            self.send_error(404)
            return

        query = parse_qs(url.query, keep_blank_values=True)
        pa_text = query.get("pa", [""])[0]
        # Pressing Enter in the PA box (no button) behaves like WIDE.
        mode = query.get("mode", ["WIDE"])[0]
        scale = SCALES.get(mode, SCALES["WIDE"]) if query else None

        result, error = None, ""
        if query:  # WIDE or NARROW was pressed
            try:
                result = calculate(pa_text, scale)
            except ValueError:
                error = "This has to be a number!"

        page = render(pa_text, scale, result, error).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(page)))
        self.end_headers()
        self.wfile.write(page)

    def log_message(self, *args):
        pass  # keep the terminal quiet


def main():
    server = HTTPServer(("127.0.0.1", 0), Handler)  # port 0 = any free port
    address = "http://127.0.0.1:%d/" % server.server_address[1]
    print("KOALA offsets calculator running at " + address)
    print("Press Ctrl+C here to quit.")
    webbrowser.open(address)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nBye.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
