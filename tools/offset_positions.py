#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# OFFSETS [in arcsec] between 2 positions


import html
import math
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse, parse_qs

decimals = 2

FIELDS = ["ra1h", "ra1m", "ra1s", "dec1d", "dec1m", "dec1s",
          "ra2h", "ra2m", "ra2s", "dec2d", "dec2m", "dec2s"]


def to_degrees(d, m, s):
    """Degrees/minutes/seconds -> decimal degrees (handles -0 degrees too)."""
    sign = -1.0 if d.strip().startswith("-") else 1.0
    return sign * (abs(float(d)) + float(m) / 60. + float(s) / 3600.)


def calculate(v):
    """Return (|dRA|, RA direction, |dDec|, Dec direction) for 1 -> 2."""
    ra1 = float(v["ra1h"]) + float(v["ra1m"]) / 60. + float(v["ra1s"]) / 3600.
    ra2 = float(v["ra2h"]) + float(v["ra2m"]) / 60. + float(v["ra2s"]) / 3600.
    dec1 = to_degrees(v["dec1d"], v["dec1m"], v["dec1s"])
    dec2 = to_degrees(v["dec2d"], v["dec2m"], v["dec2s"])

    avdec = (dec1 + dec2) / 2

    deltadec = round(3600. * (dec2 - dec1), decimals)
    deltara = round(15 * 3600. * (ra2 - ra1) * math.cos(math.radians(avdec)),
                    decimals)

    ra_dir = "West" if deltara < 0 else "East"
    dec_dir = "South" if deltadec < 0 else "North"
    return math.fabs(deltara), ra_dir, math.fabs(deltadec), dec_dir


STYLE = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif;
       margin: 24px; color: #111; background: #ececec; }
h1 { font-size: 16px; margin: 0 0 16px 0; }
table { border-collapse: separate; border-spacing: 6px 8px; }
td { white-space: nowrap; }
input { width: 4.5em; font-size: 14px; padding: 3px 4px; }
input.sec { width: 6em; }
button { font-size: 14px; padding: 5px 14px; margin: 10px 0 14px 6px; }
.y { background: yellow; padding: 4px 10px; min-width: 5em; }
.c { background: cyan; padding: 4px 10px; min-width: 5em; }
.err { color: #b00020; font-weight: bold; margin-left: 6px; }
"""

OPPOSITE = {"East": "West", "West": "East", "North": "South", "South": "North"}


def render(v, result, error):
    def box(name, cls=""):
        value = html.escape(v.get(name, ""), quote=True)
        return '<input name="%s" value="%s" class="%s" autocomplete="off">' % (
            name, value, cls)

    def position_row(n):
        return (
            "<tr><td>POSITION %d:</td><td>RA:</td>"
            "<td>%s</td><td>%s</td><td>%s</td>"
            "<td>DEC:</td><td>%s</td><td>%s</td><td>%s</td></tr>"
            % (n,
               box("ra%dh" % n), box("ra%dm" % n), box("ra%ds" % n, "sec"),
               box("dec%dd" % n), box("dec%dm" % n), box("dec%ds" % n, "sec")))

    if result:
        dra, ra_dir, ddec, dec_dir = result
        cells = (dra, ra_dir, ddec, dec_dir,
                 dra, OPPOSITE[ra_dir], ddec, OPPOSITE[dec_dir])
    else:
        cells = ("",) * 8

    message = '<span class="err">%s</span>' % html.escape(error) if error else ""

    return """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>OFFSETS [in arcsec] between 2 positions</title>
<style>%s</style></head>
<body>
<h1>OFFSETS [in arcsec] between 2 positions</h1>
<form method="get" action="/">
<table>
%s
%s
</table>
<button type="submit">CALCULATE!</button>%s
<table>
<tr><td>1 -&gt; 2 :</td>
    <td class="y">%s</td><td class="y">%s</td>
    <td class="y">%s</td><td class="y">%s</td></tr>
<tr><td>2 -&gt; 1 :</td>
    <td class="c">%s</td><td class="c">%s</td>
    <td class="c">%s</td><td class="c">%s</td></tr>
</table>
</form>
</body></html>""" % ((STYLE, position_row(1), position_row(2), message) + cells)


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        url = urlparse(self.path)
        if url.path != "/":
            self.send_error(404)
            return

        query = parse_qs(url.query, keep_blank_values=True)
        values = dict((name, query.get(name, [""])[0]) for name in FIELDS)

        result, error = None, ""
        if query:  # the CALCULATE! button was pressed
            try:
                result = calculate(values)
            except ValueError:
                error = "This has to be a number!"

        page = render(values, result, error).encode("utf-8")
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
    print("Offsets calculator running at " + address)
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
