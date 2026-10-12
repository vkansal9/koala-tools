#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Python script to estimate AAT focus position using 2 FWHM values
# By Angel R. Lopez-Sanchez (AAO/MQU)
# Version 1.1 - 27 Feb 2018
#


import html
import math
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse, parse_qs

NROWS = 12

# name -> (button label, scale, message printed in the terminal)
SCALES = {
    "KOALA": ("KOALA", 1.0, "KOALA scale selected, FWHM given in arcsec"),
    "BIN1": ("TV 1x1", 0.16,
             "TV DIRECT binning 1 scale selected, FWHM given in 0.16 arcsec / pix "),
    "BIN2": ("TV 2x2", 0.32,
             "TV DIRECT binning 2 scale selected, FWHM given in 0.32 arcsec / pix "),
}


# ---------------------------------------------------------------- the maths

def polyfit2(x, y):
    """Least-squares parabola. Returns f(x) (same result as numpy.polyfit(x, y, 2))."""
    n = len(x)
    xm = sum(x) / n                      # centre x: better numerical behaviour
    u = [xi - xm for xi in x]
    s1, s2 = sum(u), sum(ui ** 2 for ui in u)
    s3, s4 = sum(ui ** 3 for ui in u), sum(ui ** 4 for ui in u)
    t0 = sum(y)
    t1 = sum(ui * yi for ui, yi in zip(u, y))
    t2 = sum(ui * ui * yi for ui, yi in zip(u, y))

    def det(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    A = [[n, s1, s2], [s1, s2, s3], [s2, s3, s4]]
    rhs = [t0, t1, t2]
    d = det(A)
    if abs(d) < 1e-12 * max(1.0, s4) ** 2:
        raise ArithmeticError("singular fit")
    c = []
    for k in range(3):
        m = [row[:] for row in A]
        for i in range(3):
            m[i][k] = rhs[i]
        c.append(det(m) / d)
    return lambda xv: c[0] + c[1] * (xv - xm) + c[2] * (xv - xm) ** 2


def compute(rows, scale):
    """rows: list of (focus, fwhm_b, fwhm_r) as typed. Returns a dict of results."""
    focus = [r[0] for r in rows]
    fwhm1 = [r[1] * scale for r in rows]
    fwhm2 = [r[2] * scale for r in rows]
    fwhm = [(b + r) / 2 for b, r in zip(fwhm1, fwhm2)]

    xmax = max(focus) + 0.4
    xmin = min(focus) - 0.4
    xx = [i / 100. * (xmax - xmin) + xmin for i in range(100)]

    out = {"focus": focus, "fwhm": fwhm, "xx": xx, "xmin": xmin, "xmax": xmax}
    for key, values in (("", fwhm), ("b", fwhm1), ("r", fwhm2)):
        f = polyfit2(focus, values)
        fit = [f(xv) for xv in xx]
        best = min(range(len(fit)), key=fit.__getitem__)
        out["fit" + key] = fit
        out["seeing" + key] = round(fit[best], 2)
        out["bestfocus" + key] = round(xx[best], 2)
    return out


# ----------------------------------------------------------------- the plot

def ticks(lo, hi, target=7):
    raw = (hi - lo) / target
    mag = 10 ** math.floor(math.log10(raw))
    step = min((s for s in (1, 2, 2.5, 5, 10) if s * mag >= raw),
               default=10) * mag
    t = math.ceil(lo / step - 1e-9) * step
    out = []
    while t <= hi + 1e-9:
        out.append(round(t, 10))
        t += step
    return out


def fmt(v):
    return ("%.3f" % v).rstrip("0").rstrip(".") if v else "0"


def plot_svg(res):
    W, H = 760, 560
    L, R, T, B = 70, 20, 20, 55
    if res is None:
        return ('<svg viewBox="0 0 %d %d" width="%d" height="%d">'
                '<rect width="100%%" height="100%%" fill="white"/></svg>'
                % (W, H, W, H))

    x0, x1 = res["xmin"], res["xmax"]
    ys = res["fwhm"] + res["fit"] + res["fitb"] + res["fitr"]
    y0, y1 = min(ys), max(ys)
    pad = (y1 - y0) * 0.08 or 0.1
    y0, y1 = y0 - pad, y1 + pad

    def px(x):
        return L + (x - x0) / (x1 - x0) * (W - L - R)

    def py(y):
        return H - B - (y - y0) / (y1 - y0) * (H - T - B)

    def line(yv, color, opacity):
        pts = " ".join("%.1f,%.1f" % (px(a), py(b)) for a, b in zip(res["xx"], yv))
        return ('<polyline points="%s" fill="none" stroke="%s" '
                'stroke-opacity="%s" stroke-width="2"/>' % (pts, color, opacity))

    def square(x, y, color):
        return ('<rect x="%.1f" y="%.1f" width="8" height="8" fill="%s"/>'
                % (px(x) - 4, py(y) - 4, color))

    s = ['<svg viewBox="0 0 %d %d" width="%d" height="%d" '
         'font-family="Helvetica, Arial, sans-serif" font-size="13">' % (W, H, W, H),
         '<rect width="100%" height="100%" fill="white"/>',
         '<clipPath id="c"><rect x="%d" y="%d" width="%d" height="%d"/></clipPath>'
         % (L, T, W - L - R, H - T - B)]

    for t in ticks(x0, x1):
        s.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="black"/>'
                 % (px(t), H - B, px(t), H - B + 5))
        s.append('<text x="%.1f" y="%d" text-anchor="middle">%s</text>'
                 % (px(t), H - B + 20, fmt(t)))
    for t in ticks(y0, y1):
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="black"/>'
                 % (L - 5, py(t), L, py(t)))
        s.append('<text x="%d" y="%.1f" text-anchor="end">%s</text>'
                 % (L - 9, py(t) + 4, fmt(t)))

    s.append('<g clip-path="url(#c)">')
    s.append(line(res["fitb"], "blue", 0.3))
    s.append(line(res["fitr"], "red", 0.3))
    for a, b in zip(res["focus"], res["fwhm"]):
        s.append('<circle cx="%.1f" cy="%.1f" r="10" fill="#1f77b4"/>' % (px(a), py(b)))
    s.append(line(res["fit"], "#ff7f0e", 1))
    s.append(square(res["bestfocusb"], res["seeingb"], "blue"))
    s.append(square(res["bestfocusr"], res["seeingr"], "red"))
    s.append(square(res["bestfocus"], res["seeing"], "green"))
    s.append('</g>')

    s.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="black"/>'
             % (L, T, W - L - R, H - T - B))
    s.append('<text x="%.1f" y="%d" text-anchor="middle">Focus value [mm]</text>'
             % ((L + W - R) / 2., H - 12))
    s.append('<text transform="translate(18 %.1f) rotate(-90)" text-anchor="middle">'
             'FWHM ["]</text>' % ((T + H - B) / 2.))
    s.append('</svg>')
    return "\n".join(s)


# ----------------------------------------------------------------- the page

STYLE = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif;
       margin: 20px; color: #111; background: #ececec; }
h1 { font-size: 16px; margin: 0 0 12px 0; }
.wrap { display: flex; gap: 14px; align-items: flex-start; flex-wrap: wrap; }
table { border-collapse: collapse; }
th { font-weight: normal; padding: 2px 0 6px 0; }
td { padding: 1px; }
td input { width: 5.2em; font-size: 14px; padding: 4px; background: yellow;
           border: 1px solid #999; box-sizing: border-box; }
.scales { display: flex; margin-top: 6px; }
.scales label { flex: 1; }
.scales input { position: absolute; opacity: 0; }
.scales span { display: block; text-align: center; padding: 6px 0; font-size: 13px;
               background: white; border: 1px solid #999; cursor: pointer; }
.scales input:checked + span { background: #2e9e2e; color: white; }
.scales input:focus-visible + span { outline: 2px solid #06c; }
button { width: 100%; font-size: 14px; padding: 7px 0; margin-top: 6px; }
.result { background: lightblue; padding: 8px 12px; margin-top: 6px;
          white-space: pre-wrap; min-height: 1.2em; }
.err { color: #b00020; font-weight: bold; }
svg { max-width: 100%; height: auto; border: 1px solid #999; display: block; }
"""


def render(values, mode, res, error):
    def box(name):
        return '<td><input name="%s" value="%s" autocomplete="off"></td>' % (
            name, html.escape(values.get(name, ""), quote=True))

    rows = "\n".join("<tr>%s%s%s</tr>" % (box("f%d" % i), box("b%d" % i), box("r%d" % i))
                     for i in range(1, NROWS + 1))

    radios = "".join(
        '<label><input type="radio" name="scale" value="%s"%s><span>%s</span></label>'
        % (key, " checked" if key == mode else "", SCALES[key][0])
        for key in ("KOALA", "BIN1", "BIN2"))

    if error:
        result = '<span class="err">%s</span>' % html.escape(error)
    elif res:
        result = 'Focus position: %s mm              Best seeing : %s"' % (
            res["bestfocus"], res["seeing"])
        if SCALES[mode][1] == 1.0:
            result += '             Best RED seeing : %s"' % res["seeingr"]
        result = html.escape(result)
    else:
        result = ""

    return """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Focusing KOALA @ AAT</title>
<style>%s</style></head>
<body>
<h1>Focusing KOALA @ AAT</h1>
<form method="get" action="/">
<div class="wrap">
<div>
<table>
<tr><th>Focus</th><th>FWHM b</th><th>FWHM r</th></tr>
%s
</table>
<div class="scales">%s</div>
<button type="submit" name="go" value="1">Compute!</button>
</div>
<div>
%s
<div class="result">%s</div>
</div>
</div>
</form>
</body></html>""" % (STYLE, rows, radios, plot_svg(res), result)


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        url = urlparse(self.path)
        if url.path != "/":
            self.send_error(404)
            return

        query = parse_qs(url.query, keep_blank_values=True)
        values = dict((k, v[0]) for k, v in query.items())
        mode = values.get("scale", "KOALA")
        if mode not in SCALES:
            mode = "KOALA"

        res, error = None, ""
        if query:  # Compute! was pressed
            print(SCALES[mode][2])
            rows = []
            for i in range(1, NROWS + 1):
                try:
                    rows.append(tuple(float(values.get("%s%d" % (c, i), ""))
                                      for c in "fbr"))
                except ValueError:
                    print("No data in row %d" % i)
            try:
                if len(set(r[0] for r in rows)) < 3:
                    raise ArithmeticError
                res = compute(rows, SCALES[mode][1])
            except (ArithmeticError, ValueError, OverflowError):
                error = ("Need at least 3 complete rows (Focus, FWHM b, FWHM r) "
                         "with different focus values.")
                print(error)
            else:
                print(" Focus values = %s" % res["focus"])
                print(" FWHM values = %s" % res["fwhm"])
                print(' Best seeing    = %s "     b = %s "     r = %s "'
                      % (res["seeing"], res["seeingb"], res["seeingr"]))
                print(" Focus position = %s mm    b = %s mm    r = %s mm"
                      % (res["bestfocus"], res["bestfocusb"], res["bestfocusr"]))

        page = render(values, mode, res, error).encode("utf-8")
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
    print("KOALA focus tool running at " + address)
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
