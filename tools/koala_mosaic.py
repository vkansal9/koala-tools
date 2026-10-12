#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# KOALA: mosaic observing sequence for a given position angle
#

import html
import math
import re
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse, parse_qs, urlencode, quote
from urllib.request import urlopen

decimals = 2

SCALES = {"WIDE": 1.25, "NARROW": 0.7}  # arcsec/spaxel

# Spaxel layout: COLS along the long axis, ROWS along the short axis.
# FOV: 50.6" x 27.4" (WIDE) and 28.3" x 15.3" (NARROW).
COLS, ROWS = 40, 25

# Sky angle (N through E) of the long axis at rotator PA = 0.  90 = East-West,
# which matches the offsets of koala_offsets.py.
PA_ZERO = 90

# Sky images shown behind the mosaic (needs internet).
LEGACY_SERVICE = "https://www.legacysurvey.org/viewer/cutout.jpg"
LEGACY_LAYER = "ls-dr10"
LEGACY_VIEWER = "https://www.legacysurvey.org/viewer/"
HIPS_SERVICE = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
HIPS_SURVEY = "CDS/P/DSS2/color"
SESAME = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?"

# Positions (RA, Dec in degrees, J2000) that work without internet.
KNOWN = {"NGC151": (8.5117, -9.7053)}

R3 = math.sqrt(3)
MAX_OVERLAP_LONG = COLS // 2
MAX_OVERLAP_SHORT = (ROWS + 1) // 2
MAX_BLOCK = 15
MAX_MAP_LABELS = 16  # offset labels on the map up to this many fields


# geometry

def frame_angle(pa_deg):
    """Angle (N through E) of the long axis of the field on the sky."""
    return math.radians(pa_deg + PA_ZERO)


def to_sky(x, y, th):
    """Instrument (long, short) -> sky (N, E)."""
    return (x * math.cos(th) - y * math.sin(th), x * math.sin(th) + y * math.cos(th))


def to_inst(n, e, th):
    """Sky (N, E) -> instrument (Wide, Narrow)."""
    return (n * math.cos(th) + e * math.sin(th), -n * math.sin(th) + e * math.cos(th))


def field_size(mode):
    """Full size of one field in arcsec: (Wide, Narrow)."""
    p = SCALES[mode]
    return COLS * p + p / 2, (ROWS - 1) * p * R3 / 2 + 2 * p / R3


def field_outline(pa_deg, mode):
    """Corners (N, E) of one field, relative to its centre."""
    th = frame_angle(pa_deg)
    half_x, half_y = (v / 2 for v in field_size(mode))
    return [to_sky(sx * half_x, sy * half_y, th)
            for sx, sy in ((1, 1), (1, -1), (-1, -1), (-1, 1))]


def at_centre(n, e):
    return abs(n) < 1e-6 and abs(e) < 1e-6


def has_centre_field(fields):
    return any(at_centre(n, e) for n, e in fields)


def galaxy_point(u, v, gal):
    """Map unit-disc coordinates to (N, E) on the galaxy ellipse."""
    a, b, g = gal["a"], gal["b"], math.radians(gal["pa"])
    return (u * a * math.cos(g) - v * b * math.sin(g),
            u * a * math.sin(g) + v * b * math.cos(g))


def spiral_key(cell):
    """Centre first, then ring by ring; each ring starts at +long axis."""
    ci, rj = cell
    angle = math.atan2(-rj, ci)
    if angle < -1e-9:
        angle += 2 * math.pi
    return (round(max(abs(ci), abs(rj)), 3), round(angle, 6))


def snake(cells):
    """Row by row, back and forth."""
    out = []
    for k, rj in enumerate(sorted({rj for _, rj in cells})):
        row = sorted((c for c in cells if c[1] == rj), reverse=True)
        if k % 2:
            row.reverse()
        out += row
    return out


def mosaic(pa_deg, mode, gal, overlap, block=None, order="spiral", grid=160):
    """Fields to observe, in order, and the percent of the galaxy they cover.

    overlap = (spaxels along the long axis, rows along the short axis).
    block = (columns, rows) gives a fixed block of fields centred on the
    galaxy; block = None gives every field needed to cover the galaxy.
    Field centres are (N, E) in arcsec from the galaxy centre.
    """
    p = SCALES[mode]
    th = frame_angle(pa_deg)
    ov_long, ov_short = overlap
    step_x = (COLS - ov_long) * p
    step_y = (ROWS - ov_short) * p * R3 / 2
    odd = (ROWS - ov_short) % 2  # then every other row shifts by half a spaxel

    def shift(j):
        return p / 2 if (odd and j % 2) else 0.0

    samples = []
    for i in range(grid):
        u = -1 + (2 * i + 1) / grid
        for j in range(grid):
            v = -1 + (2 * j + 1) / grid
            if u * u + v * v <= 1:
                samples.append((u, v))
    for i in range(720):  # the edge of the galaxy itself
        t = math.radians(i / 2)
        samples.append((math.cos(t), math.sin(t)))
    points = [to_inst(*galaxy_point(u, v, gal), th) for u, v in samples]

    # Each cell is (column, row) counted from the galaxy centre.
    centre_of = {}
    if block:
        ncols, nrows = block
        for r in range(nrows):
            j = r - (nrows - 1) // 2
            for c in range(ncols):
                cell = (c - (ncols - 1) / 2, r - (nrows - 1) / 2)
                centre_of[cell] = (cell[0] * step_x + shift(j),
                                   cell[1] * step_y)
    else:
        for x, y in points:
            j = round(y / step_y)
            i = round((x - shift(j)) / step_x)
            centre_of[(i, j)] = (i * step_x + shift(j), j * step_y)

    if order == "snake":
        cells = snake(list(centre_of))
    else:
        cells = sorted(centre_of, key=spiral_key)
    centres = [centre_of[c] for c in cells]

    inside = sum(any(abs(x - cx) <= step_x / 2 and abs(y - cy) <= step_y / 2
                     for cx, cy in centres) for x, y in points)
    return ([to_sky(cx, cy, th) for cx, cy in centres],
            100.0 * inside / len(points))


def sequence(fields, sky, sky_every, t_obj, t_sky):
    """The list of exposures: fields, with sky frames in between."""
    steps = []
    for k, (n, e) in enumerate(fields):
        steps.append({"kind": "field", "num": k + 1, "n": n, "e": e,
                      "exp": t_obj})
        last = k == len(fields) - 1
        if sky_every and ((k + 1) % sky_every == 0 or last):
            steps.append({"kind": "sky", "n": sky[0], "e": sky[1],
                          "exp": t_sky})
    return steps


# small (dither) offsets

DITHER_LABELS = ("c ---&gt; a", "a ---&gt; b", "b ---&gt; c", "b ---&gt; d")


def dither_offsets(pa_deg, mode):
    """(N, E) offsets in arcsec for c->a, a->b, b->c, b->d.

    The same formulas as koala_offsets.py.
    """
    s = SCALES[mode]
    pa = math.radians(pa_deg)
    pa60 = math.radians(60 - pa_deg)
    ab = (s * math.sin(pa), -s * math.cos(pa))
    bc = (-s * (math.sin(pa60) + math.sin(pa)),
          s * (math.cos(pa) - math.cos(pa60)))
    ca = (s * math.sin(pa60), s * math.cos(pa60))
    bd = (-s * R3 * math.cos(pa), -s * R3 * math.sin(pa))
    return [ca, ab, bc, bd]


def dither_positions(pa_deg, mode):
    """Pointings a, b, c, d as (N, E) in arcsec relative to a."""
    ca, ab, bc, bd = dither_offsets(pa_deg, mode)
    b = ab
    return {"a": (0.0, 0.0), "b": b,
            "c": (b[0] + bc[0], b[1] + bc[1]),
            "d": (b[0] + bd[0], b[1] + bd[1])}


# coordinates

def parse_angle(text, hours):
    """Decimal degrees, or sexagesimal ('00 34 02.8' = hours if hours)."""
    parts = [x for x in re.split(r"[:\shmsd\u00b0'\"]+", text.strip()) if x]
    if not parts:
        raise ValueError("empty")
    if len(parts) == 1:
        return float(parts[0])  # one number = decimal degrees
    sign = -1 if parts[0].startswith("-") else 1
    vals = [abs(float(x)) for x in parts[:3]] + [0.0, 0.0]
    value = vals[0] + vals[1] / 60 + vals[2] / 3600
    return sign * value * (15 if hours else 1)


def sexagesimal(deg, hours):
    value = (deg % 360) / 15 if hours else deg
    sign = "" if hours else ("-" if value < 0 else "+")
    places = 2 if hours else 1
    total = round(abs(value) * 3600, places)
    a = int(total // 3600)
    b = int((total - a * 3600) // 60)
    c = total - a * 3600 - b * 60
    return "%s%02d %02d %0*.*f" % (sign, a, b, places + 3, places, c)


def coord_text(ra, dec, degrees):
    """RA and Dec as text, in degrees or as hh mm ss / dd mm ss."""
    if degrees:
        return "%.5f" % (ra % 360), "%+.5f" % dec
    return sexagesimal(ra, True), sexagesimal(dec, False)


def offset_to_radec(centre, n, e):
    return (centre[0] + e / 3600 / math.cos(math.radians(centre[1])),
            centre[1] + n / 3600)


def radec_to_offset(centre, ra, dec):
    """Offset (N, E) in arcsec of a position from the centre."""
    dra = (ra - centre[0] + 180) % 360 - 180
    return ((dec - centre[1]) * 3600,
            dra * math.cos(math.radians(centre[1])) * 3600)


_resolved = {}


def resolve(name):
    """Look a name up: built-in list first, then the CDS Sesame service."""
    key = re.sub(r"\s+", "", name).upper()
    if key in KNOWN:
        return KNOWN[key]
    if key not in _resolved:
        try:
            with urlopen(SESAME + quote(name.strip()), timeout=8) as reply:
                text = reply.read().decode("utf-8", "replace")
        except OSError:
            raise ValueError("Could not look up '%s' (no internet?). "
                             "Enter RA and Dec by hand." % name)
        found = re.search(r"^%J\s+([-+\d.]+)\s+([-+\d.]+)", text, re.M)
        if not found:
            raise ValueError("'%s' was not found. Enter RA and Dec by hand."
                             % name)
        _resolved[key] = (float(found.group(1)), float(found.group(2)))
    return _resolved[key]


# drawing

SIZE = 680  # pixels


def fmt(x):
    x = round(x, decimals)
    return "%.*f" % (decimals, x if x != 0 else 0.0)  # avoid "-0.00"


def image_urls(centre, extent, image, size=SIZE):
    """Sky image address, and a second one to try if the first fails."""
    legacy = LEGACY_SERVICE + "?" + urlencode({
        "ra": "%.6f" % centre[0], "dec": "%.6f" % centre[1],
        "layer": LEGACY_LAYER, "pixscale": "%.4f" % (2 * extent / size),
        "size": size})
    dss = HIPS_SERVICE + "?" + urlencode({
        "hips": HIPS_SURVEY, "width": size, "height": size,
        "projection": "TAN", "coordsys": "icrs",
        "fov": "%.6f" % (2 * extent / 3600),
        "ra": "%.6f" % centre[0], "dec": "%.6f" % centre[1], "format": "jpg"})
    return (dss, legacy) if image == "dss" else (legacy, dss)


def sky_image(centre, extent, image, size):
    first, second = image_urls(centre, extent, image, size)
    return ('<image href="%s" data-alt="%s" x="0" y="0" width="%d" '
            'height="%d" onerror="if(this.dataset.alt){this.setAttribute('
            '\'href\',this.dataset.alt);delete this.dataset.alt}'
            'else{this.remove()}"/>'
            % (html.escape(first, quote=True),
               html.escape(second, quote=True), size, size))


def hexagon_corners(p, th, k):
    """Corners of one spaxel in pixels (North up, East left)."""
    radius = p / R3
    out = []
    for i in range(6):
        a = math.radians(30 + 60 * i)
        n, e = to_sky(radius * math.cos(a), radius * math.sin(a), th)
        out.append((-e * k, -n * k))
    return out


def hexagon_path(mode, th, k):
    """All spaxels of one field as an SVG path, relative to its centre."""
    p = SCALES[mode]
    corners = hexagon_corners(p, th, k)

    def pixel(x, y):
        n, e = to_sky(x, y, th)
        return -e * k, -n * k  # North up, East left

    sides = "".join("l%.2f %.2f" % (corners[i + 1][0] - corners[i][0],
                                    corners[i + 1][1] - corners[i][1])
                    for i in range(5)) + "z"
    parts = []
    for r in range(ROWS):
        for c in range(COLS):
            x = (c - (COLS - 1) / 2) * p + (r % 2) * p / 2 - p / 4
            y = (r - (ROWS - 1) / 2) * p * R3 / 2
            px, py = pixel(x, y)
            parts.append("M%.2f %.2f%s" % (px + corners[0][0],
                                           py + corners[0][1], sides))
    return "".join(parts)


def compass(size_marks=True):
    return ('<g stroke="#fff" fill="#fff" font-size="12">'
            '<line x1="46" y1="46" x2="46" y2="16"/>'
            '<line x1="46" y1="46" x2="16" y2="46"/>'
            '<text x="42" y="13" stroke="none">N</text>'
            '<text x="4" y="50" stroke="none">E</text></g>')


def scale_bar(size, k, extent, choices, unit="&quot;"):
    bar = min(choices, key=lambda v: abs(v - extent / 3))
    x2, y = size - 16, size - 16
    return ('<line x1="%.1f" y1="%d" x2="%d" y2="%d" stroke="#fff" '
            'stroke-width="2"/><text x="%.1f" y="%d" font-size="12" '
            'fill="#fff" text-anchor="end">%g%s</text>'
            % (x2 - bar * k, y, x2, y, x2, y - 6, bar, unit))


def centre_marks(xy, outline):
    """Dotted outline of a field on the galaxy centre, and a cross."""
    pts = " ".join("%.1f,%.1f" % xy(dn, de) for dn, de in outline)
    cx, cy = xy(0, 0)
    return ('<polygon points="%s" fill="none" stroke="#fff" '
            'stroke-width="1.6" stroke-dasharray="2 4"/>'
            '<path d="M%.1f %.1fh14M%.1f %.1fv14" stroke="#ff4d4d" '
            'stroke-width="1.8"/>'
            '<text x="%.1f" y="%.1f" font-size="11" fill="#fff" '
            'stroke="#000" stroke-width="3" paint-order="stroke">centre</text>'
            % (pts, cx - 7, cy, cx, cy - 7, cx + 6, cy - 6))


def sky_svg(fields, pa_deg, mode, gal, centre, image):
    th = frame_angle(pa_deg)
    outline = field_outline(pa_deg, mode)
    centred = has_centre_field(fields)

    reach = [gal["a"]] + [math.hypot(dn, de) for dn, de in outline]
    for n, e in fields:
        reach += [math.hypot(n + dn, e + de) for dn, de in outline]
    extent = max(reach) * 1.08  # arcsec from the centre to the map edge
    k = SIZE / 2 / extent

    def xy(n, e):
        return SIZE / 2 - e * k, SIZE / 2 - n * k

    parts = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" '
             'style="max-width:100%%;height:auto">' % ((SIZE,) * 4),
             '<rect width="100%" height="100%" fill="#0b0b12"/>']
    if centre:
        parts.append(sky_image(centre, extent, image, SIZE))

    parts.append('<defs><path id="spaxels" d="%s"/>'
                 '<marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" '
                 'markerWidth="5" markerHeight="5" orient="auto">'
                 '<path d="M0,0 L10,5 L0,10 z" fill="#fff"/></marker></defs>'
                 % hexagon_path(mode, th, k))
    colours = ("#35d0ff", "#ffd23f")  # only to tell neighbours apart
    for idx, (n, e) in enumerate(fields):
        x, y = xy(n, e)
        colour = colours[idx % 2]
        parts.append('<use href="#spaxels" x="%.2f" y="%.2f" fill="%s" '
                     'fill-opacity="0.13" stroke="%s" stroke-opacity="0.75" '
                     'stroke-width="0.35"/>' % (x, y, colour, colour))

    pts = []
    for i in range(120):
        t = math.radians(i * 3)
        pts.append("%.1f,%.1f" % xy(*galaxy_point(math.cos(t),
                                                  math.sin(t), gal)))
    parts.append('<polygon points="%s" fill="none" stroke="#ff4d4d" '
                 'stroke-width="1.6" stroke-dasharray="6 4"/>' % " ".join(pts))

    if not centred:
        parts.append(centre_marks(xy, outline))

    # the path of the telescope from field to field
    for (n1, e1), (n2, e2) in zip(fields, fields[1:]):
        x1, y1 = xy(n1, e1)
        x2, y2 = xy(n2, e2)
        d = math.hypot(x2 - x1, y2 - y1) or 1
        t = min(0.45, 16 / d)  # leave room for the numbers
        parts.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                     'stroke="#fff" stroke-width="1.4" stroke-opacity="0.9" '
                     'marker-end="url(#arr)"/>'
                     % (x1 + (x2 - x1) * t, y1 + (y2 - y1) * t,
                        x2 - (x2 - x1) * t, y2 - (y2 - y1) * t))

    label_offsets = not centred and len(fields) <= MAX_MAP_LABELS
    for idx, (n, e) in enumerate(fields):
        x, y = xy(n, e)
        parts.append('<text x="%.1f" y="%.1f" font-size="15" '
                     'font-weight="bold" text-anchor="middle" fill="#fff" '
                     'stroke="#000" stroke-width="3" paint-order="stroke">'
                     '%d</text>' % (x, y + 5, idx + 1))
        if label_offsets:
            parts.append('<text x="%.1f" y="%.1f" font-size="11" '
                         'text-anchor="middle" fill="#fff" stroke="#000" '
                         'stroke-width="3" paint-order="stroke">'
                         '%s N, %s E</text>' % (x, y + 20, fmt(n), fmt(e)))

    parts.append(compass())
    parts.append(scale_bar(SIZE, k, extent,
                           (1, 2, 5, 10, 20, 30, 60, 120, 300)))
    parts.append("</svg>")
    return "".join(parts)


WIDE_SIZE = 520  # pixels, second (wide) view


def viewer_link(centre, n, e):
    """Address of the interactive Legacy Survey viewer at an offset."""
    dec = centre[1] + n / 3600
    ra = centre[0] + e / 3600 / math.cos(math.radians(centre[1]))
    return LEGACY_VIEWER + "?" + urlencode({
        "ra": "%.4f" % ra, "dec": "%.4f" % dec, "layer": LEGACY_LAYER,
        "zoom": 13})


def wide_svg(fields, pa_deg, mode, gal, centre, sky, image):
    """Wider sky view: outline of every field and of the sky frame."""
    outline = field_outline(pa_deg, mode)
    frames = [(n, e, "#35d0ff", "") for n, e in fields]
    if sky:
        frames.append((sky[0], sky[1], "#ffd23f", "sky"))

    reach = [gal["a"]]
    for n, e, _, _ in frames:
        for dn, de in outline:
            reach += [abs(n + dn), abs(e + de)]
    extent = max(reach) * 1.12
    size = WIDE_SIZE
    k = size / 2 / extent

    def xy(n, e):
        return size / 2 - e * k, size / 2 - n * k

    parts = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" '
             'style="max-width:100%%;height:auto">' % ((size,) * 4),
             '<rect width="100%" height="100%" fill="#0b0b12"/>']
    if centre:
        parts.append(sky_image(centre, extent, image, size))
    for n, e, colour, label in frames:
        pts = " ".join("%.1f,%.1f" % xy(n + dn, e + de) for dn, de in outline)
        parts.append('<polygon points="%s" fill="none" stroke="%s" '
                     'stroke-width="1.2"%s/>'
                     % (pts, colour,
                        ' stroke-dasharray="5 3"' if label else ""))
        if label:
            x, y = xy(n, e)
            parts.append('<text x="%.1f" y="%.1f" font-size="13" '
                         'font-weight="bold" text-anchor="middle" '
                         'fill="%s" stroke="#000" stroke-width="3" '
                         'paint-order="stroke">%s</text>'
                         % (x, y + 4, colour, label))
    pts = []
    for i in range(120):
        t = math.radians(i * 3)
        pts.append("%.1f,%.1f" % xy(*galaxy_point(math.cos(t),
                                                  math.sin(t), gal)))
    parts.append('<polygon points="%s" fill="none" stroke="#ff4d4d" '
                 'stroke-width="1.2" stroke-dasharray="6 4"/>' % " ".join(pts))
    if not has_centre_field(fields):
        parts.append(centre_marks(xy, outline))

    parts.append(compass())
    parts.append(scale_bar(size, k, extent,
                           (10, 20, 30, 60, 120, 300, 600)))
    parts.append("</svg>")
    return "".join(parts)


DITHER_SIZE = 300  # pixels, small-offset view


def dither_svg(pa_deg, mode):
    """Close-up of the spaxels with the dither pointings a, b, c, d."""
    p = SCALES[mode]
    th = frame_angle(pa_deg)
    pos = dither_positions(pa_deg, mode)
    cen_n = sum(n for n, _ in pos.values()) / 4
    cen_e = sum(e for _, e in pos.values()) / 4
    size = DITHER_SIZE
    extent = 2.6 * p
    k = size / 2 / extent

    def xy(n, e):
        return size / 2 - (e - cen_e) * k, size / 2 - (n - cen_n) * k

    parts = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" '
             'style="max-width:100%%;height:auto">' % ((size,) * 4),
             '<rect width="100%" height="100%" fill="#0b0b12"/>',
             '<defs><marker id="darr" viewBox="0 0 10 10" refX="9" '
             'refY="5" markerWidth="6" markerHeight="6" orient="auto">'
             '<path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/>'
             '</marker></defs>']

    # the spaxels at pointing a (a spaxel centre at a)
    corners = hexagon_corners(p, th, k)
    for r in range(-6, 7):
        for c in range(-6, 7):
            n, e = to_sky(c * p + (r % 2) * p / 2, r * p * R3 / 2, th)
            if math.hypot(n - cen_n, e - cen_e) > 1.6 * extent:
                continue
            x, y = xy(n, e)
            parts.append('<polygon points="%s" fill="#35d0ff" '
                         'fill-opacity="0.10" stroke="#35d0ff" '
                         'stroke-opacity="0.6" stroke-width="1"/>'
                         % " ".join("%.1f,%.1f" % (x + dx, y + dy)
                                    for dx, dy in corners))

    def arrow(p1, p2, colour, dash=""):
        x1, y1 = xy(*pos[p1])
        x2, y2 = xy(*pos[p2])
        d = math.hypot(x2 - x1, y2 - y1) or 1
        t = min(0.4, 9 / d)
        return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                'stroke-width="2"%s marker-end="url(#darr)"/>'
                % (x1 + (x2 - x1) * t, y1 + (y2 - y1) * t,
                   x2 - (x2 - x1) * t, y2 - (y2 - y1) * t, colour, dash))

    parts.append(arrow("a", "b", "#fff"))
    parts.append(arrow("b", "c", "#fff"))
    parts.append(arrow("c", "a", "#fff"))
    parts.append(arrow("b", "d", "#ffd23f", ' stroke-dasharray="5 3"'))
    for name, colour in (("a", "#ff4d4d"), ("b", "#fff"), ("c", "#fff"),
                         ("d", "#ffd23f")):
        x, y = xy(*pos[name])
        parts.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>'
                     '<text x="%.1f" y="%.1f" font-size="15" '
                     'font-weight="bold" fill="%s" stroke="#000" '
                     'stroke-width="3" paint-order="stroke">%s</text>'
                     % (x, y, colour, x + 7, y - 7, colour, name))
    parts.append(compass())
    parts.append(scale_bar(size, k, extent, (p,),
                           '&quot; = 1 spaxel'))
    parts.append("</svg>")
    return "".join(parts)


DITHER_COLOURS = {"a": "#ff4d4d", "b": "#ffffff", "c": "#35d0ff", "d": "#ffd23f"}


def dither_map_svg(field_no, fields, pa_deg, mode, gal, centre, image):
    """One mosaic field at the small-offset pointings a, b, c, d, drawn
    over the galaxy at a scale where the small offsets can be seen."""
    th = frame_angle(pa_deg)
    outline = field_outline(pa_deg, mode)
    fn, fe = fields[field_no - 1]
    pos = {name: (fn + n, fe + e)
           for name, (n, e) in dither_positions(pa_deg, mode).items()}
    cen_n = sum(n for n, _ in pos.values()) / 4
    cen_e = sum(e for _, e in pos.values()) / 4

    reach = []
    for n, e in pos.values():
        for dn, de in outline:
            reach += [abs(n + dn - cen_n), abs(e + de - cen_e)]
    extent = max(reach) * 1.15
    size = SIZE
    k = size / 2 / extent

    def xy(n, e):
        return size / 2 - (e - cen_e) * k, size / 2 - (n - cen_n) * k

    parts = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" '
             'style="max-width:100%%;height:auto">' % ((size,) * 4),
             '<rect width="100%" height="100%" fill="#0b0b12"/>']
    if centre:
        view = offset_to_radec(centre, cen_n, cen_e)
        parts.append(sky_image(view, extent, image, size))
    parts.append('<defs><path id="dspaxels" d="%s"/>'
                 '<marker id="darr2" viewBox="0 0 10 10" refX="9" refY="5" '
                 'markerWidth="5" markerHeight="5" orient="auto">'
                 '<path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/>'
                 '</marker></defs>' % hexagon_path(mode, th, k))

    # neighbouring mosaic fields, faint, for context
    for idx, (n, e) in enumerate(fields):
        if idx == field_no - 1:
            continue
        pts = " ".join("%.1f,%.1f" % xy(n + dn, e + de) for dn, de in outline)
        parts.append('<polygon points="%s" fill="none" stroke="#35d0ff" '
                     'stroke-opacity="0.45" stroke-width="1"/>' % pts)
        x, y = xy(n, e)
        if -40 < x < size + 40 and -40 < y < size + 40:
            parts.append('<text x="%.1f" y="%.1f" font-size="14" '
                         'text-anchor="middle" fill="#9fe6ff" stroke="#000" '
                         'stroke-width="3" paint-order="stroke">%d</text>'
                         % (x, y + 5, idx + 1))

    # the spaxels while at a
    x, y = xy(*pos["a"])
    parts.append('<use href="#dspaxels" x="%.2f" y="%.2f" fill="#ff4d4d" '
                 'fill-opacity="0.10" stroke="#ff4d4d" stroke-opacity="0.45" '
                 'stroke-width="0.4"/>' % (x, y))

    # field of view at a, b, c, d
    for name in "abcd":
        n, e = pos[name]
        pts = " ".join("%.1f,%.1f" % xy(n + dn, e + de) for dn, de in outline)
        parts.append('<polygon points="%s" fill="none" stroke="%s" '
                     'stroke-width="1.5"%s/>'
                     % (pts, DITHER_COLOURS[name],
                        ' stroke-dasharray="6 3"' if name == "d" else ""))

    # galaxy outline and centre
    pts = []
    for i in range(120):
        t = math.radians(i * 3)
        pts.append("%.1f,%.1f" % xy(*galaxy_point(math.cos(t),
                                                  math.sin(t), gal)))
    parts.append('<polygon points="%s" fill="none" stroke="#ff4d4d" '
                 'stroke-width="1.6" stroke-dasharray="6 4"/>' % " ".join(pts))
    cx, cy = xy(0, 0)
    parts.append('<path d="M%.1f %.1fh14M%.1f %.1fv14" stroke="#ff4d4d" '
                 'stroke-width="1.8"/>' % (cx - 7, cy, cx, cy - 7))

    # order of the small offsets
    def arrow(p1, p2, colour, dash=""):
        x1, y1 = xy(*pos[p1])
        x2, y2 = xy(*pos[p2])
        d = math.hypot(x2 - x1, y2 - y1) or 1
        t = min(0.3, 5 / d)
        return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                'stroke-width="1.6"%s marker-end="url(#darr2)"/>'
                % (x1 + (x2 - x1) * t, y1 + (y2 - y1) * t,
                   x2 - (x2 - x1) * t, y2 - (y2 - y1) * t, colour, dash))

    parts.append(arrow("a", "b", "#fff"))
    parts.append(arrow("b", "c", "#fff"))
    parts.append(arrow("c", "a", "#fff"))
    parts.append(arrow("b", "d", "#ffd23f", ' stroke-dasharray="4 2"'))
    cx, cy = xy(cen_n, cen_e)
    for name in "abcd":
        x, y = xy(*pos[name])
        dx, dy = x - cx, y - cy
        dl = math.hypot(dx, dy) or 1
        parts.append('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>'
                     '<text x="%.1f" y="%.1f" font-size="15" '
                     'font-weight="bold" text-anchor="middle" fill="%s" '
                     'stroke="#000" stroke-width="3" paint-order="stroke">'
                     '%s</text>'
                     % (x, y, DITHER_COLOURS[name], x + dx / dl * 16,
                        y + dy / dl * 16 + 5, DITHER_COLOURS[name], name))

    parts.append(compass())
    parts.append(scale_bar(size, k, extent, (1, 2, 5, 10, 20, 30)))
    parts.append("</svg>")
    return "".join(parts)


def dither_positions_table(field_no, result, mode, form):
    """Positions of a, b, c, d of one field, from the galaxy centre."""
    fields, centre = result["fields"], result["centre"]
    degrees = form.get("coords") == "deg"
    fn, fe = fields[field_no - 1]
    pos = dither_positions(result["pa"], mode)
    rows = ['<tr><th></th><th colspan="2" class="sep">offset from centre'
            '</th>%s</tr><tr><th>at</th><th class="sep">N</th><th>E</th>%s'
            '</tr>' % (('<th colspan="2" class="sep">position (J2000)</th>',
                        '<th class="sep">RA</th><th>Dec</th>') if centre
                       else ("", ""))]
    for name in "abcd":
        n, e = fn + pos[name][0], fe + pos[name][1]
        radec = ""
        if centre:
            radec = '<td class="sep">%s</td><td>%s</td>' % coord_text(
                *offset_to_radec(centre, n, e), degrees)
        rows.append('<tr><td><span class="dot" style="background:%s"></span>'
                    '<b>%s</b></td><td class="sep">%s</td><td>%s</td>%s</tr>'
                    % (DITHER_COLOURS[name], name, fmt(n), fmt(e), radec))
    return '<table class="data">%s</table>' % "".join(rows)


def dither_map_card(field_no, result, mode, form):
    return """<div class="card">
<h2>Small offsets over the galaxy &middot; field %d</h2>
%s
<p class="note">Field of view of field %d at the small-offset pointings
<span class="key" style="color:#d33">a</span> (red),
<b>b</b> (white), <span class="key" style="color:#1a9fd6">c</span> (blue)
and <span class="key" style="color:#c99a00">d</span> (yellow dashes). Red
hexagons = spaxels while at a; faint blue = the neighbouring mosaic fields,
with their numbers. Arrows give the order a &rarr; b &rarr; c &rarr; a; the
dashed arrow b &rarr; d is the alternative step from b. Red dashes = galaxy
outline, red cross = galaxy centre. North up, East left.</p></div>""" % (
        field_no,
        dither_map_svg(field_no, result["fields"], result["pa"], mode,
                       result["gal"], result["centre"],
                       form.get("image", "legacy")),
        field_no)


def dither_card(pa_deg, mode, positions=""):
    rows = "".join('<tr><td class="l">Offset %s :</td><td class="sep">%s'
                   '</td><td>%s</td></tr>' % (label, fmt(n), fmt(e))
                   for label, (n, e) in zip(DITHER_LABELS,
                                            dither_offsets(pa_deg, mode)))
    if positions:
        positions = ('<h3>Positions of a, b, c, d (arcsec from the galaxy '
                     'centre)</h3>' + positions)
    return """<div class="card">
<h2>Small offsets &middot; PA %g &middot; %s (%s&quot;/spaxel)</h2>
<div class="pair">
 <div><h3>Offsets (arcsec)</h3>
  <table class="data"><tr><th></th><th class="sep">N</th><th>E</th></tr>%s
  </table>%s</div>
 <div>%s
  <p class="note" style="max-width:%dpx">Arcsec, the same as
  koala_offsets.py. Dots: where the telescope points at a, b, c and d;
  hexagons: the spaxels while at a. The dashed yellow arrow b &rarr; d is
  the alternative step from b. North up, East left.</p></div>
</div></div>""" % (pa_deg, mode, SCALES[mode], rows, positions,
                   dither_svg(pa_deg, mode), DITHER_SIZE)


# ------------------------------------------------------------------- page

STYLE = """
* { box-sizing: border-box; }
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       margin: 0; color: #1d2330; background: #e9ecf1; font-size: 14px; }
header { background: #141a2b; color: #fff; padding: 14px 24px; }
header h1 { font-size: 18px; margin: 0; font-weight: 600; }
header p { margin: 3px 0 0 0; color: #aab3c8; font-size: 12px; }
main { padding: 20px 24px 40px 24px; }
.row { display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-start;
       margin-bottom: 20px; }
.card { background: #fff; border-radius: 10px; padding: 16px 18px;
        box-shadow: 0 1px 3px rgba(20,26,43,.12); max-width: 100%; }
.card.wide { flex: 1 1 100%; }
.row.two { display: grid; grid-template-columns: minmax(0, 1fr) auto; }
.row.two > .card { min-width: 0; overflow-x: auto; }
.row.two > .card:only-child { grid-column: 1 / -1; justify-self: start; }
@media (max-width: 1150px) {
  .row.two { grid-template-columns: minmax(0, 1fr); }
}
h2 { font-size: 14px; margin: 0 0 10px 0; font-weight: 600; color: #141a2b; }
h3 { font-size: 12px; margin: 12px 0 2px 0; font-weight: 600; color: #555;
     text-transform: uppercase; letter-spacing: .04em; }
.pair { display: flex; flex-wrap: wrap; gap: 18px; align-items: flex-start; }
table { border-collapse: separate; border-spacing: 10px 5px; }
td, th { white-space: nowrap; text-align: center; font-weight: normal; }
td.l { text-align: left; }
th { color: #667; font-size: 12px; }
td.sep, th.sep { border-left: 1px solid #d5d9e2; padding-left: 12px; }
table.data { border-spacing: 0; margin-left: -6px; }
table.data td, table.data th { padding: 4px 8px; }
table.data tr:nth-child(even) td { background: #f5f7fa; }
tr.sky td { color: #778; font-style: italic; }
tr.ref td { color: #a00; }
table.form { border-spacing: 8px 5px; margin: -5px -8px; }
table.form td.l { color: #333; }
input, select { width: 7.5em; font-size: 14px; padding: 4px 6px;
                border: 1px solid #c3c9d4; border-radius: 5px;
                background: #fff; }
input:focus, select:focus { outline: 2px solid #35a3ff; border-color: #35a3ff; }
button { font-size: 14px; padding: 5px 16px; border-radius: 5px;
         border: 1px solid #1d4ed8; background: #2563eb; color: #fff;
         cursor: pointer; font-weight: 600; }
button:hover { background: #1d4ed8; }
button.alt { background: #fff; color: #1d4ed8; }
button.alt:hover { background: #eef3ff; }
.y { background: #ffe14d; padding: 3px 8px; border-radius: 4px;
     font-weight: 600; }
.err { color: #b00020; font-weight: bold; margin: 0; }
.warn { color: #a00; margin: 10px 0 2px 0; }
ul.offs { margin: 4px 0 0 0; columns: 4; }
.note { font-size: 12px; color: #556; margin: 8px 0 0 0; max-width: 680px;
        line-height: 1.45; }
.opt { color: #778; font-size: 12px; }
.summary p { margin: 0 0 4px 0; line-height: 1.5; }
.dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%;
       border: 1px solid #333; margin-right: 6px; vertical-align: -1px; }
.key { font-weight: 700; }
.empty { width: 520px; height: 200px; display: flex; align-items: center;
         justify-content: center; color: #889; }
"""


def duration(seconds):
    minutes = int(round(seconds / 60))
    return "%d h %02d min" % divmod(minutes, 60)


def sequence_table(steps, centre, degrees, reference_row):
    head = ('<tr><th></th><th></th><th colspan="2" class="sep">offset from '
            'centre</th><th colspan="2" class="sep">move from previous</th>'
            '%s<th class="sep"></th></tr>'
            '<tr><th>#</th><th>target</th><th class="sep">N</th><th>E</th>'
            '<th class="sep">N</th><th>E</th>%s<th class="sep">exp (s)</th>'
            '</tr>'
            % (('<th colspan="2" class="sep">position (J2000%s)</th>'
                % (", degrees" if degrees else ""),
                '<th class="sep">RA</th><th>Dec</th>') if centre
               else ("", "")))
    rows = [head]
    if reference_row:
        radec = ""
        if centre:
            radec = '<td class="sep">%s</td><td>%s</td>' % coord_text(
                *centre, degrees)
        rows.append('<tr class="ref"><td>&ndash;</td><td class="l">Galaxy '
                    'centre (start here, no exposure)</td>'
                    '<td class="sep">%s</td><td>%s</td>'
                    '<td class="sep" colspan="2">start here</td>%s'
                    '<td class="sep">&ndash;</td></tr>'
                    % (fmt(0), fmt(0), radec))
    prev_n = prev_e = 0.0  # the telescope starts on the galaxy centre
    for idx, s in enumerate(steps):
        n, e = s["n"], s["e"]
        if s["kind"] == "sky":
            target = "Sky"
        else:
            target = "Field %d" % s["num"]
            if at_centre(n, e):
                target += " (centre)"
        if idx == 0 and at_centre(n, e):
            move = '<td class="sep" colspan="2">start here</td>'
        else:
            move = '<td class="sep">%s</td><td>%s</td>' % (
                fmt(n - prev_n), fmt(e - prev_e))
        radec = ""
        if centre:
            radec = '<td class="sep">%s</td><td>%s</td>' % coord_text(
                *offset_to_radec(centre, n, e), degrees)
        rows.append('<tr class="%s"><td>%d</td><td class="l">%s</td>'
                    '<td class="sep">%s</td><td>%s</td>%s%s'
                    '<td class="sep">%g</td></tr>'
                    % (s["kind"], idx + 1, target, fmt(n), fmt(e), move,
                       radec, s["exp"]))
        prev_n, prev_e = n, e
    return "<table>%s</table>" % "\n".join(rows)


def result_parts(form, mode, result, show_dither):
    """The pieces of the results page, each an HTML card."""
    fields, centre, steps, sky = (result[k] for k in
                                  ("fields", "centre", "steps", "sky"))
    p = SCALES[mode]
    pa = result["pa"]
    long_axis, short_axis = field_size(mode)
    ov_long, ov_short = result["overlap"]
    degrees = form.get("coords") == "deg"
    image = form.get("image", "legacy")
    centred = has_centre_field(fields)

    where = ""
    if centre:
        where = " Centre: RA %s, Dec %s = %s, %s deg (J2000)." % (
            coord_text(*centre, False) + coord_text(*centre, True))
    n_sky = sum(s["kind"] == "sky" for s in steps)
    if sky:
        where += "<br>Sky frame: %s N, %s E of the centre (arcsec)" % (
            fmt(sky[0]), fmt(sky[1]))
        if centre:
            at = offset_to_radec(centre, *sky)
            where += ", at RA %s, Dec %s = %s, %s deg" % (
                coord_text(*at, False) + coord_text(*at, True))
        where += "."
    links = ""
    if centre:
        links = (' <a href="%s" target="_blank" rel="noopener">Open the '
                 'Legacy Survey viewer on the galaxy</a>'
                 % html.escape(viewer_link(centre, 0, 0), quote=True))
        if sky:
            links += (' &middot; <a href="%s" target="_blank" '
                      'rel="noopener">on the sky position</a>'
                      % html.escape(viewer_link(centre, *sky), quote=True))
    travel = prev_n = prev_e = 0.0
    for n, e in fields:
        travel += math.hypot(n - prev_n, e - prev_e)
        prev_n, prev_e = n, e

    overlap = "no overlap"
    if ov_long or ov_short:
        overlap = ("overlap %d spaxels (%.1f&quot;) along the long axis, "
                   "%d rows (%.1f&quot;) along the short axis"
                   % (ov_long, ov_long * p, ov_short, ov_short * p * R3 / 2))

    summary = """<div class="card wide summary">
<p>%s mode, %s arcsec/spaxel, field %.1f&quot; x %.1f&quot;, rotator PA %g
   (long axis along sky PA %g&deg;): <b>%s</b>, %s, covering %.0f%% of the
   galaxy outline.%s</p>
<p>%d field exposures + %d sky frames = <b>%s</b> of exposure time
   (slews and read-out not included). Field-to-field travel: %.0f&quot;.</p>
""" % (mode, p, long_axis, short_axis, pa, (pa + PA_ZERO) % 180,
       result["what"], overlap, result["percent"], where,
       len(fields), n_sky, duration(sum(s["exp"] for s in steps)), travel)
    if not centred:
        summary += (
            '<p class="warn"><b>No field is on the galaxy centre.</b> Point '
            'at the centre first (dotted outline on the maps), then offset '
            'to each field. Offsets from the centre (arcsec):</p>'
            '<ul class="offs">%s</ul>'
            % "".join("<li>Field %d: %s N, %s E</li>" % (i + 1, fmt(n),
                                                          fmt(e))
                      for i, (n, e) in enumerate(fields)))
    summary += "</div>"

    centre_note = ""
    if not centred:
        centre_note = (" Dotted white outline and red cross = a field on "
                       "the galaxy centre (not observed); each field is "
                       "labelled with its offset from it.")

    parts = {"summary": summary}
    parts["wide"] = """<div class="card">
<h2>Sky view</h2>%s
<p class="note" style="max-width:520px">Blue = outline of each field, yellow
dashes = the sky frame, red dashes = galaxy outline. Check that the sky
frame falls on empty sky.%s</p></div>""" % (
        wide_svg(fields, pa, mode, result["gal"], centre, sky, image), links)
    parts["mosaic"] = """<div class="card">
<h2>Observing order</h2>%s
<p class="note">Each field is drawn as its %d hexagonal spaxels; blue and
yellow only tell neighbouring fields apart. Numbers and arrows give the
observing order. Red dashes = galaxy outline. North up, East left.%s%s</p>
</div>""" % (sky_svg(fields, pa, mode, result["gal"], centre, image),
             COLS * ROWS, centre_note,
             "" if centre else " Give a name or RA/Dec to see the sky image.")
    parts["table"] = """<div class="card">
<h2>Sequence: fields and sky</h2>%s
<p class="note">Offsets in arcsec. Negative N = South, negative E = West.
&quot;Offset from centre&quot; is measured from the galaxy centre, where the
telescope starts; &quot;move from previous&quot; is the same thing as a step
from the last position. Use whichever your offset command expects, not
both.</p></div>""" % sequence_table(steps, centre, degrees, not centred)
    if show_dither:
        parts["dmap"] = dither_map_card(result["dfield"], result, mode, form)
        parts["dclose"] = dither_card(
            pa, mode, dither_positions_table(result["dfield"], result, mode,
                                             form))
    return parts


def render(form, mode, result, error):
    def val(name):
        return html.escape(form.get(name, ""), quote=True)

    def options(name, choices):
        return "".join('<option value="%s"%s>%s</option>'
                       % (v, " selected" if form.get(name) == v else "", text)
                       for v, text in choices)

    show_dither = form.get("dither", "show") == "show"
    parts = {}
    if result:
        parts = result_parts(form, mode, result, show_dither)
    elif show_dither:  # the small offsets need only the PA
        try:
            pa = float(form.get("pa", "").strip())
            if math.isfinite(pa):
                parts["dclose"] = dither_card(pa, mode)
        except ValueError:
            pass

    error_html = ('<p class="err">%s</p>' % html.escape(error)
                  if error else "")
    button = ('<button type="submit" name="mode" value="%s"%s>%s</button>')
    form_card = """<div class="card">
<h2>Observation</h2>
<form method="get" action="/">
<table class="form">
<tr><td class="l"><span class="y">Position Angle (PA):</span></td>
    <td><input name="pa" value="%s" autocomplete="off" autofocus></td>
    <td colspan="2" class="l">%s %s</td></tr>
<tr><td></td><td colspan="3" class="l opt">rotator PA; 0 = long axis East-West</td></tr>
<tr><td class="l">Galaxy name:</td>
    <td><input name="name" value="%s" autocomplete="off"></td>
    <td colspan="2" class="l opt">used to find RA/Dec for the sky image</td></tr>
<tr><td class="l">or RA, Dec (J2000):</td>
    <td><input name="ra" value="%s" autocomplete="off"></td>
    <td><input name="dec" value="%s" autocomplete="off"></td>
    <td class="l opt">degrees, or "00 34 02.8"</td></tr>
<tr><td class="l">Galaxy major x minor axis (arcmin):</td>
    <td><input name="gmaj" value="%s" autocomplete="off"></td>
    <td><input name="gmin" value="%s" autocomplete="off"></td>
    <td class="l opt">full diameters</td></tr>
<tr><td class="l">Galaxy PA (deg, N through E):</td>
    <td><input name="gpa" value="%s" autocomplete="off"></td></tr>
<tr><td class="l">Fields: columns x rows:</td>
    <td><input name="ncols" value="%s" autocomplete="off"></td>
    <td><input name="nrows" value="%s" autocomplete="off"></td>
    <td class="l opt">blank = cover the whole galaxy</td></tr>
<tr><td class="l">Order:</td>
    <td><select name="order">%s</select></td>
    <td colspan="2" class="l opt">spiral starts on the galaxy centre</td></tr>
<tr><td class="l">Overlap (spaxels): long, short axis:</td>
    <td><input name="ovlong" value="%s" autocomplete="off"></td>
    <td><input name="ovshort" value="%s" autocomplete="off"></td>
    <td class="l opt">0 = none; up to %d of %d, %d of %d</td></tr>
<tr><td class="l">Exposure: field, sky (s):</td>
    <td><input name="texp" value="%s" autocomplete="off"></td>
    <td><input name="tsky" value="%s" autocomplete="off"></td></tr>
<tr><td class="l">Sky frame after every ... fields:</td>
    <td><input name="skyevery" value="%s" autocomplete="off"></td>
    <td colspan="2" class="l opt">0 = no sky frames</td></tr>
<tr><td class="l">Sky offset from centre N, E (arcmin):</td>
    <td><input name="skyn" value="%s" autocomplete="off"></td>
    <td><input name="skye" value="%s" autocomplete="off"></td>
    <td class="l opt">-4, 0 = 4' South</td></tr>
<tr><td class="l">or sky RA, Dec (J2000):</td>
    <td><input name="skyra" value="%s" autocomplete="off"></td>
    <td><input name="skydec" value="%s" autocomplete="off"></td>
    <td class="l opt">replaces the offset</td></tr>
<tr><td class="l">Small offsets (a, b, c, d), on field:</td>
    <td><select name="dither">%s</select></td>
    <td><input name="dfield" value="%s" autocomplete="off"></td>
    <td class="l opt">as in koala_offsets.py</td></tr>
<tr><td class="l">Sky image:</td>
    <td><select name="image">%s</select></td>
    <td colspan="2" class="l opt">Legacy does not cover the whole sky</td></tr>
<tr><td class="l">Show coordinates as:</td>
    <td><select name="coords">%s</select></td></tr>
</table>
</form>
%s</div>""" % (val("pa"),
               button % ("WIDE", "" if mode == "WIDE" else ' class="alt"',
                         "WIDE"),
               button % ("NARROW", "" if mode == "NARROW" else ' class="alt"',
                         "NARROW"),
               val("name"), val("ra"), val("dec"),
               val("gmaj"), val("gmin"), val("gpa"), val("ncols"),
               val("nrows"),
               options("order", (("spiral", "spiral"),
                                 ("snake", "row by row"))),
               val("ovlong"), val("ovshort"),
               MAX_OVERLAP_LONG, COLS, MAX_OVERLAP_SHORT, ROWS,
               val("texp"), val("tsky"),
               val("skyevery"), val("skyn"), val("skye"),
               val("skyra"), val("skydec"),
               options("dither", (("show", "show"), ("hide", "hide"))),
               val("dfield"),
               options("image", (("legacy", "Legacy"), ("dss", "DSS2"))),
               options("coords", (("sex", "hh mm ss"),
                                  ("deg", "degrees"))),
               error_html)

    def row(*keys):
        cards = [parts[k] for k in keys if parts.get(k)]
        return ('<div class="row two">%s</div>' % "".join(cards)
                if cards else "")

    top = form_card + parts.get("wide", "")
    body = ('<div class="row">%s</div>' % top
            + ('<div class="row">%s</div>' % parts["summary"]
               if "summary" in parts else "")
            + row("mosaic", "table")
            + row("dmap", "dclose"))

    return """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>KOALA: mosaic observing sequence</title>
<style>%s</style></head>
<body>
<header><h1>KOALA: mosaic observing sequence for a given position angle</h1>
<p>Mosaic fields, observing order, sky frames and the small offsets
a, b, c, d &middot; offsets in arcsec, North up, East left</p></header>
<main>%s</main>
</body></html>""" % (STYLE, body)


def number(form, name, label):
    try:
        value = float(form.get(name, "").strip())
    except ValueError:
        raise ValueError("%s has to be a number!" % label)
    if not math.isfinite(value):
        raise ValueError("%s has to be a number!" % label)
    return value


def whole(form, name, label, low, high):
    value = number(form, name, label)
    if value != int(value) or not low <= value <= high:
        raise ValueError("%s has to be a whole number from %d to %d!"
                         % (label, low, high))
    return int(value)


def position(ra_text, dec_text, label):
    """RA, Dec in degrees from text in degrees or hh mm ss / dd mm ss."""
    try:
        ra = parse_angle(ra_text, True)
        dec = parse_angle(dec_text, False)
    except (ValueError, IndexError):
        raise ValueError("%s RA and Dec both have to be given, in degrees "
                         "or as hh mm ss and dd mm ss!" % label)
    if not (0 <= ra <= 360 and -90 <= dec <= 90):
        raise ValueError("%s RA has to be 0-360 deg (0-24 h) and Dec "
                         "-90 to +90!" % label)
    return ra, dec


def compute(form, mode):
    pa = number(form, "pa", "Position Angle")
    gmaj = number(form, "gmaj", "Galaxy major axis")
    gmin = number(form, "gmin", "Galaxy minor axis")
    gpa = number(form, "gpa", "Galaxy PA")
    if gmaj <= 0 or gmin <= 0:
        raise ValueError("Galaxy axes have to be greater than zero!")
    if max(gmaj, gmin) > 30:
        raise ValueError("Galaxy axes are in arcmin; 30 is the limit here.")
    overlap = (whole(form, "ovlong", "Overlap along the long axis",
                     0, MAX_OVERLAP_LONG),
               whole(form, "ovshort", "Overlap along the short axis",
                     0, MAX_OVERLAP_SHORT))
    gmaj, gmin = max(gmaj, gmin), min(gmaj, gmin)
    gal = {"a": gmaj * 30, "b": gmin * 30, "pa": gpa}  # half axes, arcsec

    block = None
    if form.get("ncols", "").strip() or form.get("nrows", "").strip():
        block = (whole(form, "ncols", "Columns", 1, MAX_BLOCK),
                 whole(form, "nrows", "Rows", 1, MAX_BLOCK))

    t_obj = number(form, "texp", "Field exposure")
    t_sky = number(form, "tsky", "Sky exposure")
    if t_obj < 0 or t_sky < 0:
        raise ValueError("Exposure times cannot be negative!")
    sky_every = whole(form, "skyevery", "Sky frame frequency", 0, 99)
    centre = None
    ra_text, dec_text = form.get("ra", "").strip(), form.get("dec", "").strip()
    if ra_text or dec_text:
        centre = position(ra_text, dec_text, "Galaxy")
    elif form.get("name", "").strip():
        centre = resolve(form["name"])

    sky_ra = form.get("skyra", "").strip()
    sky_dec = form.get("skydec", "").strip()
    if sky_ra or sky_dec:
        if not centre:
            raise ValueError("Sky RA/Dec need the galaxy position too: "
                             "give a galaxy name or its RA and Dec.")
        sky = radec_to_offset(centre, *position(sky_ra, sky_dec, "Sky"))
        if math.hypot(*sky) > 3600:
            raise ValueError("The sky position is more than 1 degree from "
                             "the galaxy. Check the sky RA and Dec.")
    else:
        sky = (number(form, "skyn", "Sky offset N") * 60,
               number(form, "skye", "Sky offset E") * 60)

    order = "snake" if form.get("order") == "snake" else "spiral"
    fields, percent = mosaic(pa, mode, gal, overlap, block, order)
    dfield = 1
    if form.get("dfield", "").strip():
        dfield = min(whole(form, "dfield", "Field for the small offsets",
                           1, 999), len(fields))
    what = "%d fields" % len(fields)
    if block:
        what = "%d x %d = %s" % (block[0], block[1], what)
    return {"fields": fields, "percent": percent, "what": what,
            "steps": sequence(fields, sky, sky_every, t_obj, t_sky),
            "sky": sky if sky_every else None, "overlap": overlap,
            "dfield": dfield,
            "pa": pa, "gal": gal, "centre": centre}


DEFAULTS = {"pa": "0", "name": "NGC 151", "ra": "", "dec": "",
            "gmaj": "3.7", "gmin": "1.7", "gpa": "75",
            "ovlong": "0", "ovshort": "0",
            "ncols": "3", "nrows": "3", "order": "spiral",
            "texp": "1200", "tsky": "300", "skyevery": "2",
            "skyn": "-4", "skye": "0", "skyra": "", "skydec": "",
            "dither": "show", "dfield": "1", "image": "legacy", "coords": "sex"}


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        url = urlparse(self.path)
        if url.path != "/":
            self.send_error(404)
            return

        query = parse_qs(url.query, keep_blank_values=True)
        form = dict(DEFAULTS)
        form.update({k: v[0] for k, v in query.items()})
        # Pressing Enter in a box (no button) behaves like WIDE.
        mode = form.get("mode", "WIDE")
        if mode not in SCALES:
            mode = "WIDE"

        result, error = None, ""
        try:
            result = compute(form, mode)
        except ValueError as exc:
            error = str(exc)

        page = render(form, mode, result, error).encode("utf-8")
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
    print("KOALA mosaic calculator running at " + address)
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
