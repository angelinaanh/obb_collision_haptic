"""
Builds the SVG charts of the README from the CSV written by obb_benchmark.

Usage:  python docs/make_charts.py <benchmark.csv> [output_dir]

No dependency beyond the standard library: the charts are written as SVG.
"""

import csv
import math
import sys
from pathlib import Path

# reference palette (categorical slots 1-3, validated all-pairs on the light surface)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'

W, H = 780, 420
LEFT, RIGHT, TOP, BOTTOM = 72, 190, 92, 56


def fmt_number(v):
    if v >= 1:
        return f"{v:,.0f}"
    return f"{v:g}"


def fmt_triangles(v):
    return {1e3: "1K", 1e4: "10K", 1e5: "100K", 1e6: "1M"}.get(v, f"{v:,.0f}")


def log_chart(path, title, subtitle, y_label, x_values, series):
    """series: list of (name, [y values]) plotted on log-log axes."""
    xs = [math.log10(x) for x in x_values]
    ys_all = [y for _, values in series for y in values if y > 0]
    y_lo = math.floor(math.log10(min(ys_all)))
    y_hi = math.ceil(math.log10(max(ys_all)))
    if y_hi == y_lo:
        y_hi += 1
    x_lo, x_hi = math.floor(min(xs)), math.ceil(max(xs))

    plot_w, plot_h = W - LEFT - RIGHT, H - TOP - BOTTOM

    def px(lx):
        return LEFT + (lx - x_lo) / (x_hi - x_lo) * plot_w

    def py(ly):
        return TOP + (y_hi - ly) / (y_hi - y_lo) * plot_h

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
           f'font-family=\'{FONT}\' role="img" aria-label="{title}">',
           f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
           f'<text x="{LEFT}" y="30" font-size="17" font-weight="600" fill="{INK}">{title}</text>',
           f'<text x="{LEFT}" y="52" font-size="13" fill="{INK_2}">{subtitle}</text>']

    # legend (always present for two or more series)
    lx = LEFT
    for i, (name, _) in enumerate(series):
        out.append(f'<rect x="{lx}" y="66" width="14" height="4" rx="2" fill="{SERIES[i]}"/>')
        out.append(f'<text x="{lx + 20}" y="72" font-size="12" fill="{INK_2}">{name}</text>')
        lx += 28 + 6.8 * len(name)

    # gridlines and y ticks (powers of ten)
    for k in range(y_lo, y_hi + 1):
        y = py(k)
        out.append(f'<line x1="{LEFT}" x2="{LEFT + plot_w}" y1="{y:.1f}" y2="{y:.1f}" stroke="{GRID}" stroke-width="1"/>')
        out.append(f'<text x="{LEFT - 8}" y="{y + 4:.1f}" font-size="11" fill="{MUTED}" text-anchor="end">{fmt_number(10.0 ** k)}</text>')
    out.append(f'<text transform="translate(18,{TOP + plot_h / 2:.1f}) rotate(-90)" font-size="12" fill="{MUTED}" text-anchor="middle">{y_label}</text>')

    # x axis
    base = TOP + plot_h
    out.append(f'<line x1="{LEFT}" x2="{LEFT + plot_w}" y1="{base}" y2="{base}" stroke="{AXIS}" stroke-width="1"/>')
    for x in x_values:
        out.append(f'<text x="{px(math.log10(x)):.1f}" y="{base + 20}" font-size="11" fill="{MUTED}" text-anchor="middle">{fmt_triangles(x)}</text>')
    out.append(f'<text x="{LEFT + plot_w / 2:.1f}" y="{H - 12}" font-size="12" fill="{MUTED}" text-anchor="middle">Number of object triangles (log scale)</text>')

    # lines, end markers with a surface ring
    ends = []
    for i, (name, values) in enumerate(series):
        points = [(px(math.log10(x)), py(math.log10(y))) for x, y in zip(x_values, values) if y > 0]
        d = " ".join(f"{'M' if j == 0 else 'L'}{p[0]:.1f},{p[1]:.1f}" for j, p in enumerate(points))
        out.append(f'<path d="{d}" fill="none" stroke="{SERIES[i]}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
        for p in points:
            out.append(f'<circle cx="{p[0]:.1f}" cy="{p[1]:.1f}" r="4" fill="{SERIES[i]}" stroke="{SURFACE}" stroke-width="2"/>')
        ends.append((points[-1][1], values[-1], i))

    # direct end labels, only when they do not collide
    ends.sort()
    if all(b[0] - a[0] >= 16 for a, b in zip(ends, ends[1:])):
        for y, value, i in ends:
            label = f"{value:.2f}" if value < 100 else f"{value:,.0f}"
            out.append(f'<text x="{LEFT + plot_w + 10}" y="{y + 4:.1f}" font-size="12" fill="{INK_2}">'
                       f'{label} - {series[i][0]}</text>')

    out.append("</svg>")
    Path(path).write_text("\n".join(out), encoding="utf-8")
    print("wrote", path)


def main():
    rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8")))
    out_dir = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images")
    out_dir.mkdir(parents=True, exist_ok=True)

    # x = nominal sizes, so that the axis ticks sit on the data points
    nominal = [1e3, 1e4, 1e5, 1e6][:len(rows)]
    col = lambda name: [float(r[name]) for r in rows]

    log_chart(out_dir / "build_time.svg",
              "Tree construction time",
              "OBB tree (SAH, tight OBB per node) vs CHAI3D AABB tree, in seconds",
              "Seconds (log scale)", nominal,
              [("OBB tree (SAH)", col("build_sah_s")),
               ("CHAI3D AABB", col("build_aabb_s"))])

    log_chart(out_dir / "segment_query.svg",
              "Segment query time (haptic tool)",
              "Nearest collision of one segment, called through the CHAI3D API, in microseconds",
              "Microseconds (log scale)", nominal,
              [("OBB tree (SAH)", col("seg_obb_chai3d_us")),
               ("CHAI3D AABB", col("seg_aabb_chai3d_us")),
               ("Brute force", col("seg_brute_us"))])

    log_chart(out_dir / "hand_object_query.svg",
              "Hand vs object query time",
              "All intersecting triangle pairs, hand of 3,132 triangles; reference = uniform grid + exact tests, in microseconds",
              "Microseconds (log scale)", nominal,
              [("OBB tree (SAH)", col("hand_all_us")),
               ("OBB tree (median)", col("hand_all_median_us")),
               ("Reference (grid)", col("hand_reference_us"))])


if __name__ == "__main__":
    main()
