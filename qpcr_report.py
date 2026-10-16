#!/usr/bin/env python3
"""
Visual report for qpcr_design.py.

Runs the design and writes ONE self-contained HTML file: a map of where every
assay sits on the target, a Tm plot, and the full table. No dependencies beyond
the standard library -- open the HTML in a browser or in VS Code's Live Preview.

Usage:
    python qpcr_report.py --fasta pathogen.fasta --out report.html
"""

from __future__ import annotations

import argparse
import html
import sys
from dataclasses import asdict
from pathlib import Path
from typing import List, Sequence

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))      # find qpcr_design.py beside this file

import qpcr_design as qd


# --- palette: categorical slots 1-3, validated for all-pairs CVD separation ---
CSS = """
:root {
  color-scheme: light;
  --surface-0: #f6f6f4;
  --surface-1: #fcfcfb;
  --border:    #dcdcd6;
  --text-1:    #0b0b0b;
  --text-2:    #52514e;
  --text-3:    #86857e;
  --forward:   #2a78d6;
  --probe:     #eb6834;
  --reverse:   #1baf7a;
  --track:     #e7e7e1;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --surface-0: #121211;
    --surface-1: #1a1a19;
    --border:    #34342f;
    --text-1:    #ffffff;
    --text-2:    #c3c2b7;
    --text-3:    #8f8e85;
    --forward:   #3987e5;
    --probe:     #d95926;
    --reverse:   #199e70;
    --track:     #2b2b27;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 32px 24px 64px;
  background: var(--surface-0); color: var(--text-1);
  font: 15px/1.55 ui-sans-serif, -apple-system, "Segoe UI", Roboto, sans-serif;
}
.wrap { max-width: 1040px; margin: 0 auto; }
h1 { font-size: 22px; margin: 0 0 4px; letter-spacing: -0.01em; }
h2 { font-size: 15px; margin: 40px 0 12px; letter-spacing: 0.02em;
     text-transform: uppercase; color: var(--text-2); font-weight: 600; }
.sub { color: var(--text-2); margin: 0 0 28px; font-size: 14px; }
.tiles { display: flex; flex-wrap: wrap; gap: 12px; }
.tile { flex: 1 1 140px; background: var(--surface-1); border: 1px solid var(--border);
        border-radius: 10px; padding: 14px 16px; }
.tile .n { font-size: 26px; font-weight: 600; letter-spacing: -0.02em;
           font-variant-numeric: tabular-nums; }
.tile .l { font-size: 12px; color: var(--text-2); margin-top: 2px; }
.panel { background: var(--surface-1); border: 1px solid var(--border);
         border-radius: 10px; padding: 18px 16px; overflow-x: auto; }
.legend { display: flex; gap: 18px; flex-wrap: wrap; margin: 0 0 14px;
          font-size: 13px; color: var(--text-2); }
.key { display: inline-flex; align-items: center; gap: 7px; }
.sw { width: 11px; height: 11px; border-radius: 3px; display: inline-block; }
svg { display: block; max-width: 100%; height: auto; }
svg text { font: 11px ui-sans-serif, -apple-system, "Segoe UI", sans-serif;
           fill: var(--text-2); }
svg text.row { fill: var(--text-1); font-weight: 500; }
svg text.val { font-variant-numeric: tabular-nums; }
.grid { stroke: var(--border); stroke-width: 1; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid var(--border);
         white-space: nowrap; }
th { color: var(--text-2); font-weight: 600; font-size: 12px;
     text-transform: uppercase; letter-spacing: 0.03em; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
code { font: 12.5px ui-monospace, "SF Mono", Menlo, Consolas, monospace;
       letter-spacing: 0.02em; }
.f { color: var(--forward); } .p { color: var(--probe); } .r { color: var(--reverse); }
.note { color: var(--text-3); font-size: 12.5px; margin-top: 10px; }
"""


def esc(text: object) -> str:
    return html.escape(str(text))


def target_map(assays: Sequence[qd.Assay], seq_len: int) -> str:
    """Where each assay sits on the target, drawn to scale."""
    left, right, top, row_h = 118, 26, 34, 26
    width = 940
    height = top + row_h * len(assays) + 26
    plot_w = width - left - right

    def x_of(pos: int) -> float:
        return left + plot_w * (pos - 1) / max(seq_len - 1, 1)

    out = [f'<svg viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="Assay positions along the target">']

    # scale
    step = 10 ** (len(str(seq_len // 6)) - 1) * max(1, round((seq_len / 6) / 10 ** (len(str(seq_len // 6)) - 1)))
    tick = 0
    while tick <= seq_len:
        x = x_of(max(tick, 1))
        out.append(f'<line class="grid" x1="{x:.1f}" y1="{top - 10}" x2="{x:.1f}" y2="{height - 22}"/>')
        out.append(f'<text class="val" x="{x:.1f}" y="{height - 8}" text-anchor="middle">{tick or 1}</text>')
        tick += step
    out.append(f'<text x="{left}" y="{top - 18}">position on target (bp)</text>')

    for i, a in enumerate(assays):
        y = top + row_h * i
        cy = y + row_h / 2
        amp_x0, amp_x1 = x_of(a.forward_start), x_of(a.forward_start + a.amplicon_size - 1)
        out.append(f'<text class="row" x="0" y="{cy + 4:.1f}">assay {i + 1}</text>')
        out.append(f'<text class="val" x="74" y="{cy + 4:.1f}" text-anchor="end">'
                   f'{a.amplicon_size} bp</text>')
        out.append(f'<rect x="{amp_x0:.1f}" y="{cy - 3:.1f}" width="{max(amp_x1 - amp_x0, 1):.1f}" '
                   f'height="6" rx="3" fill="var(--track)"><title>amplicon '
                   f'{a.forward_start}-{a.forward_start + a.amplicon_size - 1}</title></rect>')
        blocks = (
            ("var(--forward)", a.forward_start, len(a.forward_primer),
             f"forward {a.forward_primer} Tm {a.forward_tm} C"),
            ("var(--probe)", a.probe_start, len(a.probe),
             f"probe {a.probe} Tm {a.probe_tm} C"),
            ("var(--reverse)", a.reverse_start, len(a.reverse_primer),
             f"reverse {a.reverse_primer} Tm {a.reverse_tm} C"),
        )
        for fill, start, length, tip in blocks:
            bx0, bx1 = x_of(start), x_of(start + length - 1)
            out.append(f'<rect x="{bx0:.1f}" y="{cy - 7:.1f}" '
                       f'width="{max(bx1 - bx0, 2.5):.1f}" height="14" rx="4" fill="{fill}" '
                       f'stroke="var(--surface-1)" stroke-width="2">'
                       f'<title>{esc(tip)}</title></rect>')
    out.append("</svg>")
    return "".join(out)


def tm_plot(assays: Sequence[qd.Assay]) -> str:
    """Primer and probe Tm per assay, on one shared axis."""
    left, right, top, row_h = 118, 26, 34, 26
    width = 940
    height = top + row_h * len(assays) + 26
    plot_w = width - left - right

    values = [v for a in assays for v in (a.forward_tm, a.reverse_tm, a.probe_tm)]
    lo, hi = min(values) - 2, max(values) + 2

    def x_of(tm: float) -> float:
        return left + plot_w * (tm - lo) / (hi - lo)

    out = [f'<svg viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="Melting temperature of each oligo">']
    tick = int(lo) + (1 if lo % 1 else 0)
    while tick <= hi:
        if tick % 2 == 0:
            x = x_of(tick)
            out.append(f'<line class="grid" x1="{x:.1f}" y1="{top - 10}" x2="{x:.1f}" y2="{height - 22}"/>')
            out.append(f'<text class="val" x="{x:.1f}" y="{height - 8}" text-anchor="middle">{tick}</text>')
        tick += 1
    out.append(f'<text x="{left}" y="{top - 18}">melting temperature (&#176;C)</text>')

    for i, a in enumerate(assays):
        y = top + row_h * i
        cy = y + row_h / 2
        out.append(f'<text class="row" x="0" y="{cy + 4:.1f}">assay {i + 1}</text>')
        xs = [x_of(a.forward_tm), x_of(a.reverse_tm), x_of(a.probe_tm)]
        out.append(f'<line x1="{min(xs):.1f}" y1="{cy:.1f}" x2="{max(xs):.1f}" y2="{cy:.1f}" '
                   f'stroke="var(--track)" stroke-width="2"/>')
        for fill, tm, label in (("var(--forward)", a.forward_tm, "forward"),
                                ("var(--reverse)", a.reverse_tm, "reverse"),
                                ("var(--probe)", a.probe_tm, "probe")):
            out.append(f'<circle cx="{x_of(tm):.1f}" cy="{cy:.1f}" r="5.5" fill="{fill}" '
                       f'stroke="var(--surface-1)" stroke-width="2">'
                       f'<title>{label} Tm {tm} &#176;C</title></circle>')
        if i == 0:   # one direct label, not a number on every point
            out.append(f'<text class="val" x="{x_of(a.probe_tm) + 10:.1f}" y="{cy + 4:.1f}">'
                       f'probe {a.probe_tm}&#176;C</text>')
    out.append("</svg>")
    return "".join(out)


def table(assays: Sequence[qd.Assay]) -> str:
    head = ("#", "Forward primer", "Tm", "Reverse primer", "Tm",
            "Probe", "Tm", "Amplicon", "Position", "Score")
    rows = ["<tr>" + "".join(f"<th>{esc(h)}</th>" for h in head) + "</tr>"]
    for i, a in enumerate(assays, start=1):
        rows.append(
            f"<tr><td class='num'>{i}</td>"
            f"<td><code class='f'>{esc(a.forward_primer)}</code></td>"
            f"<td class='num'>{a.forward_tm}</td>"
            f"<td><code class='r'>{esc(a.reverse_primer)}</code></td>"
            f"<td class='num'>{a.reverse_tm}</td>"
            f"<td><code class='p'>{esc(a.probe)}</code></td>"
            f"<td class='num'>{a.probe_tm}</td>"
            f"<td class='num'>{a.amplicon_size} bp</td>"
            f"<td class='num'>{a.forward_start}&#8211;"
            f"{a.forward_start + a.amplicon_size - 1}</td>"
            f"<td class='num'>{a.score}</td></tr>")
    return "<table>" + "".join(rows) + "</table>"


def build_html(record: qd.FastaRecord, assays: List[qd.Assay],
               n_primers: int, n_assays: int) -> str:
    legend = (
        '<div class="legend">'
        '<span class="key"><span class="sw" style="background:var(--forward)"></span>forward primer</span>'
        '<span class="key"><span class="sw" style="background:var(--probe)"></span>probe</span>'
        '<span class="key"><span class="sw" style="background:var(--reverse)"></span>reverse primer</span>'
        '</div>')
    tiles = (
        f'<div class="tiles">'
        f'<div class="tile"><div class="n">{len(record.sequence):,}</div>'
        f'<div class="l">target length (bp)</div></div>'
        f'<div class="tile"><div class="n">{n_primers:,}</div>'
        f'<div class="l">primers passing filters</div></div>'
        f'<div class="tile"><div class="n">{n_assays:,}</div>'
        f'<div class="l">valid primer/probe sets</div></div>'
        f'<div class="tile"><div class="n">{len(assays)}</div>'
        f'<div class="l">reported after de-duplication</div></div>'
        f'</div>')

    body = (f'<h1>qPCR assay design &#8212; {esc(record.identifier)}</h1>'
            f'<p class="sub">{esc(record.description) or "target sequence"}</p>'
            f'{tiles}'
            f'<h2>Where the assays sit on the target</h2>'
            f'<div class="panel">{legend}{target_map(assays, len(record.sequence))}'
            f'<p class="note">Hover any block for its sequence and Tm. '
            f'Grey bar = amplicon span.</p></div>'
            f'<h2>Melting temperatures</h2>'
            f'<div class="panel">{legend}{tm_plot(assays)}'
            f'<p class="note">Primers are matched to each other; the probe sits '
            f'{qd.PROBE_TM_ABOVE_PRIMER[0]:.0f}&#8211;{qd.PROBE_TM_ABOVE_PRIMER[1]:.0f}&#176;C above them '
            f'so it binds before extension reaches it.</p></div>'
            f'<h2>All reported assays</h2>'
            f'<div class="panel">{table(assays)}'
            f'<p class="note">Specificity is not checked here &#8212; BLAST every oligo '
            f'against the host and near-neighbour genomes before ordering.</p></div>')

    return ("<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width, initial-scale=1'>"
            f"<title>qPCR assay design &#8212; {esc(record.identifier)}</title>"
            f"<style>{CSS}</style></head><body><div class='wrap'>{body}</div></body></html>")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Visual report for qpcr_design.py")
    parser.add_argument("--fasta", default=str(HERE / "pathogen.fasta"))
    parser.add_argument("--out", default=str(HERE / "qpcr_report.html"))
    parser.add_argument("--top", type=int, default=8)
    args = parser.parse_args(argv)

    records = qd.read_fasta(args.fasta)
    record = records[0]
    if len(records) > 1:
        print(f"note: {len(records)} records found, reporting on {record.identifier}")

    primers = qd.find_candidate_primers(record.sequence)
    assays = qd.design_assays(record)
    selected = qd.select_non_redundant(assays, limit=args.top)
    if not selected:
        print("No assay met the criteria; nothing to plot.")
        return 1

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(build_html(record, selected, len(primers), len(assays)))
    print(f"{record.identifier}: {len(primers)} primers, {len(assays)} assays, "
          f"{len(selected)} reported")
    print(f"Wrote {args.out} -- open it in a browser.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
