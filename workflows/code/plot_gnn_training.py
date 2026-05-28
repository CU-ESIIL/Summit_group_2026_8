#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


COLORS = {
    "blue": "#2563eb",
    "orange": "#f97316",
    "green": "#16a34a",
    "red": "#dc2626",
    "purple": "#7c3aed",
    "teal": "#0d9488",
    "gray": "#64748b",
    "black": "#111827",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SVG diagnostic plots from a GNN training run directory."
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Training run directory containing metrics.csv and test_metrics.json.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Directory for plots. Defaults to <run-dir>/plots.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional human-readable title prefix for plot headings.",
    )
    return parser.parse_args()


def read_metrics(path: Path) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: Dict[str, float] = {}
            for key, value in raw.items():
                if value is None or value == "":
                    row[key] = math.nan
                    continue
                try:
                    row[key] = float(value)
                except ValueError:
                    row[key] = math.nan
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def read_json(path: Path) -> Dict:
    if not path.exists():
        return {}
    with path.open() as handle:
        return json.load(handle)


def finite(values: Sequence[float]) -> List[float]:
    return [v for v in values if isinstance(v, (int, float)) and math.isfinite(v)]


def column(rows: Sequence[Dict[str, float]], name: str) -> List[float]:
    return [row.get(name, math.nan) for row in rows]


def best_epoch(rows: Sequence[Dict[str, float]]) -> int:
    vals = [(row.get("val_loss", math.nan), int(row.get("epoch", i + 1))) for i, row in enumerate(rows)]
    vals = [(loss, epoch) for loss, epoch in vals if math.isfinite(loss)]
    if not vals:
        return int(rows[-1].get("epoch", len(rows)))
    return min(vals, key=lambda item: item[0])[1]


def fmt(value: Optional[float], digits: int = 3) -> str:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return "NA"
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.{digits}f}"


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def write_text(
    parts: List[str],
    x: float,
    y: float,
    text: object,
    size: int = 12,
    fill: str = COLORS["black"],
    anchor: str = "start",
    weight: str = "400",
) -> None:
    parts.append(
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" '
        f'font-family="Arial, Helvetica, sans-serif" fill="{fill}" '
        f'text-anchor="{anchor}" font-weight="{weight}">{esc(text)}</text>'
    )


def nice_ticks(min_value: float, max_value: float, n_ticks: int = 5) -> List[float]:
    if not math.isfinite(min_value) or not math.isfinite(max_value):
        return [0, 1]
    if min_value == max_value:
        pad = abs(min_value) * 0.1 or 1.0
        min_value -= pad
        max_value += pad
    raw_step = (max_value - min_value) / max(1, n_ticks - 1)
    power = 10 ** math.floor(math.log10(raw_step))
    step = min([1, 2, 2.5, 5, 10], key=lambda m: abs(raw_step - m * power)) * power
    start = math.floor(min_value / step) * step
    end = math.ceil(max_value / step) * step
    ticks = []
    current = start
    while current <= end + step * 0.5:
        ticks.append(current)
        current += step
    return ticks


def line_path(
    xs: Sequence[float],
    ys: Sequence[float],
    x_to_px,
    y_to_px,
) -> str:
    commands: List[str] = []
    drawing = False
    for x, y in zip(xs, ys):
        if not (math.isfinite(x) and math.isfinite(y)):
            drawing = False
            continue
        px = x_to_px(x)
        py = y_to_px(y)
        if drawing:
            commands.append(f"L {px:.2f} {py:.2f}")
        else:
            commands.append(f"M {px:.2f} {py:.2f}")
            drawing = True
    return " ".join(commands)


def draw_line_plot(
    out_path: Path,
    title: str,
    ylabel: str,
    series: Sequence[Tuple[str, Sequence[float], Sequence[float], str, str]],
    best: Optional[int] = None,
    y_limits: Optional[Tuple[float, float]] = None,
    note: Optional[str] = None,
    width: int = 980,
    height: int = 620,
) -> None:
    margin_left = 82
    margin_right = 38
    margin_top = 82
    margin_bottom = 92
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    all_x = finite([x for _, xs, _, _, _ in series for x in xs])
    all_y = finite([y for _, _, ys, _, _ in series for y in ys])
    if not all_x or not all_y:
        raise ValueError(f"No finite data available for {out_path}")

    x_min, x_max = min(all_x), max(all_x)
    if x_min == x_max:
        x_min -= 1
        x_max += 1

    if y_limits:
        y_min, y_max = y_limits
    else:
        y_min, y_max = min(all_y), max(all_y)
        pad = (y_max - y_min) * 0.08 or 0.1
        y_min -= pad
        y_max += pad
        if y_min > 0 and min(all_y) >= 0:
            y_min = 0

    x_ticks = nice_ticks(x_min, x_max)
    y_ticks = nice_ticks(y_min, y_max)
    y_min, y_max = min(y_ticks), max(y_ticks)

    def x_to_px(x: float) -> float:
        return margin_left + (x - x_min) / (x_max - x_min) * plot_w

    def y_to_px(y: float) -> float:
        return margin_top + (y_max - y) / (y_max - y_min) * plot_h

    parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    write_text(parts, margin_left, 34, title, size=24, weight="700")
    if note:
        write_text(parts, margin_left, 58, note, size=13, fill=COLORS["gray"])

    parts.append(
        f'<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>'
    )

    for tick in y_ticks:
        y = y_to_px(tick)
        parts.append(
            f'<line x1="{margin_left}" x2="{margin_left + plot_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e2e8f0"/>'
        )
        write_text(parts, margin_left - 10, y + 4, fmt(tick), size=11, fill=COLORS["gray"], anchor="end")

    for tick in x_ticks:
        if tick < x_min - 1e-9 or tick > x_max + 1e-9:
            continue
        x = x_to_px(tick)
        parts.append(
            f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{margin_top}" y2="{margin_top + plot_h}" stroke="#eef2f7"/>'
        )
        write_text(parts, x, margin_top + plot_h + 24, fmt(tick, digits=0), size=11, fill=COLORS["gray"], anchor="middle")

    if best is not None and x_min <= best <= x_max:
        bx = x_to_px(float(best))
        parts.append(
            f'<line x1="{bx:.2f}" x2="{bx:.2f}" y1="{margin_top}" y2="{margin_top + plot_h}" '
            'stroke="#111827" stroke-width="1.6" stroke-dasharray="5 5"/>'
        )
        write_text(parts, bx + 7, margin_top + 18, f"best epoch {best}", size=12, fill=COLORS["black"])

    for label, xs, ys, color, dash in series:
        path = line_path(xs, ys, x_to_px, y_to_px)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(
            f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.6" stroke-linejoin="round"{dash_attr}/>'
        )

    write_text(parts, margin_left + plot_w / 2, height - 30, "epoch", size=13, fill=COLORS["black"], anchor="middle")
    write_text(
        parts,
        24,
        margin_top + plot_h / 2,
        ylabel,
        size=13,
        fill=COLORS["black"],
        anchor="middle",
    )
    parts[-1] = parts[-1].replace(
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}"',
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}" transform="rotate(-90 24 {margin_top + plot_h / 2:.2f})"',
    )

    legend_x = margin_left
    legend_y = height - 62
    for i, (label, _, _, color, dash) in enumerate(series):
        x = legend_x + i * 220
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(
            f'<line x1="{x}" x2="{x + 28}" y1="{legend_y}" y2="{legend_y}" '
            f'stroke="{color}" stroke-width="3"{dash_attr}/>'
        )
        write_text(parts, x + 36, legend_y + 4, label, size=12, fill=COLORS["black"])

    parts.append("</svg>")
    out_path.write_text("\n".join(parts))


def draw_bar_plot(
    out_path: Path,
    title: str,
    categories: Sequence[str],
    values_a: Sequence[float],
    label_a: str,
    values_b: Optional[Sequence[float]] = None,
    label_b: Optional[str] = None,
    ylabel: str = "value",
    note: Optional[str] = None,
    width: int = 980,
    height: int = 560,
) -> None:
    margin_left = 86
    margin_right = 34
    margin_top = 86
    margin_bottom = 112
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    all_values = finite(list(values_a) + list(values_b or []))
    if not all_values:
        raise ValueError(f"No finite bar values available for {out_path}")
    y_min = 0
    y_max = max(all_values) * 1.12 or 1.0
    y_ticks = nice_ticks(y_min, y_max)
    y_min, y_max = min(y_ticks), max(y_ticks)

    def y_to_px(y: float) -> float:
        return margin_top + (y_max - y) / (y_max - y_min) * plot_h

    n = len(categories)
    group_w = plot_w / max(1, n)
    bar_w = group_w * (0.28 if values_b is not None else 0.42)

    parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    write_text(parts, margin_left, 34, title, size=24, weight="700")
    if note:
        write_text(parts, margin_left, 58, note, size=13, fill=COLORS["gray"])
    parts.append(
        f'<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="#cbd5e1"/>'
    )

    for tick in y_ticks:
        y = y_to_px(tick)
        parts.append(
            f'<line x1="{margin_left}" x2="{margin_left + plot_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e2e8f0"/>'
        )
        write_text(parts, margin_left - 10, y + 4, fmt(tick), size=11, fill=COLORS["gray"], anchor="end")

    for i, category in enumerate(categories):
        center = margin_left + group_w * (i + 0.5)
        if values_b is None:
            xs = [center - bar_w / 2]
            vals = [values_a[i]]
            colors = [COLORS["blue"]]
        else:
            xs = [center - bar_w - 3, center + 3]
            vals = [values_a[i], values_b[i]]
            colors = [COLORS["blue"], COLORS["orange"]]

        for x, val, color in zip(xs, vals, colors):
            if not math.isfinite(val):
                continue
            y = y_to_px(val)
            bar_h = margin_top + plot_h - y
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" '
                f'fill="{color}" rx="2"/>'
            )
            write_text(parts, x + bar_w / 2, y - 7, fmt(val), size=11, fill=COLORS["black"], anchor="middle")

        label = category
        if len(label) > 16:
            label = label.replace(" ", "\n")
        label_lines = label.split("\n")
        for j, line in enumerate(label_lines):
            write_text(
                parts,
                center,
                margin_top + plot_h + 24 + j * 14,
                line,
                size=11,
                fill=COLORS["gray"],
                anchor="middle",
            )

    write_text(
        parts,
        24,
        margin_top + plot_h / 2,
        ylabel,
        size=13,
        fill=COLORS["black"],
        anchor="middle",
    )
    parts[-1] = parts[-1].replace(
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}"',
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}" transform="rotate(-90 24 {margin_top + plot_h / 2:.2f})"',
    )

    legend_y = height - 34
    parts.append(f'<rect x="{margin_left}" y="{legend_y - 12}" width="14" height="14" fill="{COLORS["blue"]}"/>')
    write_text(parts, margin_left + 22, legend_y, label_a, size=12)
    if values_b is not None and label_b:
        parts.append(
            f'<rect x="{margin_left + 190}" y="{legend_y - 12}" width="14" height="14" fill="{COLORS["orange"]}"/>'
        )
        write_text(parts, margin_left + 212, legend_y, label_b, size=12)

    parts.append("</svg>")
    out_path.write_text("\n".join(parts))


def metric_at_best(rows: Sequence[Dict[str, float]], metric: str, best: int) -> Optional[float]:
    for row in rows:
        if int(row.get("epoch", -1)) == best:
            value = row.get(metric)
            return value if value is not None and math.isfinite(value) else None
    return None


def write_summary(
    out_path: Path,
    run_dir: Path,
    rows: Sequence[Dict[str, float]],
    test_metrics: Dict,
    run_config: Dict,
    best: int,
    plot_files: Sequence[Path],
) -> None:
    args = run_config.get("args", {})
    feature_manifest = run_config.get("feature_manifest", {})
    n_epochs = len(rows)
    val_loss = metric_at_best(rows, "val_loss", best)
    val_auc = metric_at_best(rows, "val_presence_auroc", best)
    val_ap = metric_at_best(rows, "val_presence_average_precision", best)
    val_rmse_present = metric_at_best(rows, "val_log_rmse_present", best)

    lines = [
        "# GNN Training Diagnostic Plots",
        "",
        f"Run directory: `{run_dir}`",
        f"Dataset: `{args.get('dataset_dir', 'unknown')}`",
        f"Kingdom: `{args.get('kingdom', 'unknown')}`",
        f"Model: `{args.get('model', 'unknown')}`",
        f"Edge scope: `{feature_manifest.get('edge_scope', 'unknown')}`",
        f"Epochs recorded: `{n_epochs}`",
        f"Best epoch by validation loss: `{best}`",
        "",
        "## Best Validation Epoch",
        "",
        f"- Validation loss: `{fmt(val_loss)}`",
        f"- Validation presence AUROC: `{fmt(val_auc)}`",
        f"- Validation presence average precision: `{fmt(val_ap)}`",
        f"- Validation log-RMSE on present taxa: `{fmt(val_rmse_present)}`",
        "",
        "## Held-Out Test Metrics",
        "",
        f"- Presence rate: `{fmt(test_metrics.get('presence_rate'))}`",
        f"- Predicted presence rate at threshold 0.5: `{fmt(test_metrics.get('presence_predicted_rate'))}`",
        f"- Presence AUROC: `{fmt(test_metrics.get('presence_auroc'))}`",
        f"- Presence average precision: `{fmt(test_metrics.get('presence_average_precision'))}`",
        f"- GNN log-RMSE on present taxa: `{fmt(test_metrics.get('log_rmse_present'))}`",
        f"- Resident baseline log-RMSE on present taxa: `{fmt(test_metrics.get('resident_log_rmse_present'))}`",
        "",
        "## Sanity Checks For 'Too Good To Be True'",
        "",
        "- Prefer this resident-edge run over any `resident_final` run for strict prediction, because final-derived edges can leak target information.",
        "- AUROC is high, but average precision is the better sparse-target check; compare it to the test presence rate.",
        "- If predicted presence rate is much higher than true presence rate, the ranking can be good while the 0.5 threshold is poorly calibrated.",
        "- Compare present-taxon abundance error to the resident baseline; compare all-taxon abundance error carefully because zeros dominate that metric.",
        "- Next controls: repeat with `edge_scope=none`, try multiple seeds, and compare an edge-aware model such as `gine`.",
        "",
        "## Plot Files",
        "",
    ]
    for plot_file in plot_files:
        lines.append(f"- `{plot_file.name}`")
    out_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    out_dir = Path(args.out_dir) if args.out_dir else run_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_path = run_dir / "metrics.csv"
    test_metrics_path = run_dir / "test_metrics.json"
    run_config_path = run_dir / "run_config.json"

    rows = read_metrics(metrics_path)
    test_metrics = read_json(test_metrics_path)
    run_config = read_json(run_config_path)
    best = best_epoch(rows)
    epochs = column(rows, "epoch")

    title_prefix = args.title or run_dir.name
    plot_files: List[Path] = []

    loss_plot = out_dir / "01_loss_curves.svg"
    draw_line_plot(
        loss_plot,
        f"{title_prefix}: Loss Curves",
        "loss",
        [
            ("train total", epochs, column(rows, "train_loss"), COLORS["blue"], ""),
            ("val total", epochs, column(rows, "val_loss"), COLORS["orange"], ""),
            ("train abundance", epochs, column(rows, "train_abundance_loss"), COLORS["teal"], "6 4"),
            ("val abundance", epochs, column(rows, "val_abundance_loss"), COLORS["purple"], "6 4"),
            ("train presence", epochs, column(rows, "train_presence_loss"), COLORS["green"], "3 4"),
            ("val presence", epochs, column(rows, "val_presence_loss"), COLORS["red"], "3 4"),
        ],
        best=best,
        note="Early stopping saves the checkpoint at the minimum validation loss.",
    )
    plot_files.append(loss_plot)

    rank_plot = out_dir / "02_presence_ranking_metrics.svg"
    draw_line_plot(
        rank_plot,
        f"{title_prefix}: Final-Presence Ranking",
        "score",
        [
            ("validation AUROC", epochs, column(rows, "val_presence_auroc"), COLORS["blue"], ""),
            ("validation average precision", epochs, column(rows, "val_presence_average_precision"), COLORS["orange"], ""),
        ],
        best=best,
        y_limits=(0, 1),
        note="Average precision is usually the sharper check when final-present taxa are rare.",
    )
    plot_files.append(rank_plot)

    abundance_plot = out_dir / "03_abundance_error.svg"
    draw_line_plot(
        abundance_plot,
        f"{title_prefix}: Log-Abundance Error",
        "log-scale error",
        [
            ("RMSE present taxa", epochs, column(rows, "val_log_rmse_present"), COLORS["blue"], ""),
            ("MAE present taxa", epochs, column(rows, "val_log_mae_present"), COLORS["teal"], ""),
            ("RMSE all taxa", epochs, column(rows, "val_log_rmse_all"), COLORS["orange"], "6 4"),
            ("MAE all taxa", epochs, column(rows, "val_log_mae_all"), COLORS["purple"], "6 4"),
        ],
        best=best,
        note="Present-taxon error reflects abundance prediction after establishment/persistence.",
    )
    plot_files.append(abundance_plot)

    rate_plot = out_dir / "04_presence_rate_calibration.svg"
    draw_line_plot(
        rate_plot,
        f"{title_prefix}: Presence Rate At Threshold 0.5",
        "fraction of taxon targets",
        [
            ("true validation presence rate", epochs, column(rows, "val_presence_rate"), COLORS["black"], ""),
            ("predicted validation presence rate", epochs, column(rows, "val_presence_predicted_rate"), COLORS["red"], ""),
        ],
        best=best,
        y_limits=(0, max(0.2, max(finite(column(rows, "val_presence_predicted_rate")) + [0.0]) * 1.15)),
        note="A strong ranker can still need threshold calibration.",
    )
    plot_files.append(rate_plot)

    baseline_plot = out_dir / "05_test_abundance_vs_resident_baseline.svg"
    categories = ["RMSE present", "MAE present", "RMSE all", "MAE all"]
    draw_bar_plot(
        baseline_plot,
        f"{title_prefix}: Held-Out Abundance Error",
        categories,
        [
            test_metrics.get("log_rmse_present", math.nan),
            test_metrics.get("log_mae_present", math.nan),
            test_metrics.get("log_rmse_all", math.nan),
            test_metrics.get("log_mae_all", math.nan),
        ],
        "GNN",
        [
            test_metrics.get("resident_log_rmse_present", math.nan),
            test_metrics.get("resident_log_mae_present", math.nan),
            test_metrics.get("resident_log_rmse_all", math.nan),
            test_metrics.get("resident_log_mae_all", math.nan),
        ],
        "resident baseline",
        ylabel="log-scale error",
        note="The present-taxa comparison asks whether the model improves over final = resident.",
    )
    plot_files.append(baseline_plot)

    test_rank_plot = out_dir / "06_test_presence_summary.svg"
    draw_bar_plot(
        test_rank_plot,
        f"{title_prefix}: Held-Out Presence Summary",
        ["AUROC", "Avg precision", "true rate", "predicted rate"],
        [
            test_metrics.get("presence_auroc", math.nan),
            test_metrics.get("presence_average_precision", math.nan),
            test_metrics.get("presence_rate", math.nan),
            test_metrics.get("presence_predicted_rate", math.nan),
        ],
        "test metric",
        ylabel="score or fraction",
        note="Average precision should be interpreted relative to the true presence rate.",
    )
    plot_files.append(test_rank_plot)

    write_summary(
        out_dir / "README.md",
        run_dir=run_dir,
        rows=rows,
        test_metrics=test_metrics,
        run_config=run_config,
        best=best,
        plot_files=plot_files,
    )

    print(f"Wrote {len(plot_files)} SVG plots to {out_dir}")
    print(f"Best epoch by validation loss: {best}")
    print(f"Summary: {out_dir / 'README.md'}")


if __name__ == "__main__":
    main()
