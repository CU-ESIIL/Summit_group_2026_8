#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


COLORS = [
    "#2563eb",
    "#f97316",
    "#16a34a",
    "#dc2626",
    "#7c3aed",
    "#0d9488",
    "#9333ea",
    "#ca8a04",
]
BLACK = "#111827"
GRAY = "#64748b"
GRID = "#e2e8f0"
PANEL = "#f8fafc"


@dataclass
class Run:
    label: str
    path: Path
    color: str
    metrics: List[Dict[str, float]]
    test_metrics: Dict[str, float]
    config: Dict
    best_epoch: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay GNN training curves from multiple run directories."
    )
    parser.add_argument(
        "--run",
        action="append",
        default=[],
        metavar="LABEL=RUN_DIR",
        help="Run to compare. Can be repeated, for example sage_resident=workflows/output/...",
    )
    parser.add_argument(
        "--out-dir",
        default="workflows/output/gnn_training/model_comparison",
        help="Directory for comparison SVGs and summary CSV/README.",
    )
    parser.add_argument(
        "--title",
        default="Fungi GNN Model Comparison",
        help="Title prefix for plots.",
    )
    return parser.parse_args()


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


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


def parse_run_arg(raw: str) -> Tuple[str, Path]:
    if "=" not in raw:
        raise ValueError(f"--run must look like LABEL=RUN_DIR, got: {raw}")
    label, path = raw.split("=", 1)
    label = label.strip()
    if not label:
        raise ValueError(f"Empty run label in --run {raw!r}")
    return label, Path(path)


def finite(values: Sequence[float]) -> List[float]:
    return [v for v in values if isinstance(v, (int, float)) and math.isfinite(v)]


def best_epoch(rows: Sequence[Dict[str, float]]) -> int:
    candidates = []
    for i, row in enumerate(rows):
        loss = row.get("val_loss", math.nan)
        if math.isfinite(loss):
            candidates.append((loss, int(row.get("epoch", i + 1))))
    if not candidates:
        return int(rows[-1].get("epoch", len(rows)))
    return min(candidates, key=lambda item: item[0])[1]


def value_at_epoch(rows: Sequence[Dict[str, float]], metric: str, epoch: int) -> float:
    for row in rows:
        if int(row.get("epoch", -1)) == epoch:
            return row.get(metric, math.nan)
    return math.nan


def fmt(value: Optional[float], digits: int = 3) -> str:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return "NA"
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.{digits}f}"


def column(rows: Sequence[Dict[str, float]], name: str) -> List[float]:
    return [row.get(name, math.nan) for row in rows]


def nice_ticks(min_value: float, max_value: float, n_ticks: int = 5) -> List[float]:
    if not math.isfinite(min_value) or not math.isfinite(max_value):
        return [0.0, 1.0]
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


def write_text(
    parts: List[str],
    x: float,
    y: float,
    text: object,
    size: int = 12,
    fill: str = BLACK,
    anchor: str = "start",
    weight: str = "400",
) -> None:
    parts.append(
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" '
        f'font-family="Arial, Helvetica, sans-serif" fill="{fill}" '
        f'text-anchor="{anchor}" font-weight="{weight}">{esc(text)}</text>'
    )


def line_path(xs: Sequence[float], ys: Sequence[float], x_to_px, y_to_px) -> str:
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


def draw_comparison_line_plot(
    out_path: Path,
    title: str,
    ylabel: str,
    runs: Sequence[Run],
    metric: str,
    note: str,
    y_limits: Optional[Tuple[float, float]] = None,
    width: int = 1040,
    height: int = 640,
) -> None:
    margin_left = 82
    margin_right = 38
    margin_top = 88
    margin_bottom = 116
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    all_epochs = finite([x for run in runs for x in column(run.metrics, "epoch")])
    all_values = finite([v for run in runs for v in column(run.metrics, metric)])
    if not all_epochs or not all_values:
        raise ValueError(f"No finite values for metric {metric}")

    x_min, x_max = min(all_epochs), max(all_epochs)
    if x_min == x_max:
        x_min -= 1
        x_max += 1

    if y_limits:
        y_min, y_max = y_limits
    else:
        y_min, y_max = min(all_values), max(all_values)
        pad = (y_max - y_min) * 0.08 or 0.1
        y_min -= pad
        y_max += pad
        if min(all_values) >= 0 and y_min > 0:
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
    write_text(parts, margin_left, 60, note, size=13, fill=GRAY)
    parts.append(
        f'<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" fill="{PANEL}" stroke="#cbd5e1"/>'
    )

    for tick in y_ticks:
        y = y_to_px(tick)
        parts.append(f'<line x1="{margin_left}" x2="{margin_left + plot_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="{GRID}"/>')
        write_text(parts, margin_left - 10, y + 4, fmt(tick), size=11, fill=GRAY, anchor="end")

    for tick in x_ticks:
        if tick < x_min - 1e-9 or tick > x_max + 1e-9:
            continue
        x = x_to_px(tick)
        parts.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{margin_top}" y2="{margin_top + plot_h}" stroke="#eef2f7"/>')
        write_text(parts, x, margin_top + plot_h + 24, fmt(tick, digits=0), size=11, fill=GRAY, anchor="middle")

    for run in runs:
        xs = column(run.metrics, "epoch")
        ys = column(run.metrics, metric)
        path = line_path(xs, ys, x_to_px, y_to_px)
        parts.append(
            f'<path d="{path}" fill="none" stroke="{run.color}" stroke-width="2.8" stroke-linejoin="round"/>'
        )
        best_y = value_at_epoch(run.metrics, metric, run.best_epoch)
        if math.isfinite(best_y):
            parts.append(
                f'<circle cx="{x_to_px(run.best_epoch):.2f}" cy="{y_to_px(best_y):.2f}" r="5" '
                f'fill="#ffffff" stroke="{run.color}" stroke-width="2.5"/>'
            )

    write_text(parts, margin_left + plot_w / 2, height - 34, "epoch", size=13, anchor="middle")
    write_text(parts, 24, margin_top + plot_h / 2, ylabel, size=13, anchor="middle")
    parts[-1] = parts[-1].replace(
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}"',
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}" transform="rotate(-90 24 {margin_top + plot_h / 2:.2f})"',
    )

    legend_y = height - 78
    legend_x = margin_left
    for i, run in enumerate(runs):
        row = i // 3
        col = i % 3
        x = legend_x + col * 300
        y = legend_y + row * 24
        parts.append(f'<line x1="{x}" x2="{x + 28}" y1="{y}" y2="{y}" stroke="{run.color}" stroke-width="3"/>')
        parts.append(f'<circle cx="{x + 14}" cy="{y}" r="4.5" fill="#ffffff" stroke="{run.color}" stroke-width="2"/>')
        write_text(parts, x + 36, y + 4, f"{run.label} (best {run.best_epoch})", size=12)

    parts.append("</svg>")
    out_path.write_text("\n".join(parts))


def draw_test_bar_plot(
    out_path: Path,
    title: str,
    runs: Sequence[Run],
    metrics: Sequence[Tuple[str, str]],
    ylabel: str,
    note: str,
    width: int = 1100,
    height: int = 660,
) -> None:
    margin_left = 86
    margin_right = 36
    margin_top = 88
    margin_bottom = 144
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    values = [
        run.test_metrics.get(metric, math.nan)
        for run in runs
        for metric, _ in metrics
    ]
    all_values = finite(values)
    if not all_values:
        raise ValueError("No finite test metric values for bar plot")

    y_min = 0.0
    y_max = max(all_values) * 1.16 or 1.0
    if max(all_values) <= 1.0:
        y_max = 1.0
    y_ticks = nice_ticks(y_min, y_max)
    y_min, y_max = min(y_ticks), max(y_ticks)

    def y_to_px(y: float) -> float:
        return margin_top + (y_max - y) / (y_max - y_min) * plot_h

    n_groups = len(metrics)
    group_w = plot_w / max(1, n_groups)
    bar_gap = 4
    bar_w = min(40, (group_w * 0.72 - bar_gap * (len(runs) - 1)) / max(1, len(runs)))

    parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    write_text(parts, margin_left, 34, title, size=24, weight="700")
    write_text(parts, margin_left, 60, note, size=13, fill=GRAY)
    parts.append(
        f'<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" fill="{PANEL}" stroke="#cbd5e1"/>'
    )

    for tick in y_ticks:
        y = y_to_px(tick)
        parts.append(f'<line x1="{margin_left}" x2="{margin_left + plot_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="{GRID}"/>')
        write_text(parts, margin_left - 10, y + 4, fmt(tick), size=11, fill=GRAY, anchor="end")

    for group_i, (metric, label) in enumerate(metrics):
        center = margin_left + group_w * (group_i + 0.5)
        total_bar_w = len(runs) * bar_w + (len(runs) - 1) * bar_gap
        start = center - total_bar_w / 2
        for run_i, run in enumerate(runs):
            value = run.test_metrics.get(metric, math.nan)
            if not math.isfinite(value):
                continue
            x = start + run_i * (bar_w + bar_gap)
            y = y_to_px(value)
            h = margin_top + plot_h - y
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{h:.2f}" '
                f'fill="{run.color}" rx="2"/>'
            )
            write_text(parts, x + bar_w / 2, y - 7, fmt(value), size=10, anchor="middle")

        label_lines = label.split("\\n")
        for j, line in enumerate(label_lines):
            write_text(parts, center, margin_top + plot_h + 26 + j * 14, line, size=11, fill=GRAY, anchor="middle")

    write_text(parts, 24, margin_top + plot_h / 2, ylabel, size=13, anchor="middle")
    parts[-1] = parts[-1].replace(
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}"',
        f'x="24.00" y="{margin_top + plot_h / 2:.2f}" transform="rotate(-90 24 {margin_top + plot_h / 2:.2f})"',
    )

    legend_y = height - 56
    for i, run in enumerate(runs):
        row = i // 3
        col = i % 3
        x = margin_left + col * 310
        y = legend_y + row * 24
        parts.append(f'<rect x="{x}" y="{y - 12}" width="14" height="14" fill="{run.color}"/>')
        write_text(parts, x + 22, y, run.label, size=12)

    parts.append("</svg>")
    out_path.write_text("\n".join(parts))


def load_runs(run_args: Sequence[str]) -> List[Run]:
    if not run_args:
        raise ValueError("Provide at least one --run LABEL=RUN_DIR.")
    runs: List[Run] = []
    for i, raw in enumerate(run_args):
        label, path = parse_run_arg(raw)
        metrics_path = path / "metrics.csv"
        if not metrics_path.exists():
            raise FileNotFoundError(f"Missing metrics file: {metrics_path}")
        metrics = read_metrics(metrics_path)
        runs.append(
            Run(
                label=label,
                path=path,
                color=COLORS[i % len(COLORS)],
                metrics=metrics,
                test_metrics=read_json(path / "test_metrics.json"),
                config=read_json(path / "run_config.json"),
                best_epoch=best_epoch(metrics),
            )
        )
    return runs


def write_summary_csv(out_path: Path, runs: Sequence[Run]) -> None:
    fields = [
        "label",
        "run_dir",
        "dataset_dir",
        "edge_scope",
        "model",
        "kingdom",
        "best_epoch",
        "best_val_loss",
        "best_val_presence_auroc",
        "best_val_presence_average_precision",
        "best_val_log_rmse_present",
        "test_loss",
        "test_presence_auroc",
        "test_presence_average_precision",
        "test_log_rmse_present",
        "test_log_mae_present",
        "test_presence_rate",
        "test_presence_predicted_rate",
        "resident_log_rmse_present",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            args = run.config.get("args", {})
            features = run.config.get("feature_manifest", {})
            writer.writerow(
                {
                    "label": run.label,
                    "run_dir": str(run.path),
                    "dataset_dir": args.get("dataset_dir", ""),
                    "edge_scope": features.get("edge_scope", ""),
                    "model": args.get("model", ""),
                    "kingdom": args.get("kingdom", ""),
                    "best_epoch": run.best_epoch,
                    "best_val_loss": value_at_epoch(run.metrics, "val_loss", run.best_epoch),
                    "best_val_presence_auroc": value_at_epoch(run.metrics, "val_presence_auroc", run.best_epoch),
                    "best_val_presence_average_precision": value_at_epoch(
                        run.metrics, "val_presence_average_precision", run.best_epoch
                    ),
                    "best_val_log_rmse_present": value_at_epoch(run.metrics, "val_log_rmse_present", run.best_epoch),
                    "test_loss": run.test_metrics.get("loss", math.nan),
                    "test_presence_auroc": run.test_metrics.get("presence_auroc", math.nan),
                    "test_presence_average_precision": run.test_metrics.get("presence_average_precision", math.nan),
                    "test_log_rmse_present": run.test_metrics.get("log_rmse_present", math.nan),
                    "test_log_mae_present": run.test_metrics.get("log_mae_present", math.nan),
                    "test_presence_rate": run.test_metrics.get("presence_rate", math.nan),
                    "test_presence_predicted_rate": run.test_metrics.get("presence_predicted_rate", math.nan),
                    "resident_log_rmse_present": run.test_metrics.get("resident_log_rmse_present", math.nan),
                }
            )


def write_readme(out_path: Path, title: str, runs: Sequence[Run], plot_files: Sequence[Path]) -> None:
    lines = [
        f"# {title}",
        "",
        "These plots compare training curves and held-out test metrics across GNN run directories.",
        "",
        "## Runs",
        "",
    ]
    for run in runs:
        args = run.config.get("args", {})
        features = run.config.get("feature_manifest", {})
        lines.extend(
            [
                f"- `{run.label}`",
                f"  - run directory: `{run.path}`",
                f"  - model: `{args.get('model', 'unknown')}`",
                f"  - dataset: `{args.get('dataset_dir', 'unknown')}`",
                f"  - edge scope: `{features.get('edge_scope', 'unknown')}`",
                f"  - best epoch: `{run.best_epoch}`",
                f"  - test AP: `{fmt(run.test_metrics.get('presence_average_precision'))}`",
                f"  - test present-taxon RMSE: `{fmt(run.test_metrics.get('log_rmse_present'))}`",
            ]
        )

    lines.extend(
        [
            "",
            "## Plot Files",
            "",
        ]
    )
    for plot_file in plot_files:
        lines.append(f"- `{plot_file.name}`")
    lines.extend(
        [
            "",
            "## Reading The Comparison",
            "",
            "- Similar SAGE performance with resident edges and no edges means donor/resident node features are carrying most of the signal.",
            "- If GINE with resident edges improves over SAGE/no-edge runs, edge weights and signs may be adding useful information.",
            "- Average precision is more informative than AUROC for sparse final-presence targets.",
            "- Predicted presence rate should be compared against true presence rate before using a fixed 0.5 threshold biologically.",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    runs = load_runs(args.run)

    plots: List[Path] = []
    plot_specs = [
        (
            "01_validation_loss.svg",
            "Validation Loss",
            "val_loss",
            "loss",
            None,
            "Lower is better; circles mark each run's best validation-loss epoch.",
        ),
        (
            "02_validation_presence_average_precision.svg",
            "Validation Final-Presence Average Precision",
            "val_presence_average_precision",
            "average precision",
            (0, 1),
            "Higher is better; this is the key sparse-target ranking metric.",
        ),
        (
            "03_validation_presence_auroc.svg",
            "Validation Final-Presence AUROC",
            "val_presence_auroc",
            "AUROC",
            (0, 1),
            "Higher is better, but AUROC can look optimistic when positives are rare.",
        ),
        (
            "04_validation_present_taxa_rmse.svg",
            "Validation Log-Abundance RMSE On Present Taxa",
            "val_log_rmse_present",
            "log RMSE",
            None,
            "Lower is better for taxa that are actually present in final communities.",
        ),
        (
            "05_validation_predicted_presence_rate.svg",
            "Validation Predicted Presence Rate At Threshold 0.5",
            "val_presence_predicted_rate",
            "fraction predicted present",
            None,
            "Compare with the true validation presence rate; high ranking can still need threshold calibration.",
        ),
        (
            "06_training_loss.svg",
            "Training Loss",
            "train_loss",
            "loss",
            None,
            "Useful for seeing whether models are still optimizing or overfitting.",
        ),
    ]

    for filename, title, metric, ylabel, limits, note in plot_specs:
        path = out_dir / filename
        draw_comparison_line_plot(
            path,
            f"{args.title}: {title}",
            ylabel,
            runs,
            metric,
            note,
            y_limits=limits,
        )
        plots.append(path)

    test_presence_plot = out_dir / "07_test_presence_metrics.svg"
    draw_test_bar_plot(
        test_presence_plot,
        f"{args.title}: Held-Out Presence Metrics",
        runs,
        [
            ("presence_auroc", "AUROC"),
            ("presence_average_precision", "avg\\nprecision"),
            ("presence_rate", "true\\npresence rate"),
            ("presence_predicted_rate", "predicted\\npresence rate"),
        ],
        "score or fraction",
        "Average precision should be interpreted relative to the true presence rate.",
    )
    plots.append(test_presence_plot)

    test_abundance_plot = out_dir / "08_test_abundance_metrics.svg"
    draw_test_bar_plot(
        test_abundance_plot,
        f"{args.title}: Held-Out Abundance Error",
        runs,
        [
            ("log_rmse_present", "GNN RMSE\\npresent"),
            ("log_mae_present", "GNN MAE\\npresent"),
            ("resident_log_rmse_present", "resident RMSE\\npresent"),
            ("resident_log_mae_present", "resident MAE\\npresent"),
        ],
        "log-scale error",
        "Lower is better; resident columns are the final = resident baseline.",
    )
    plots.append(test_abundance_plot)

    write_summary_csv(out_dir / "comparison_summary.csv", runs)
    write_readme(out_dir / "README.md", args.title, runs, plots)

    print(f"Wrote {len(plots)} comparison plots to {out_dir}")
    print(f"Summary CSV: {out_dir / 'comparison_summary.csv'}")
    print(f"README: {out_dir / 'README.md'}")


if __name__ == "__main__":
    main()
