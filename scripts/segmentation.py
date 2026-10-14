#!/usr/bin/env python3
"""
segmentation.py - Per-channel cell segmentation using Mesmer (DeepCell)

For each targeted channel:
  - Segment cells independently using Mesmer
  - Save a label mask: [filename]_[channel]_mask.tif

Output table (segmentation_results.csv):
  - One row per detected cell per channel
  - Columns: region, cell_index, area, x, y, source_channel
    - region:         filename without extension
    - cell_index:     cell label ID within that channel's mask (1-based)
    - area:           number of pixels in the cell mask
    - x, y:           centroid coordinates (pixels, floating point)
    - source_channel: which channel this cell was segmented from

Usage:
    python segmentation.py --input_dir TMA1 --channels DRAQ5,CD45
    python segmentation.py --input_dir TMA1           # all non-blank channels
    python segmentation.py --input_dir TMA1 --channels DRAQ5 --output_dir out/
"""

import argparse
import logging
import pathlib
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-channel cell segmentation using Mesmer (DeepCell)"
    )
    parser.add_argument(
        "--input_dir", type=pathlib.Path, required=True,
        help="Directory containing .tif files and channelnames.txt"
    )
    parser.add_argument(
        "--channels", type=str, default=None,
        help="Comma-separated channel names (e.g. 'DRAQ5,CD45'). "
             "Default: all non-Blank, non-empty channels."
    )
    parser.add_argument(
        "--output_dir", type=pathlib.Path, default=None,
        help="Directory for mask TIFs and results CSV. Defaults to --input_dir."
    )
    parser.add_argument(
        "--image_mpp", type=float, default=0.5,
        help="Microns per pixel of the input images (used by Mesmer). Default: 0.5"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(output_dir: pathlib.Path) -> logging.Logger:
    logger = logging.getLogger("segmentation")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    log_path = output_dir / "segmentation.log"
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# ---------------------------------------------------------------------------
# Channel map
# ---------------------------------------------------------------------------

def read_channel_map(input_dir: pathlib.Path) -> Dict[str, List[int]]:
    """
    Parse channelnames.txt into {channel_name: [0-based_page_indices]}.

    Handles two formats:
      - Plain:          one name per line   ("DRAQ5")
      - Tab-separated:  "index\\tname"       ("72\\tDRAQ5")
    """
    channel_file = input_dir / "channelnames.txt"
    if not channel_file.exists():
        raise FileNotFoundError(f"channelnames.txt not found in {input_dir}")

    name_to_indices: Dict[str, List[int]] = defaultdict(list)
    with open(channel_file, "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            line = line.rstrip("\n\r")
            name = line.split("\t")[-1].strip() if "\t" in line else line.strip()
            name_to_indices[name].append(idx)

    return dict(name_to_indices)


def resolve_channels(
    channel_map: Dict[str, List[int]],
    requested: Optional[str],
    logger: logging.Logger,
) -> List[Tuple[str, int]]:
    """
    Returns ordered list of (channel_name, 0-based_page_index) to process.
    Duplicates: first occurrence used with a warning.
    requested=None: all channels except empty strings and 'Blank'.
    """
    EXCLUDE = {"", "Blank"}

    def first_idx(name: str) -> int:
        indices = channel_map[name]
        if len(indices) > 1:
            logger.warning(
                f"Channel '{name}' appears {len(indices)} times "
                f"(pages {indices}); using first occurrence (page {indices[0]})."
            )
        return indices[0]

    if requested is None:
        seen = set()
        result = []
        for name in channel_map:
            if name in EXCLUDE or name in seen:
                continue
            seen.add(name)
            result.append((name, first_idx(name)))
        result.sort(key=lambda x: x[1])
        return result

    result = []
    lower_map = {k.lower(): k for k in channel_map}
    for raw in (c.strip() for c in requested.split(",")):
        if not raw:
            continue
        if raw in channel_map:
            result.append((raw, first_idx(raw)))
        elif raw.lower() in lower_map:
            canonical = lower_map[raw.lower()]
            logger.warning(f"Channel '{raw}' matched '{canonical}' case-insensitively.")
            result.append((canonical, first_idx(canonical)))
        else:
            available = sorted(k for k in channel_map if k not in EXCLUDE)
            raise ValueError(
                f"Channel '{raw}' not found in channelnames.txt. "
                f"Available: {available}"
            )
    return result


# ---------------------------------------------------------------------------
# Image I/O
# ---------------------------------------------------------------------------

def read_channel_image(tif_path: pathlib.Path, page_index: int) -> np.ndarray:
    """Read one page from a multi-page TIFF. Returns 2D uint16 array."""
    with tifffile.TiffFile(tif_path) as tf:
        n_pages = len(tf.pages)
        if page_index >= n_pages:
            raise IndexError(
                f"Page {page_index} out of range for {tif_path.name} ({n_pages} pages)."
            )
        img = tf.pages[page_index].asarray()
    if img.ndim != 2:
        raise ValueError(
            f"Expected 2D from page {page_index}, got shape {img.shape}."
        )
    return img


def save_mask(mask: np.ndarray, output_path: pathlib.Path, logger: logging.Logger) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(output_path), mask.astype(np.uint16), photometric="minisblack")
    logger.info(f"  Saved mask: {output_path.name} ({int(mask.max())} cells)")


def save_overlay(
    img: np.ndarray,
    mask: np.ndarray,
    output_path: pathlib.Path,
    logger: logging.Logger,
) -> None:
    """
    Save an RGB TIFF with cell boundary lines drawn in red over the grayscale channel image.
    The channel image is normalized to uint8 for display.
    """
    from skimage.segmentation import find_boundaries

    gray = normalize_to_uint8(img)
    rgb = np.stack([gray, gray, gray], axis=-1)  # (H, W, 3) uint8

    borders = find_boundaries(mask, mode="inner")  # True where cell edge pixels are
    rgb[borders, 0] = 255  # red channel max
    rgb[borders, 1] = 0
    rgb[borders, 2] = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(output_path), rgb, photometric="rgb")
    logger.info(f"  Saved overlay: {output_path.name}")


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

def normalize_to_uint8(img: np.ndarray) -> np.ndarray:
    """
    Normalize a 16-bit image to uint8 [0, 255] for Mesmer input.
    Uses percentile clipping to handle outlier bright pixels.
    """
    p1 = np.percentile(img, 1)
    p99 = np.percentile(img, 99)
    if p99 <= p1:
        return np.zeros(img.shape, dtype=np.uint8)
    clipped = np.clip(img, p1, p99)
    scaled = (clipped - p1) / (p99 - p1) * 255.0
    return scaled.astype(np.uint8)


# ---------------------------------------------------------------------------
# Mesmer segmentation
# ---------------------------------------------------------------------------

def run_mesmer(
    img: np.ndarray,
    app,
    image_mpp: float,
    logger: logging.Logger,
) -> np.ndarray:
    """
    Run Mesmer on a single 2D channel image.
    Mesmer expects input shape (batch, H, W, 2): [nuclear_channel, membrane_channel].
    We pass the target channel as nuclear and a zero array as membrane.

    Returns a 2D uint16 label mask.
    """
    h, w = img.shape
    nuclear = normalize_to_uint8(img).astype(np.float32)
    membrane = np.zeros((h, w), dtype=np.float32)

    # Mesmer input: (1, H, W, 2)
    inp = np.stack([nuclear, membrane], axis=-1)[np.newaxis]

    result = app.predict(inp, image_mpp=image_mpp)
    # result shape: (1, H, W, 2) — channel 0 = whole-cell, channel 1 = nuclear
    # We use whole-cell mask (channel 0)
    mask = result[0, :, :, 0].astype(np.uint16)
    return mask


# ---------------------------------------------------------------------------
# Cell measurements
# ---------------------------------------------------------------------------

def measure_cells(mask: np.ndarray, stem: str, ch_name: str) -> List[Dict]:
    """
    Extract area and centroid for every cell in a label mask.
    Returns a list of dicts with keys: region, cell_index, area, x, y, source_channel.
    """
    rows = []
    for prop in regionprops(mask):
        cy, cx = prop.centroid  # regionprops returns (row, col) = (y, x)
        rows.append({
            "region": stem,
            "cell_index": int(prop.label),
            "area": int(prop.area),
            "x": round(cx, 2),
            "y": round(cy, 2),
            "source_channel": ch_name,
        })
    return rows


# ---------------------------------------------------------------------------
# Per-TIF processing
# ---------------------------------------------------------------------------

def process_tif(
    tif_path: pathlib.Path,
    channels_to_process: List[Tuple[str, int]],
    output_dir: pathlib.Path,
    image_mpp: float,
    app,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Segment all requested channels for one TIF. Returns per-cell DataFrame.

    DataFrame columns: region, cell_index, area, x, y, source_channel
    """
    logger.info(f"\n=== {tif_path.name} ===")
    stem = tif_path.stem
    rows = []

    for ch_name, page_idx in tqdm(channels_to_process, desc=stem, unit="ch", leave=False):
        mask_path = output_dir / f"{stem}_{ch_name}_mask.tif"
        overlay_path = output_dir / f"{stem}_{ch_name}_overlay.tif"

        # Always read the channel image — needed for the overlay regardless of cache
        try:
            img = read_channel_image(tif_path, page_idx)
            logger.info(
                f"  [{ch_name}] page={page_idx}, shape={img.shape}, "
                f"range=[{img.min()}, {img.max()}]"
            )
        except Exception as exc:
            logger.error(f"  [{ch_name}] image read FAILED: {exc}")
            continue

        # Load existing mask or run segmentation
        if mask_path.exists():
            logger.info(f"  [{ch_name}] mask exists, loading.")
            mask = tifffile.imread(str(mask_path)).astype(np.uint16)
        else:
            try:
                mask = run_mesmer(img, app, image_mpp, logger)
                save_mask(mask, mask_path, logger)
            except Exception as exc:
                logger.error(f"  [{ch_name}] segmentation FAILED: {exc}")
                continue

        # Save overlay (always regenerate so it reflects the current mask)
        save_overlay(img, mask, overlay_path, logger)

        n_cells = int(mask.max())
        logger.info(f"  [{ch_name}] {n_cells} cells — measuring properties.")
        rows.extend(measure_cells(mask, stem, ch_name))

    if not rows:
        logger.warning(f"  No cells detected in any channel for {tif_path.name}.")
        return pd.DataFrame(columns=["region", "cell_index", "area", "x", "y", "source_channel"])

    return pd.DataFrame(rows)[["region", "cell_index", "area", "x", "y", "source_channel"]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if args.output_dir:
        if args.output_dir.exists():
            print(
                f"ERROR: Output directory already exists: {args.output_dir}\n"
                "Remove or rename it, or choose a different --output_dir.",
                file=sys.stderr,
            )
            sys.exit(1)
        output_dir = args.output_dir
        output_dir.mkdir(parents=True)
    else:
        output_dir = args.input_dir

    output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(output_dir)
    logger.info(f"Input dir : {args.input_dir}")
    logger.info(f"Output dir: {output_dir}")
    logger.info(f"Image MPP : {args.image_mpp}")

    # Load Mesmer (import here so errors surface clearly if deepcell not installed)
    logger.info("Loading Mesmer model (first run downloads weights ~500 MB)...")
    try:
        from deepcell.applications import Mesmer
        app = Mesmer()
        logger.info("Mesmer loaded.")
    except ImportError:
        logger.error(
            "deepcell is not installed. "
            "Run: conda env create -f environment.yml && conda activate segmentation"
        )
        sys.exit(1)

    # Channel map
    channel_map = read_channel_map(args.input_dir)
    total_entries = sum(len(v) for v in channel_map.values())
    logger.info(
        f"channelnames.txt: {total_entries} entries, {len(channel_map)} unique names"
    )

    # Resolve targeted channels
    channels_to_process = resolve_channels(channel_map, args.channels, logger)
    channel_names = [ch for ch, _ in channels_to_process]
    logger.info(f"Channels to segment ({len(channels_to_process)}): {channel_names}")

    # Discover source TIF files (exclude already-created masks)
    tif_files = sorted(
        p for p in args.input_dir.glob("*.tif")
        if not p.name.endswith("_mask.tif") and not p.name.endswith("_overlay.tif")
    )
    if not tif_files:
        logger.error(f"No TIF files found in {args.input_dir}")
        sys.exit(1)
    logger.info(f"TIF files ({len(tif_files)}): {[p.name for p in tif_files]}")

    # Process each TIF and collect results
    all_dfs = []
    for tif_path in tqdm(tif_files, desc="TIF files", unit="file"):
        df = process_tif(
            tif_path=tif_path,
            channels_to_process=channels_to_process,
            output_dir=output_dir,
            image_mpp=args.image_mpp,
            app=app,
            logger=logger,
        )
        all_dfs.append(df)

    # Save combined results CSV
    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        csv_path = output_dir / "segmentation_results.csv"
        combined.to_csv(csv_path, index=False)
        logger.info(f"\nResults CSV saved: {csv_path}")
        logger.info(f"Total rows: {len(combined)}")
        logger.info(f"\n{combined.head(10).to_string()}")
    else:
        logger.warning("No results to save.")

    logger.info("Done.")


if __name__ == "__main__":
    main()
