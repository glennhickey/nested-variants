#!/usr/bin/env python3
"""Compose multiple PNG panels into a labeled grid summary figure.

Usage:
    python3 scripts/compose-summary.py \
        --output OUT.png \
        --title "Figure Title" \
        --cols 2 \
        --panels "Segment Lengths:path/a.png" "Density:path/b.png" ...

Missing panel files are silently skipped (grid contracts).
"""

import argparse
import os
import string
import sys

from PIL import Image, ImageDraw, ImageFont

# Layout constants (pixels at 300 dpi)
GUTTER = 40
TITLE_HEIGHT = 80
LABEL_PAD = 10
CANVAS_WIDTH = 4800  # 16 inches at 300 dpi


def load_font(size):
    """Try system fonts, fall back to Pillow default."""
    for path in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
    ]:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def main():
    parser = argparse.ArgumentParser(description="Compose PNG panels into grid")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--title", default="Summary", help="Figure title")
    parser.add_argument("--cols", type=int, default=2, help="Grid columns")
    parser.add_argument("--panels", nargs="+", required=True,
                        help="Panel specs as 'Label:filepath'")
    args = parser.parse_args()

    # Parse and filter panels to those that exist
    panels = []
    for spec in args.panels:
        if ":" not in spec:
            continue
        label, path = spec.split(":", 1)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            panels.append((label.strip(), path))
        else:
            print(f"  Skipping missing: {path}", file=sys.stderr)

    if not panels:
        print("No panels found — creating placeholder.", file=sys.stderr)
        img = Image.new("RGB", (800, 200), "white")
        draw = ImageDraw.Draw(img)
        draw.text((20, 80), "No panels available", fill="grey")
        img.save(args.output)
        return

    # Load panel images
    images = []
    for label, path in panels:
        images.append((label, Image.open(path).convert("RGB")))

    # Compute cell width
    cols = min(args.cols, len(images))
    cell_w = (CANVAS_WIDTH - (cols + 1) * GUTTER) // cols

    # Scale panels and compute row heights
    rows = []
    row = []
    for i, (label, img) in enumerate(images):
        scale = cell_w / img.width
        new_h = int(img.height * scale)
        scaled = img.resize((cell_w, new_h), Image.LANCZOS)
        row.append((label, scaled))
        if len(row) == cols:
            rows.append(row)
            row = []
    if row:
        rows.append(row)

    # Compute canvas height
    row_heights = [max(img.height for _, img in r) for r in rows]
    canvas_h = (TITLE_HEIGHT + GUTTER
                + sum(row_heights)
                + (len(rows) - 1) * GUTTER
                + GUTTER)

    # Create canvas
    canvas = Image.new("RGB", (CANVAS_WIDTH, canvas_h), "white")
    draw = ImageDraw.Draw(canvas)

    # Draw title
    title_font = load_font(48)
    bbox = draw.textbbox((0, 0), args.title, font=title_font)
    tw = bbox[2] - bbox[0]
    draw.text(((CANVAS_WIDTH - tw) // 2, GUTTER // 2), args.title,
              fill="black", font=title_font)

    # Draw panels
    label_font = load_font(36)
    caption_font = load_font(28)
    letter_idx = 0
    y = TITLE_HEIGHT + GUTTER
    for ri, row in enumerate(rows):
        x = GUTTER
        for label, img in row:
            # Paste panel
            canvas.paste(img, (x, y))

            # Panel letter label (A, B, C, ...)
            letter = string.ascii_uppercase[letter_idx]
            letter_idx += 1
            draw.text((x + LABEL_PAD, y + LABEL_PAD),
                      letter, fill="black", font=label_font)

            # Caption below letter
            draw.text((x + LABEL_PAD + 40, y + LABEL_PAD + 4),
                      label, fill="grey", font=caption_font)

            x += cell_w + GUTTER
        y += row_heights[ri] + GUTTER

    # Close source images
    for _, img in images:
        img.close()

    canvas.save(args.output, dpi=(300, 300))
    print(f"Saved: {args.output} ({CANVAS_WIDTH}x{canvas_h} px, {len(panels)} panels)")


if __name__ == "__main__":
    main()
