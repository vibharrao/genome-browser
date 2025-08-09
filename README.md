# Gene Locus Visualization

This project is a Python-based plotting script that generates a multi-panel genomic visualization of a specified gene locus from annotation (GTF) and alignment (PSL) files.  
It replicates a provided template figure at 2400 DPI using only matplotlib's rectangle primitives.

---

## Overview

Given a chromosome region (via the `-c` flag), the script:

1. **Parses the GTF annotation file** (`gencode.vM12.annotation.gtf`) to identify and plot gene models for that locus, including:
   - **Exons** (short rectangles)
   - **Coding Sequences (CDS)** (taller rectangles)
   - **Introns** (thin connecting bars between exons/CDS)

2. **Parses PSL alignment files**:
   - `BME163_Input_Data_5.psl` → Plotted in **orange** (second panel from top), sorted by alignment **end position**.
   - `BME163_Input_Data_6.psl` → Plotted twice in **blue** (third and fourth panels), sorted by alignment **start position**.

3. **Stacks reads efficiently** within each panel to minimize vertical space usage.

4. **Generates a read coverage histogram** (bottom panel) from the alignments in the panel above, using 1 bp bins.

---

## Figure Layout

From top to bottom:

1. **Panel 3.9–5.4 inches**: Gene model track from GTF.
2. **Panel 2.2–3.7 inches**: PSL reads (orange), sorted by end coordinate.
3. **Panel 0.5–2.0 inches**: PSL reads (blue), sorted by start coordinate.
4. **Panel 0.1–0.4 inches**: Read coverage histogram (blue bars, no outline).

---

## Features

- Accepts any locus in `chr:start-end` format with `-c` flag.
- Maintains correct proportions for exon/CDS/intron blocks.
- Colors and line widths match the assignment specification:
  - Gene model: grey outlines, linewidth 0.25
  - Orange reads: (230/255, 87/255, 43/255), linewidth default
  - Blue reads: (88/255, 85/255, 120/255), linewidth 0.05
- Automatic stacking to avoid overlaps in each track.

---

## Usage

```bash
python Rao_Vibha_BME163_Assignment_Final.py \
    -p5 /path/to/BME163_Input_Data_5.psl \
    -p6 /path/to/BME163_Input_Data_6.psl \
    -g /path/to/gencode.vM12.annotation.gtf \
    -c chr7:45232000-45241000 \
    -o output.png
