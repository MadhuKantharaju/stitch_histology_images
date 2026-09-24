# Stitch Histology Images

Utilities for diagnosing and correcting misaligned histology-image tile registrations when Fiji Grid/Collection Stitching produces inaccurate results.

The repository contains two complementary tools:

- `fix_sequential_jumps.py` detects unusually large positional jumps between sequentially acquired tiles and realigns misplaced tiles or contiguous blocks by translation. Ambiguous cases can be reviewed interactively and adjusted with translation and rotation controls.
- `topology_graph.py` analyzes registered tile positions, checks the topology of the scan path, and creates an interactive Plotly HTML graph showing trusted, rejected, and suspicious links.

## Features

- Detects unexpected jumps between sequential tiles using the numeric suffix in each filename.
- Corrects misplaced tiles or blocks in a Fiji `TileConfiguration.registered.txt` file.
- Supports whole-block translation, rotation, and optional per-tile nudging during interactive review.
- Checks for disconnected components, unexpected spatial neighbors, acquisition-order violations, and crossing edges.
- Produces interactive topology visualizations and CSV reports of detected issues.
- Supports open paths and closed-loop acquisitions.

## Requirements

Python 3.8 or newer is recommended.

For topology visualization:

```bash
pip install plotly
```

For the interactive correction interface:

```bash
pip install matplotlib
```

`tkinter` is also required for the graphical review windows. It is included with many Python installations, but on some Linux distributions it must be installed separately, for example:

```bash
sudo apt-get install python3-tk
```

The scripts can still run in a text-based mode when the graphical dependencies are unavailable.

## Correct misplaced tiles or blocks

### Using a tile path

If the tile is inside a sequence directory containing `TileConfiguration.registered.txt`, the script can derive the input file automatically:

```bash
python fix_sequential_jumps.py \
  --tile-path "/path/to/project/slices/MySequence/Snap-1.czi"
```

The corrected configuration is written as `TileConfiguration.jumpfixed.txt` in the sequence directory.

### Using a sequence directory

```bash
python fix_sequential_jumps.py \
  --sequence-dir "/path/to/project/slices/MySequence" \
  --out "TileConfiguration.jumpfixed.txt"
```

### Using an explicit registered configuration

```bash
python fix_sequential_jumps.py \
  --registered "/path/to/TileConfiguration.registered.txt" \
  --out "/path/to/TileConfiguration.jumpfixed.txt"
```

### Interactive review

By default, blocks with insufficient or conflicting evidence are reported but are not automatically changed. Use `--interactive` to review those blocks in a graphical interface:

```bash
python fix_sequential_jumps.py \
  --registered "/path/to/TileConfiguration.registered.txt" \
  --interactive
```

Useful options:

- `--jump-threshold-multiplier`: Flag a sequential step when its distance exceeds this multiple of the typical step distance. Default: `3.0`.
- `--interactive`: Review ambiguous blocks manually.
- `--loop`: Mark the acquisition as a closed loop. This option is reserved for circular handling in the correction script.

## Analyze a registered topology

### Analyze one sequence from a tile path

```bash
python topology_graph.py \
  --tile-path "/path/to/project/slices/MySequence/Snap-1.czi" \
  --out "topology.html"
```

The script derives the project, sequence, log, and registered-configuration paths from the tile path. It writes an interactive HTML graph and a companion issues CSV file.

### Analyze using explicit paths

```bash
python topology_graph.py \
  --log-dir "/path/to/logs" \
  --registered "/path/to/TileConfiguration.registered.txt" \
  --out "topology.html"
```

### Analyze all sequences in a project

```bash
python topology_graph.py \
  --project-root "/path/to/project" \
  --out-dir "/path/to/output"
```

Batch mode expects the project to contain sequence folders and log files in the layout used by the script. Review the path definitions in `run_batch()` if your project uses different folder names.

Useful options:

- `--shift-discrepancy-threshold`: Pixel tolerance for comparing logged shifts with registered positions. Default: `50`.
- `--max-scan-jump`: Maximum expected numeric gap between linked filenames. Default: `3`.
- `--loop`: Treat the acquisition order as a closed loop.
- `--spatial-threshold-multiplier`: Multiplier used to detect unexpected spatially close tiles. Default: `1.5`.
- `--no-cluster-colors`: Disable cluster-based node coloring in the graph.

## Input format

The scripts read Fiji-style tile configuration files such as:

```text
# Define the number of dimensions we are working on
dim = 2

# Define the image coordinates
Snap-1.czi; ; (100.0000, 200.0000)
Snap-2.czi; ; (512.0000, 205.0000)
```

Tile filenames should include a trailing numeric acquisition or scan index, such as `Snap-17149.czi`. That number is used to determine the expected acquisition order.

Topology analysis also reads log lines containing tile pairs, shifts, and correlation values. Links are classified as:

- **Green solid**: trusted and geometrically consistent.
- **Orange dashed**: high-correlation but contradicted by the final registration.
- **Gray dotted**: rejected or below the correlation threshold.

## Outputs

`fix_sequential_jumps.py` produces:

- `TileConfiguration.jumpfixed.txt`: corrected tile positions.
- Console reports describing realigned and manually reviewed blocks.

`topology_graph.py` produces:

- An interactive Plotly HTML topology graph.
- An `_issues.csv` file containing detected structural and geometric problems.
- In batch mode, a `summary.csv` file containing one verdict per sequence.

## Workflow

1. Run Fiji stitching and generate the registered tile configuration.
2. Run `topology_graph.py` to inspect the registered geometry and identify suspicious jumps or crossings.
3. Run `fix_sequential_jumps.py` to correct confidently identified misplaced tiles or blocks.
4. Review ambiguous cases with `--interactive`.
5. Re-run the topology analysis to verify the corrected configuration.
6. Use the resulting `TileConfiguration.jumpfixed.txt` for downstream stitching or inspection.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
