"""
fix_sequential_jumps.py

Uses the same index-guided trace idea from topology_graph.py (tiles ordered
by their filename's scan number, not by blind nearest-distance) to find
tiles that are likely MISPLACED, then corrects just those tiles by
interpolating from their scan-order neighbors.

Key idea for telling WHICH tile is actually wrong:
A "jump" is flagged on the EDGE between two sequential tiles (e.g. tile 35
-> tile 36). If tile 36 is the one that's actually misplaced, it will show
up in TWO jumps: 35->36 AND 36->37. A tile that's only involved in one jump
is ambiguous (either it or its neighbor is the problem) and is reported for
manual review rather than auto-corrected.

For every tile flagged on both sides, its corrected position is the
midpoint of its two scan-order neighbors' actual positions (tile 35 and
tile 37, in the example above) -- i.e. "where it should sit given the
tiles that come immediately before and after it in the scan."

Usage (simplest):
    python fix_sequential_jumps.py --tile-path "/path/to/project/slices/MySeq/Snap-1.czi"

Or manually:
    python fix_sequential_jumps.py --registered "/path/to/TileConfiguration.registered.txt" --out "TileConfiguration.jumpfixed.txt"
"""

import argparse
import glob
import math
import os
import re
import sys

try:
    import tkinter as tk
    TKINTER_AVAILABLE = True
except ImportError:
    TKINTER_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def extract_scan_index(filename):
    """Pull the trailing acquisition number out of a filename like 'Snap-17149.czi' -> 17149."""
    m = re.search(r"(\d+)(?=\.[^.]+$)", filename)
    return int(m.group(1)) if m else None


def parse_tile_configuration(path):
    positions = {}
    if not path or not os.path.exists(path):
        return positions
    with open(path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith("dim"):
                continue
            parts = line.split(";")
            if len(parts) < 3:
                continue
            fname = parts[0].strip()
            coord_str = parts[2].strip().strip("()")
            try:
                coords = [float(c.strip()) for c in coord_str.split(",")]
            except ValueError:
                continue
            if len(coords) >= 2:
                positions[os.path.basename(fname)] = (coords[0], coords[1])
    return positions


def find_registered_file(seq_dir):
    matches = glob.glob(os.path.join(seq_dir, "**", "TileConfiguration.registered.txt"), recursive=True)
    return matches[0] if matches else None


def dist(a, b):
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


def build_ordered_tiles(registered, expect_loop=False):
    tiles_with_idx = [(t, extract_scan_index(t)) for t in registered]
    tiles_with_idx = [(t, i) for t, i in tiles_with_idx if i is not None]
    tiles_with_idx.sort(key=lambda x: x[1])
    ordered = [t for t, _ in tiles_with_idx]
    skipped = [t for t in registered if extract_scan_index(t) is None]
    if skipped:
        print(f"WARNING: {len(skipped)} tile(s) have no parseable scan-order number and were excluded from the trace: {skipped}")
    return ordered


def find_jump_edges(ordered, registered, jump_threshold_multiplier=3.0):
    """Return (jump_edge_indices, threshold, median_step). jump_edge_indices[k]=True means
    the step from ordered[k] to ordered[k+1] is a flagged jump."""
    n = len(ordered)
    step_dists = [dist(registered[ordered[k]], registered[ordered[k + 1]]) for k in range(n - 1)]
    if not step_dists:
        return [], None, None
    sorted_steps = sorted(step_dists)
    median_step = sorted_steps[len(sorted_steps) // 2]
    threshold = median_step * jump_threshold_multiplier
    flags = [d > threshold for d in step_dists]
    return flags, threshold, median_step


def split_into_segments(flags):
    """
    Split the tile sequence into segments (index ranges) at every flagged
    step. Each segment is a maximal run of tiles connected by NORMAL
    (non-flagged) steps -- this could be a single misplaced tile (segment
    length 1) or a whole block of tiles that's internally fine but shifted
    as a unit relative to its neighbors (segment length > 1). Both cases
    are handled the same way: realign the whole segment with one rigid
    translation.

    Returns a list of (start_idx, end_idx) inclusive index ranges into the
    ordered tile list.
    """
    n = len(flags) + 1  # number of tiles
    segments = []
    seg_start = 0
    for k, is_jump in enumerate(flags):
        if is_jump:
            segments.append((seg_start, k))
            seg_start = k + 1
    segments.append((seg_start, n - 1))
    return segments


def estimate_extrapolated_position(ordered, registered, from_idx, direction, n_steps_for_trend=3):
    """
    Look at the last few tiles ending at `from_idx` in the given direction
    (+1 = moving forward, -1 = moving backward) and extrapolate one more
    step in that same direction, using the AVERAGE step vector over the
    last n_steps_for_trend steps (falls back to fewer if not enough tiles).
    Returns the predicted (x, y) for the tile immediately beyond from_idx.
    """
    trend_points = []
    idx = from_idx
    for _ in range(n_steps_for_trend + 1):
        if idx < 0 or idx >= len(ordered):
            break
        trend_points.append(registered[ordered[idx]])
        idx -= direction

    if len(trend_points) < 2:
        return None  # not enough data to extrapolate

    # average step vector between consecutive trend points (in the direction we walked)
    steps = []
    for i in range(len(trend_points) - 1):
        p_later, p_earlier = trend_points[i], trend_points[i + 1]
        steps.append((p_later[0] - p_earlier[0], p_later[1] - p_earlier[1]))
    avg_dx = sum(s[0] for s in steps) / len(steps)
    avg_dy = sum(s[1] for s in steps) / len(steps)

    anchor = registered[ordered[from_idx]]
    return (anchor[0] + avg_dx, anchor[1] + avg_dy)


def rotate_point(x, y, pivot_x, pivot_y, angle_deg):
    """Rotate (x, y) by angle_deg (counter-clockwise, standard math convention) around a pivot point."""
    angle_rad = math.radians(angle_deg)
    dx, dy = x - pivot_x, y - pivot_y
    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
    return (pivot_x + dx * cos_a - dy * sin_a, pivot_y + dx * sin_a + dy * cos_a)


def ask_user_to_confirm_shift(seg_tiles, left_offset, right_offset, offset, disagreement, registered, ordered, start_idx, end_idx, other_unresolved_tiles=None, context_n=None, parent=None):
    """
    Show an actual plot: the surrounding tiles for context (gray), the
    floating block's CURRENT position (orange x's, dashed outline), and
    where it would move TO if shifted (green circles, solid outline) --
    connected by arrows so the proposed move is visible, not just numbers.

    Two levels of control:
      - Whole-block offset: moves every tile in the block together (starts
        pre-filled with the computed average).
      - Per-tile nudge: pick one specific tile from the block and adjust
        it individually on TOP of the block offset, for fine-tuning a
        single tile that needs a bit more/less than the rest.

    Returns a dict {tile: (dx, dy)} with the final approved offset for
    every tile in the block, or None if skipped entirely.
    Falls back to a plain text version if matplotlib/tkinter aren't available.
    """
    tiles_desc = seg_tiles[0] if len(seg_tiles) == 1 else f"{seg_tiles[0]} .. {seg_tiles[-1]} ({len(seg_tiles)} tiles)"

    if not (MATPLOTLIB_AVAILABLE and TKINTER_AVAILABLE):
        print(
            f"Block: {tiles_desc}\n"
            f"Left estimate: ({left_offset[0]:.1f}, {left_offset[1]:.1f})  "
            f"Right estimate: ({right_offset[0]:.1f}, {right_offset[1]:.1f})  "
            f"Disagree by {disagreement:.1f}px\n"
            f"Computed shift: ({offset[0]:.1f}, {offset[1]:.1f})\n"
            f"[No display/matplotlib available -- can't show a plot]"
        )
        answer = input("Apply this shift? [y]es / [n]o / or type your own whole-block 'dx,dy': ").strip().lower()
        if answer in ("n", "no", ""):
            return None
        if answer in ("y", "yes"):
            block_offset = offset
        else:
            try:
                dx_str, dy_str = answer.split(",")
                block_offset = (float(dx_str.strip()), float(dy_str.strip()))
            except ValueError:
                print("Couldn't parse that as 'dx,dy' -- treating as No.")
                return None

        rot_answer = input("Rotate the whole block too? Enter degrees (counter-clockwise), or blank for 0: ").strip()
        block_rotation = 0.0
        if rot_answer:
            try:
                block_rotation = float(rot_answer)
            except ValueError:
                print("Couldn't parse that as a number -- using 0 degrees.")

        pivot_x = sum(registered[t][0] for t in seg_tiles) / len(seg_tiles)
        pivot_y = sum(registered[t][1] for t in seg_tiles) / len(seg_tiles)

        per_tile = {}
        for t in seg_tiles:
            ox, oy = registered[t]
            rx, ry = rotate_point(ox, oy, pivot_x, pivot_y, block_rotation)
            fx, fy = rx + block_offset[0], ry + block_offset[1]
            per_tile[t] = (fx - ox, fy - oy)

        print(f"Whole-block offset ({block_offset[0]:.1f}, {block_offset[1]:.1f}) + rotation {block_rotation:.1f} deg set. You can now nudge individual tiles on top of this (linear dx,dy AND/OR their own extra rotation around the same block pivot).")
        while True:
            line = input("Tile to nudge (name), or blank to finish: ").strip()
            if not line:
                break
            if line not in per_tile:
                print(f"'{line}' isn't in this block ({seg_tiles}). Try again.")
                continue
            try:
                extra = input(f"Extra dx,dy to add to {line}'s current offset {per_tile[line]}: ").strip()
                edx, edy = (float(v.strip()) for v in extra.split(","))
                extra_rot_str = input(f"Extra rotation (deg) for {line} around the block pivot, or blank for 0: ").strip()
                extra_rot = float(extra_rot_str) if extra_rot_str else 0.0

                ox, oy = registered[line]
                rx, ry = rotate_point(ox, oy, pivot_x, pivot_y, block_rotation + extra_rot)
                fx, fy = rx + block_offset[0] + edx, ry + block_offset[1] + edy
                per_tile[line] = (fx - ox, fy - oy)
                print(f"{line} is now at offset {per_tile[line]}")
            except ValueError:
                print("Couldn't parse that, skipped.")
        return per_tile

    # Show the WHOLE sequence for full context, not just a local window --
    # everything before the block, everything after it. Other tiles that
    # belong to STILL-UNRESOLVED candidate blocks get pulled out and shown
    # in a distinct color, so you can see all the outstanding problems at
    # once instead of them blending into the "trusted" context.
    other_unresolved_tiles = set(other_unresolved_tiles or ())
    context_before_all = ordered[0:start_idx]
    context_after_all = ordered[end_idx + 1:]
    context_before = [t for t in context_before_all if t not in other_unresolved_tiles]
    context_after = [t for t in context_after_all if t not in other_unresolved_tiles]
    other_flagged = [t for t in (context_before_all + context_after_all) if t in other_unresolved_tiles]
    block_tiles = ordered[start_idx:end_idx + 1]

    fig, ax = plt.subplots(figsize=(9, 7))

    def xy(tiles):
        return [registered[t][0] for t in tiles], [registered[t][1] for t in tiles]

    bx, by = xy(block_tiles)

    block_offset = [offset[0], offset[1]]
    individual_deltas = {t: [0.0, 0.0] for t in block_tiles}
    block_rotation = [0.0]
    individual_rotations = {t: 0.0 for t in block_tiles}

    pivot_x = sum(bx) / len(bx)
    pivot_y = sum(by) / len(by)

    def final_positions():
        result = {}
        for t, x, y in zip(block_tiles, bx, by):
            total_rotation = block_rotation[0] + individual_rotations[t]
            rx, ry = rotate_point(x, y, pivot_x, pivot_y, total_rotation)
            fx = rx + block_offset[0] + individual_deltas[t][0]
            fy = ry + block_offset[1] + individual_deltas[t][1]
            result[t] = (fx, fy)
        return result

    def final_offsets():
        fp = final_positions()
        return {t: (fp[t][0] - x, fp[t][1] - y) for t, x, y in zip(block_tiles, bx, by)}

    LABEL_LIMIT = 30  # avoid cluttering the plot with hundreds of labels on large context sets

    def render():
        ax.clear()
        if context_before:
            cx, cy = xy(context_before)
            ax.plot(cx, cy, "o-", color="steelblue", label="Context (before block)", markersize=6)
            if len(context_before) <= LABEL_LIMIT:
                for t, x, y in zip(context_before, cx, cy):
                    ax.annotate(t, (x, y), fontsize=6, color="steelblue")
            else:
                # still label the endpoints so you know where you are in the sequence
                ax.annotate(context_before[0], (cx[0], cy[0]), fontsize=7, color="steelblue", fontweight="bold")
                ax.annotate(context_before[-1], (cx[-1], cy[-1]), fontsize=7, color="steelblue", fontweight="bold")
        if context_after:
            cx, cy = xy(context_after)
            ax.plot(cx, cy, "o-", color="seagreen", label="Context (after block)", markersize=6)
            if len(context_after) <= LABEL_LIMIT:
                for t, x, y in zip(context_after, cx, cy):
                    ax.annotate(t, (x, y), fontsize=6, color="seagreen")
            else:
                ax.annotate(context_after[0], (cx[0], cy[0]), fontsize=7, color="seagreen", fontweight="bold")
                ax.annotate(context_after[-1], (cx[-1], cy[-1]), fontsize=7, color="seagreen", fontweight="bold")

        if other_flagged:
            ox, oy = xy(other_flagged)
            ax.scatter(ox, oy, c="red", marker="X", s=90, label="OTHER unresolved block(s)", zorder=6, edgecolors="black")
            for t, x, y in zip(other_flagged, ox, oy):
                ax.annotate(t, (x, y), fontsize=6, color="red", fontweight="bold")

        ax.plot(bx, by, "x--", color="darkorange", label="Block -- CURRENT position", markersize=10, markeredgewidth=2)
        for t, x, y in zip(block_tiles, bx, by):
            ax.annotate(t, (x, y), fontsize=6, color="darkorange")

        fp = final_positions()
        new_bx = [fp[t][0] for t in block_tiles]
        new_by = [fp[t][1] for t in block_tiles]
        colors = ["gold" if target_var.get() == t else "limegreen" for t in block_tiles]
        ax.scatter(new_bx, new_by, c=colors, edgecolors="black", s=80, zorder=5, label="Block -- PROPOSED position")
        ax.plot(new_bx, new_by, "-", color="limegreen", alpha=0.5, zorder=4)
        for x1, y1, x2, y2 in zip(bx, by, new_bx, new_by):
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", color="gray", alpha=0.6))
        ax.scatter([pivot_x], [pivot_y], marker="+", c="black", s=150, linewidths=2, label="Rotation pivot (block centroid)", zorder=7)

        target = target_var.get()
        rot_display = block_rotation[0] if target == "Whole Block" else individual_rotations[target]
        ax.set_title(f"Proposed shift for {tiles_desc}\nBlock offset: ({block_offset[0]:.1f}, {block_offset[1]:.1f})px  |  Block rotation: {block_rotation[0]:.1f}\u00b0  |  Moving: {target} (rotation: {rot_display:.1f}\u00b0)  |  L/R disagreement: {disagreement:.1f}px")
        ax.set_aspect("equal", adjustable="datalim")
        ax.invert_yaxis()
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        canvas.draw()

    result = {"offsets": None}

    # IMPORTANT: if a parent Tk root is already running (e.g. the block
    # selector window), create a Toplevel bound to it instead of a second
    # independent tk.Tk(). Two separate Tk() instances running nested
    # inside each other's callbacks is unsupported and causes exactly the
    # kind of stale/non-refreshing rendering this was fixing.
    root = tk.Toplevel(parent) if parent is not None else tk.Tk()
    root.title("Confirm block shift")

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    control_frame = tk.Frame(root)
    control_frame.pack(side=tk.TOP, fill=tk.X, pady=4)

    tk.Label(control_frame, text="Move:").pack(side=tk.LEFT, padx=(10, 2))
    target_var = tk.StringVar(value="Whole Block")
    target_menu = tk.OptionMenu(control_frame, target_var, "Whole Block", *block_tiles, command=lambda _=None: refresh_entries_and_render())
    target_menu.pack(side=tk.LEFT)

    tk.Label(control_frame, text="X:").pack(side=tk.LEFT, padx=(10, 2))
    x_var = tk.StringVar()
    x_entry = tk.Entry(control_frame, textvariable=x_var, width=10)
    x_entry.pack(side=tk.LEFT)
    tk.Label(control_frame, text="Y:").pack(side=tk.LEFT, padx=(10, 2))
    y_var = tk.StringVar()
    y_entry = tk.Entry(control_frame, textvariable=y_var, width=10)
    y_entry.pack(side=tk.LEFT)
    tk.Label(control_frame, text="Rot(\u00b0):").pack(side=tk.LEFT, padx=(10, 2))
    rot_var = tk.StringVar()
    rot_entry = tk.Entry(control_frame, textvariable=rot_var, width=8)
    rot_entry.pack(side=tk.LEFT)

    tk.Label(control_frame, text="Step:").pack(side=tk.LEFT, padx=(10, 2))
    step_var = tk.StringVar(value="10")
    tk.Entry(control_frame, textvariable=step_var, width=6).pack(side=tk.LEFT)
    tk.Label(control_frame, text="Rot Step(\u00b0):").pack(side=tk.LEFT, padx=(10, 2))
    rot_step_var = tk.StringVar(value="5")
    tk.Entry(control_frame, textvariable=rot_step_var, width=5).pack(side=tk.LEFT)

    def get_target_value():
        target = target_var.get()
        if target == "Whole Block":
            return block_offset[0], block_offset[1], block_rotation[0]
        return individual_deltas[target][0], individual_deltas[target][1], individual_rotations[target]

    def set_target_value(vx, vy, vrot):
        target = target_var.get()
        if target == "Whole Block":
            block_offset[0], block_offset[1] = vx, vy
            block_rotation[0] = vrot
        else:
            individual_deltas[target][0], individual_deltas[target][1] = vx, vy
            individual_rotations[target] = vrot

    def refresh_entries_and_render():
        vx, vy, vrot = get_target_value()
        x_var.set(f"{vx:.2f}")
        y_var.set(f"{vy:.2f}")
        rot_var.set(f"{vrot:.2f}")
        render()

    def on_apply_typed_values():
        try:
            _, _, cur_rot = get_target_value()
            new_rot = float(rot_var.get()) if rot_var.get().strip() else cur_rot
            set_target_value(float(x_var.get()), float(y_var.get()), new_rot)
        except ValueError:
            pass
        render()

    tk.Button(control_frame, text="Preview", command=on_apply_typed_values).pack(side=tk.LEFT, padx=10)

    def nudge(dx_dir, dy_dir):
        try:
            step = float(step_var.get())
        except ValueError:
            step = 10.0
        vx, vy, vrot = get_target_value()
        set_target_value(vx + dx_dir * step, vy + dy_dir * step, vrot)
        refresh_entries_and_render()

    def nudge_rotation(direction):
        try:
            rot_step = float(rot_step_var.get())
        except ValueError:
            rot_step = 5.0
        vx, vy, vrot = get_target_value()
        set_target_value(vx, vy, vrot + direction * rot_step)
        refresh_entries_and_render()

    arrow_frame = tk.Frame(root)
    arrow_frame.pack(side=tk.TOP, pady=4)
    tk.Button(arrow_frame, text="\u2191", command=lambda: nudge(0, -1), width=4).grid(row=0, column=1)
    tk.Button(arrow_frame, text="\u2190", command=lambda: nudge(-1, 0), width=4).grid(row=1, column=0)
    tk.Button(arrow_frame, text="\u2192", command=lambda: nudge(1, 0), width=4).grid(row=1, column=2)
    tk.Button(arrow_frame, text="\u2193", command=lambda: nudge(0, 1), width=4).grid(row=2, column=1)
    tk.Button(arrow_frame, text="\u21ba", command=lambda: nudge_rotation(-1), width=4).grid(row=0, column=3, padx=(15, 0))
    tk.Button(arrow_frame, text="\u21bb", command=lambda: nudge_rotation(1), width=4).grid(row=0, column=4)

    button_frame = tk.Frame(root)
    button_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=8)

    def on_yes():
        result["offsets"] = final_offsets()
        plt.close(fig)
        root.destroy()

    def on_no():
        result["offsets"] = None
        plt.close(fig)
        root.destroy()

    tk.Button(button_frame, text="Apply these offsets", command=on_yes, bg="#c8f7c5", width=20, height=2).pack(side=tk.LEFT, expand=True, padx=10)
    tk.Button(button_frame, text="Skip (leave as-is)", command=on_no, bg="#f7c5c5", width=20, height=2).pack(side=tk.RIGHT, expand=True, padx=10)

    root.protocol("WM_DELETE_WINDOW", on_no)
    refresh_entries_and_render()

    if parent is not None:
        root.transient(parent)
        root.grab_set()
        parent.wait_window(root)
    else:
        root.mainloop()

    return result["offsets"]


def compute_offset_guess(left_offset, right_offset):
    if left_offset is not None and right_offset is not None:
        return ((left_offset[0] + right_offset[0]) / 2, (left_offset[1] + right_offset[1]) / 2)
    return left_offset or right_offset or (0.0, 0.0)


def select_and_review_blocks(candidates, registered, ordered):
    """
    candidates: list of dicts with keys:
        seg_tiles, start, end, left_offset, right_offset, disagreement (or None),
        kind ('disagreement'|'one_sided'), prev_end_idx, next_start_idx

    Shows a menu of all blocks that weren't auto-corrected, lets the user
    pick ANY of them (in any order, can review more than one, can skip
    entirely). Maintains a WORKING copy of positions that updates after
    each applied fix, so reviewing one block reflects any earlier fixes
    (recomputes that block's left/right estimates fresh each time, rather
    than using stale numbers from before anything was corrected). Other
    still-unresolved blocks are shown highlighted in the plot too, so you
    can see every outstanding problem at once, not just the current one.

    Returns (results, working_registered): results = {index: per_tile_offsets_dict_or_None},
    working_registered = final positions after all applied fixes.
    """
    results = {i: None for i in range(len(candidates))}
    working_registered = dict(registered)

    def recompute_offsets(c):
        left_offset = None
        if c["prev_end_idx"] is not None:
            predicted_start_pos = estimate_extrapolated_position(ordered, working_registered, c["prev_end_idx"], direction=+1)
            if predicted_start_pos:
                actual_start_pos = working_registered[ordered[c["start"]]]
                left_offset = (predicted_start_pos[0] - actual_start_pos[0], predicted_start_pos[1] - actual_start_pos[1])
        right_offset = None
        if c["next_start_idx"] is not None:
            predicted_end_pos = estimate_extrapolated_position(ordered, working_registered, c["next_start_idx"], direction=-1)
            if predicted_end_pos:
                actual_end_pos = working_registered[ordered[c["end"]]]
                right_offset = (predicted_end_pos[0] - actual_end_pos[0], predicted_end_pos[1] - actual_end_pos[1])
        return left_offset, right_offset

    def other_unresolved_tiles(exclude_idx):
        tiles = set()
        for j, c in enumerate(candidates):
            if j != exclude_idx and results[j] is None:
                tiles.update(c["seg_tiles"])
        return tiles

    def describe(i):
        c = candidates[i]
        seg_tiles = c["seg_tiles"]
        tiles_desc = seg_tiles[0] if len(seg_tiles) == 1 else f"{seg_tiles[0]} .. {seg_tiles[-1]} ({len(seg_tiles)} tiles)"
        status = "[FIXED]" if results[i] is not None else "[unresolved]"
        if c["kind"] == "disagreement":
            detail = f"disagreement {c['disagreement']:.1f}px"
        else:
            side = "left" if c["left_offset"] is not None else "right"
            detail = f"one-sided ({side} only)"
        return f"{tiles_desc} -- {detail} {status}"

    def review(idx, parent=None):
        c = candidates[idx]
        left_offset, right_offset = recompute_offsets(c)
        disagreement = dist(left_offset, right_offset) if (left_offset is not None and right_offset is not None) else (c["disagreement"] or 0.0)
        offset_guess = compute_offset_guess(left_offset, right_offset)
        per_tile = ask_user_to_confirm_shift(
            c["seg_tiles"], left_offset or (0.0, 0.0), right_offset or (0.0, 0.0),
            offset_guess, disagreement, working_registered, ordered, c["start"], c["end"],
            other_unresolved_tiles=other_unresolved_tiles(idx), parent=parent,
        )
        results[idx] = per_tile
        if per_tile is not None:
            for t, (dx, dy) in per_tile.items():
                old_pos = registered[t]
                working_registered[t] = (old_pos[0] + dx, old_pos[1] + dy)

    if not (MATPLOTLIB_AVAILABLE and TKINTER_AVAILABLE):
        while True:
            print("\n=== Blocks available for review ===")
            for i in range(len(candidates)):
                print(f"  [{i}] {describe(i)}")
            choice = input("Enter a number to review, or blank to finish: ").strip()
            if not choice:
                break
            try:
                idx = int(choice)
                candidates[idx]
            except (ValueError, IndexError):
                print("Not a valid number, try again.")
                continue
            review(idx)
        return results, working_registered

    root = tk.Tk()
    root.title("Select a block to review")
    tk.Label(root, text="Blocks that weren't auto-corrected -- pick one to review:", font=("", 11)).pack(padx=10, pady=(10, 4))

    listbox = tk.Listbox(root, width=70, height=min(15, max(4, len(candidates))))
    listbox.pack(padx=10, pady=4, fill=tk.BOTH, expand=True)

    def refresh_listbox():
        listbox.delete(0, tk.END)
        for i in range(len(candidates)):
            listbox.insert(tk.END, describe(i))

    def on_review():
        sel = listbox.curselection()
        if not sel:
            return
        review(sel[0], parent=root)
        refresh_listbox()

    button_frame = tk.Frame(root)
    button_frame.pack(pady=8)
    tk.Button(button_frame, text="Review selected block", command=on_review, width=22, height=2).pack(side=tk.LEFT, padx=10)
    tk.Button(button_frame, text="Done", command=root.destroy, width=15, height=2, bg="#c8f7c5").pack(side=tk.LEFT, padx=10)

    refresh_listbox()
    root.mainloop()

    return results, working_registered


def identify_and_fix(registered, expect_loop=False, jump_threshold_multiplier=3.0, offset_agreement_tolerance_multiplier=1.0, interactive=False):
    ordered = build_ordered_tiles(registered, expect_loop)
    n = len(ordered)
    if n < 3:
        print("Not enough tiles to run this repair.")
        return dict(registered), []

    flags, threshold, median_step = find_jump_edges(ordered, registered, jump_threshold_multiplier)
    print(f"Typical step distance: {median_step:.1f}px, flagging jumps over {threshold:.1f}px")
    print(f"Found {sum(flags)} flagged step(s) out of {len(flags)}.")

    segments = split_into_segments(flags)
    print(f"Split into {len(segments)} segment(s) (contiguous runs with no internal jumps).")

    corrected = dict(registered)
    fixed_report = []
    candidates = []  # blocks needing manual review (disagreement or one-sided)
    agreement_tolerance = median_step * offset_agreement_tolerance_multiplier

    for seg_i, (start, end) in enumerate(segments):
        has_prev = seg_i > 0
        has_next = seg_i < len(segments) - 1

        left_offset = None
        if has_prev:
            prev_end = segments[seg_i - 1][1]
            predicted_start_pos = estimate_extrapolated_position(ordered, registered, prev_end, direction=+1)
            if predicted_start_pos:
                actual_start_pos = registered[ordered[start]]
                left_offset = (predicted_start_pos[0] - actual_start_pos[0], predicted_start_pos[1] - actual_start_pos[1])

        right_offset = None
        if has_next:
            next_start = segments[seg_i + 1][0]
            predicted_end_pos = estimate_extrapolated_position(ordered, registered, next_start, direction=-1)
            if predicted_end_pos:
                actual_end_pos = registered[ordered[end]]
                right_offset = (predicted_end_pos[0] - actual_end_pos[0], predicted_end_pos[1] - actual_end_pos[1])

        seg_tiles = [ordered[i] for i in range(start, end + 1)]

        if left_offset is None and right_offset is None:
            continue  # this is the only segment (whole trace is one block) -- nothing to align against

        if left_offset is not None and right_offset is not None:
            disagreement = dist(left_offset, right_offset)
            offset = ((left_offset[0] + right_offset[0]) / 2, (left_offset[1] + right_offset[1]) / 2)
            if disagreement > agreement_tolerance:
                candidates.append(dict(seg_tiles=seg_tiles, start=start, end=end, left_offset=left_offset, right_offset=right_offset, disagreement=disagreement, kind="disagreement", prev_end_idx=(segments[seg_i - 1][1] if has_prev else None), next_start_idx=(segments[seg_i + 1][0] if has_next else None)))
                continue
            confidence = "high (both sides agree)"
        else:
            # Only one neighboring segment exists (this is the very first or
            # very last segment of the whole trace) -- there's no independent
            # way to verify which side is actually wrong by itself, so it
            # goes to manual review instead of being auto-corrected.
            candidates.append(dict(seg_tiles=seg_tiles, start=start, end=end, left_offset=left_offset, right_offset=right_offset, disagreement=None, kind="one_sided", prev_end_idx=(segments[seg_i - 1][1] if has_prev else None), next_start_idx=(segments[seg_i + 1][0] if has_next else None)))
            continue

        if abs(offset[0]) < 1e-6 and abs(offset[1]) < 1e-6:
            continue  # nothing to correct, already aligned

        per_tile_offsets = {t: offset for t in seg_tiles}
        for t in seg_tiles:
            old_pos = registered[t]
            corrected[t] = (old_pos[0] + offset[0], old_pos[1] + offset[1])
        fixed_report.append((seg_tiles, per_tile_offsets, confidence))

    ambiguous_report = []

    if candidates and interactive:
        print(f"\n{len(candidates)} block(s) need manual review -- opening block selector...")
        review_results, working_registered = select_and_review_blocks(candidates, corrected, ordered)
        for i, c in enumerate(candidates):
            per_tile_offsets = review_results.get(i)
            if per_tile_offsets is None:
                ambiguous_report.append((c["seg_tiles"], c["left_offset"], c["right_offset"], c["disagreement"]))
                continue
            for t in c["seg_tiles"]:
                corrected[t] = working_registered[t]
            distinct = set(per_tile_offsets.values())
            if len(distinct) == 1:
                confidence = "user-reviewed via block selector"
            else:
                confidence = "user-reviewed via block selector, with per-tile nudges"
            fixed_report.append((c["seg_tiles"], per_tile_offsets, confidence))
    else:
        for c in candidates:
            ambiguous_report.append((c["seg_tiles"], c["left_offset"], c["right_offset"], c["disagreement"]))

    if fixed_report:
        print(f"\nREALIGNED {len(fixed_report)} block(s):")
        for seg_tiles, per_tile_offsets, confidence in fixed_report:
            tiles_desc = seg_tiles[0] if len(seg_tiles) == 1 else f"{seg_tiles[0]} .. {seg_tiles[-1]} ({len(seg_tiles)} tiles)"
            distinct = set(per_tile_offsets.values())
            if len(distinct) == 1:
                off = next(iter(distinct))
                print(f"    {tiles_desc}: shifted by ({off[0]:.1f}, {off[1]:.1f})  -- confidence: {confidence}")
            else:
                print(f"    {tiles_desc}: per-tile offsets applied -- confidence: {confidence}")
                for t in seg_tiles:
                    dx, dy = per_tile_offsets[t]
                    print(f"        {t}: ({dx:.1f}, {dy:.1f})")
    else:
        print("\nNo blocks were realigned.")

    if ambiguous_report:
        print(f"\n{len(ambiguous_report)} block(s) were NOT auto-corrected (insufficient independent evidence):")
        for seg_tiles, left_offset, right_offset, disagreement in ambiguous_report:
            tiles_desc = seg_tiles[0] if len(seg_tiles) == 1 else f"{seg_tiles[0]} .. {seg_tiles[-1]} ({len(seg_tiles)} tiles)"
            if disagreement is not None:
                print(f"    {tiles_desc}: left-side estimate {left_offset}, right-side estimate {right_offset} DISAGREE by {disagreement:.1f}px -- likely needs rotation, not just translation.")
            else:
                side = "left" if left_offset is not None else "right"
                only_offset = left_offset if left_offset is not None else right_offset
                print(f"    {tiles_desc}: only a {side}-side estimate available ({only_offset}) -- this is an end segment, can't cross-check, review manually.")

    return corrected, fixed_report


def write_tile_configuration(positions, out_path):
    with open(out_path, "w") as f:
        f.write("# Define the number of dimensions we are working on\n")
        f.write("dim = 2\n\n")
        f.write("# Define the image coordinates\n")
        for tile in sorted(positions):
            x, y = positions[tile]
            f.write(f"{tile}; ; ({x:.4f}, {y:.4f})\n")
    print(f"\nSaved corrected positions to {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--tile-path", help="Full path to ONE tile file inside the sequence folder -- auto-derives the registered file's location.")
    mode.add_argument("--sequence-dir", help="Full path to the sequence folder itself (e.g. '/path/to/project/slices/MySeq') -- auto-finds TileConfiguration.registered.txt inside it.")
    mode.add_argument("--registered", help="Direct path to TileConfiguration.registered.txt")

    parser.add_argument("--out", default="TileConfiguration.jumpfixed.txt", help="Output filename/path for the corrected position file")
    parser.add_argument("--loop", action="store_true", help="Structure is a closed loop (affects nothing here yet, reserved for future circular handling)")
    parser.add_argument("--jump-threshold-multiplier", type=float, default=3.0, help="A step is 'flagged' if it's more than this many times the typical step distance (default 3.0)")
    parser.add_argument("--interactive", action="store_true", help="When the left-side and right-side offset estimates disagree, pop up a window showing the proposed shift and ask for confirmation before applying it, instead of skipping automatically.")
    args = parser.parse_args()

    if args.tile_path:
        tile_path = os.path.normpath(args.tile_path)
        sequence_dir = os.path.dirname(tile_path)
        registered_path = find_registered_file(sequence_dir)
        out_path = os.path.join(sequence_dir, args.out) if not os.path.isabs(args.out) else args.out
    elif args.sequence_dir:
        sequence_dir = os.path.normpath(args.sequence_dir)
        if not os.path.isdir(sequence_dir):
            print(f"--sequence-dir path does not exist or isn't a folder: {sequence_dir}", file=sys.stderr)
            sys.exit(1)
        registered_path = find_registered_file(sequence_dir)
        out_path = os.path.join(sequence_dir, args.out) if not os.path.isabs(args.out) else args.out
    else:
        registered_path = args.registered
        out_path = args.out

    if not registered_path:
        print("No TileConfiguration.registered.txt found.", file=sys.stderr)
        sys.exit(1)

    registered = parse_tile_configuration(registered_path)
    print(f"Loaded {len(registered)} tiles from {registered_path}\n")

    corrected, fixed_report = identify_and_fix(registered, expect_loop=args.loop, jump_threshold_multiplier=args.jump_threshold_multiplier, interactive=args.interactive)
    write_tile_configuration(corrected, out_path)


if __name__ == "__main__":
    main()
