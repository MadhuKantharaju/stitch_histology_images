"""
topology_graph.py

Draws an interactive graph of your tiles at their ACTUAL registered
positions, with edges colored by trust status:
  - green solid  = trusted, geometrically-consistent link
  - orange dashed = high correlation but contradicts the final registration
                     (likely a false-positive match, e.g. repetitive tissue)
  - gray dotted  = rejected (below correlation threshold)

Useful for irregular/curved layouts (e.g. anatomical structures) where you
expect a specific topology (a closed loop, say) and want to see directly
whether nodes have collapsed on top of each other or the shape is tangled.

Usage (simplest -- point at one tile file, everything else is derived):
    python topology_graph.py --tile-path "/path/to/project/slices/MySeq/Snap-1.czi"

Or manually:
    python topology_graph.py --log-dir "/path/to/logs" --registered "/path/to/TileConfiguration.registered.txt" --out topology.html

Requires: pip install plotly
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import plotly.graph_objects as go

IMG_EXT_PATTERN = r"[\w\-.]+\.(?:tif|tiff|czi|png|jpg|jpeg)"

# Format: "image1.tif <- image2.tif: (12.3, 456.7) correlation (R)=0.891"
SHIFT_PATTERN = re.compile(
    r"(?P<file1>" + IMG_EXT_PATTERN + r")(?:\[\d+\])?\s*(?:<-|->)\s*"
    r"(?P<file2>" + IMG_EXT_PATTERN + r")(?:\[\d+\])?\s*:\s*"
    r"\(\s*(?P<dx>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*,\s*(?P<dy>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\)"
    r".*?correlation\s*\(?R?\)?\s*[:=]?\s*(?P<r>-?\d+\.?\d*(?:[eE][+-]?\d+)?)",
    re.IGNORECASE,
)

CORR_PATTERN = re.compile(
    r"(?P<file1>" + IMG_EXT_PATTERN + r")"
    r".*?"
    r"(?P<file2>" + IMG_EXT_PATTERN + r")"
    r".*?"
    r"(?:correlation\s*\(?R\)?\s*[:=]|R\s*=)\s*"
    r"(?P<r>-?\d+\.?\d*(?:[eE][+-]?\d+)?)",
    re.IGNORECASE,
)

INVALID_PATTERN = re.compile(r"(invalid|below threshold|rejected|no overlap)", re.IGNORECASE)


def parse_log_file(path):
    """Return a list of dicts: {file1, file2, correlation, dx, dy, valid}."""
    records = []
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()

    for line in lines:
        files, r_val, dx, dy = None, None, None, None

        m_shift = SHIFT_PATTERN.search(line)
        if m_shift:
            files = (m_shift.group("file1").strip(), m_shift.group("file2").strip())
            r_val = float(m_shift.group("r"))
            dx = float(m_shift.group("dx"))
            dy = float(m_shift.group("dy"))
        else:
            m = CORR_PATTERN.search(line)
            if m:
                files = (m.group("file1").strip(), m.group("file2").strip())
                r_val = float(m.group("r"))

        if r_val is not None and files is not None:
            valid = not bool(INVALID_PATTERN.search(line))
            records.append(
                {
                    "file1": os.path.basename(files[0]),
                    "file2": os.path.basename(files[1]),
                    "correlation": r_val,
                    "dx": dx,
                    "dy": dy,
                    "valid": valid,
                }
            )
    return records


def parse_tile_configuration(path):
    """Return dict: filename -> (x, y)."""
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
    matches = glob.glob(os.path.join(seq_dir, "**", "TileConfiguration.jumpfixed.txt"), recursive=True)
    return matches[0] if matches else None


def analyze_clusters(pair_records, registered, shift_discrepancy_threshold=50.0):
    """Find clusters of tiles tied together by trusted-AND-consistent links."""
    all_tiles = set()
    adjacency = defaultdict(set)
    all_links = defaultdict(list)
    suspicious_links = []

    for rec in pair_records:
        t1, t2 = rec["file1"], rec["file2"]
        all_tiles.add(t1)
        all_tiles.add(t2)

        geometrically_consistent = True
        dx, dy = rec.get("dx"), rec.get("dy")
        if dx is not None and dy is not None:
            pos1, pos2 = registered.get(t1), registered.get(t2)
            if pos1 and pos2:
                # Compare the full displacement VECTOR, not just its magnitude --
                # two points can be the same DISTANCE apart in a totally wrong
                # DIRECTION, which a magnitude-only check would wrongly accept.
                actual_dx = pos2[0] - pos1[0]
                actual_dy = pos2[1] - pos1[1]
                vector_error = ((actual_dx - dx) ** 2 + (actual_dy - dy) ** 2) ** 0.5
                if vector_error > shift_discrepancy_threshold:
                    geometrically_consistent = False

        trusted_final = rec["valid"] and geometrically_consistent
        all_links[t1].append((t2, trusted_final, rec["correlation"]))
        all_links[t2].append((t1, trusted_final, rec["correlation"]))

        if trusted_final:
            adjacency[t1].add(t2)
            adjacency[t2].add(t1)
        elif rec["valid"] and not geometrically_consistent:
            suspicious_links.append(rec)

    all_tiles |= set(registered)

    if suspicious_links:
        print(
            f"\n{len(suspicious_links)} link(s) had high correlation but disagreed with the "
            f"final registration by >{shift_discrepancy_threshold}px -- likely false-positive "
            f"matches (e.g. repetitive tissue). Treated as untrustworthy:"
        )
        for rec in suspicious_links:
            print(f"    {rec['file1']} <-> {rec['file2']}: correlation={rec['correlation']:.3f} (contradicted)")

    tile_cluster = {}
    cluster_id = 0
    for tile in sorted(all_tiles):
        if tile in tile_cluster:
            continue
        stack = [tile]
        tile_cluster[tile] = cluster_id
        while stack:
            cur = stack.pop()
            for nb in adjacency.get(cur, ()):
                if nb not in tile_cluster:
                    tile_cluster[nb] = cluster_id
                    stack.append(nb)
        cluster_id += 1

    n_clusters = cluster_id
    sizes = defaultdict(int)
    for c in tile_cluster.values():
        sizes[c] += 1
    print(f"\nFound {n_clusters} cluster(s) of tiles connected by trusted links.")
    for c in sorted(sizes, key=lambda k: -sizes[k]):
        print(f"  Cluster {c}: {sizes[c]} tiles")

    cross_links = defaultdict(list)
    for tile, links in all_links.items():
        c1 = tile_cluster.get(tile)
        for neighbor, valid, corr in links:
            c2 = tile_cluster.get(neighbor)
            if c1 is None or c2 is None or c1 == c2:
                continue
            key = tuple(sorted((c1, c2)))
            cross_links[key].append((valid, corr))

    if cross_links:
        print("\n--- How each cluster connects to the others ---")
        for (ca, cb), links in sorted(cross_links.items()):
            n_total = len(links)
            n_valid = sum(1 for v, _ in links if v)
            avg_corr = sum(c for _, c in links) / n_total
            verdict = ""
            if n_valid == 0:
                verdict = "  <-- ALL links REJECTED. Relative position is essentially a guess."
            elif n_valid <= 2:
                verdict = f"  <-- only {n_valid} trusted link(s) -- fragile anchor."
            print(f"  Cluster {ca} <-> Cluster {cb}: {n_total} link(s), {n_valid} trusted, avg correlation {avg_corr:.3f}{verdict}")

    return tile_cluster


def analyze_spatial_proximity_violations(registered, max_scan_jump=3, spatial_threshold_multiplier=1.5, corroboration_window=3, expect_loop=False):
    """
    Pure geometry check, no log/correlation data needed at all: look at
    where every tile ACTUALLY ended up (x, y) and flag any pair that landed
    spatially close together (neighbors in ANY direction) despite not being
    expected to be neighbors based on scan order. This catches problems
    even between tiles Fiji never directly compared in the log.

    'Expected neighbor' = scan-order index difference <= max_scan_jump
    (circular distance if expect_loop=True, so the legitimate ring-closing
    seam isn't wrongly flagged). 'Spatially close' = distance <= (typical
    nearest-neighbor spacing) * spatial_threshold_multiplier.
    """
    tiles = sorted(registered)
    n = len(tiles)

    def dist(a, b):
        return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5

    all_indices = [extract_scan_index(t) for t in tiles]
    all_indices = [i for i in all_indices if i is not None]
    idx_range = (max(all_indices) - min(all_indices)) if all_indices else None

    # Estimate typical tile spacing from each tile's nearest neighbor
    nn_dists = []
    for i, t in enumerate(tiles):
        best = None
        for j, u in enumerate(tiles):
            if i == j:
                continue
            d = dist(registered[t], registered[u])
            if best is None or d < best:
                best = d
        if best is not None:
            nn_dists.append(best)

    if not nn_dists:
        print("\nNot enough tiles to estimate spacing for the spatial-proximity check.")
        return {"isolated": [], "corroborated": [], "median_spacing": None}

    nn_dists.sort()
    median_spacing = nn_dists[len(nn_dists) // 2]
    spatial_threshold = median_spacing * spatial_threshold_multiplier

    print(f"\n=== Pure positional proximity check (no log/correlation data used) ===")
    print(f"Estimated typical tile spacing: {median_spacing:.1f}px -- flagging pairs closer than {spatial_threshold:.1f}px")

    candidates = []
    for i in range(n):
        for j in range(i + 1, n):
            t1, t2 = tiles[i], tiles[j]
            idx1, idx2 = extract_scan_index(t1), extract_scan_index(t2)
            if idx1 is None or idx2 is None:
                continue
            raw_diff = abs(idx1 - idx2)
            if expect_loop and idx_range:
                idx_diff = min(raw_diff, idx_range - raw_diff)
            else:
                idx_diff = raw_diff
            if idx_diff <= max_scan_jump:
                continue  # expected to be neighbors anyway -- not a candidate
            d = dist(registered[t1], registered[t2])
            if d <= spatial_threshold:
                candidates.append((t1, t2, idx_diff, d))

    if candidates:
        print(f"\nRULE VIOLATION -- {len(candidates)} pair(s) of tiles ended up spatially close (within {spatial_threshold:.1f}px) despite NOT being expected neighbors:")
        for t1, t2, idx_diff, d in sorted(candidates, key=lambda x: x[3])[:20]:
            print(f"    {t1} <-> {t2}  (scan-index gap={idx_diff}, distance={d:.1f}px)")
    else:
        print("\nNo unexpected spatial neighbors found.")

    return {"isolated": candidates, "corroborated": [], "median_spacing": median_spacing}


def segments_intersect(p1, p2, p3, p4):
    """Standard segment-intersection test (excludes shared-endpoint 'touching', only true crossings)."""
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    d1 = cross(p3, p4, p1)
    d2 = cross(p3, p4, p2)
    d3 = cross(p1, p2, p3)
    d4 = cross(p1, p2, p4)

    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    return False


def extract_scan_index(filename):
    """Pull the trailing acquisition number out of a filename like 'Snap-17149.czi' -> 17149."""
    m = re.search(r"(\d+)(?=\.[^.]+$)", filename)
    return int(m.group(1)) if m else None


def build_index_ordered_trace(registered, expect_loop=False):
    """
    Build the trace using FILENAMES (scan order), not blind nearest-distance.
    Tile 19 always connects to tile 20 -- its actual sequential neighbor --
    never to some unrelated tile 36 just because it happens to be closer in
    space. Positions are only used afterward to sanity-check this trace
    (are consecutive tiles actually close together? does it cross itself?),
    never to decide which tiles connect.
    """
    tiles_with_idx = [(t, extract_scan_index(t)) for t in registered]
    tiles_with_idx = [(t, i) for t, i in tiles_with_idx if i is not None]
    tiles_with_idx.sort(key=lambda x: x[1])

    edges = [(tiles_with_idx[k][0], tiles_with_idx[k + 1][0]) for k in range(len(tiles_with_idx) - 1)]
    if expect_loop and len(tiles_with_idx) >= 2:
        edges.append((tiles_with_idx[-1][0], tiles_with_idx[0][0]))

    return edges, [t for t, _ in tiles_with_idx]


def check_index_guided_trace(registered, expect_loop=False, jump_threshold_multiplier=3.0):
    """
    The trace itself comes purely from filename/scan order. Positions are
    only used to sanity-check it: does each consecutive-index step actually
    land close together (no unexpected jump), and does the resulting path
    ever cross itself (which would mean something is geometrically wrong
    even though the scan order itself looks fine)?
    """
    edges, ordered_tiles = build_index_ordered_trace(registered, expect_loop)
    print(f"\n=== Index-guided trace check ({len(edges)} sequential steps, built from filenames -- positions only used to sanity-check) ===")
    if not edges:
        print("Not enough tiles with parseable scan-order numbers to build a trace.")
        return {"traceable": None, "jumps": [], "crossings": []}

    def dist(a, b):
        return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5

    step_dists = []
    usable_edges = []
    for t1, t2 in edges:
        if t1 in registered and t2 in registered:
            step_dists.append(dist(registered[t1], registered[t2]))
            usable_edges.append((t1, t2))

    if not step_dists:
        print("No usable positions to sanity-check the trace.")
        return {"traceable": None, "jumps": [], "crossings": []}

    step_dists_sorted = sorted(step_dists)
    median_step = step_dists_sorted[len(step_dists_sorted) // 2]
    jump_threshold = median_step * jump_threshold_multiplier

    jumps = [(t1, t2, d) for (t1, t2), d in zip(usable_edges, step_dists) if d > jump_threshold]

    crossings = []
    for i in range(len(usable_edges)):
        a1, a2 = usable_edges[i]
        for j in range(i + 1, len(usable_edges)):
            b1, b2 = usable_edges[j]
            if len({a1, a2, b1, b2}) < 4:
                continue
            if segments_intersect(registered[a1], registered[a2], registered[b1], registered[b2]):
                crossings.append(((a1, a2), (b1, b2)))

    if jumps:
        print(f"RULE VIOLATION -- {len(jumps)} step(s) where sequential tiles ended up much farther apart than typical (>{jump_threshold:.1f}px, typical step ~{median_step:.1f}px):")
        for t1, t2, d in sorted(jumps, key=lambda x: -x[2])[:15]:
            print(f"    {t1} -> {t2}  (distance={d:.1f}px)")
    else:
        print(f"Step-distance check: OK -- every sequential step is within {jump_threshold:.1f}px (typical ~{median_step:.1f}px).")

    if crossings:
        print(f"RULE VIOLATION -- {len(crossings)} pair(s) of steps in the trace geometrically CROSS each other:")
        for (a1, a2), (b1, b2) in crossings[:15]:
            print(f"    [{a1} -> {a2}]  crosses  [{b1} -> {b2}]")
    else:
        print("Self-crossing check: OK -- the trace never crosses itself.")

    traceable = not jumps and not crossings
    if traceable:
        print(f">>> YES -- the filename-order trace is geometrically clean: consecutive tiles stay close together and the path never crosses itself. <<<")
    else:
        print(f">>> NO -- the filename-order trace has geometric problems (see above) -- something is wrong even though the scan order itself is fine. <<<")

    return {"traceable": traceable, "jumps": jumps, "crossings": crossings, "edges": usable_edges}


def analyze_graph_rules(pair_records, registered, tile_cluster=None, max_scan_jump=3, expect_loop=False):
    """
    Structural checks on the TRUSTED graph, independent of any single link's
    correlation/consistency: does the overall shape make sense for a simple
    chain/loop, or has something gone topologically wrong?

    expect_loop: if True, the acquisition-order jump check uses CIRCULAR
    distance (wrapping around from the highest-numbered tile back to the
    lowest) instead of raw difference. Without this, the one legitimate
    link that closes a ring (first tile <-> last tile) would always be
    flagged, since it has the largest possible index gap by definition.
    """
    trusted_edges = []
    for rec in pair_records:
        # 'trusted' here means: passed correlation AND (if we have cluster info)
        # ended up in the same cluster as its neighbor -- i.e. survived both
        # filters used elsewhere in this script.
        t1, t2 = rec["file1"], rec["file2"]
        if tile_cluster is not None:
            if tile_cluster.get(t1) is None or tile_cluster.get(t1) != tile_cluster.get(t2):
                continue
        if rec["valid"]:
            trusted_edges.append((t1, t2))

    trusted_edges = list(set(tuple(sorted(e)) for e in trusted_edges))  # de-dup

    degree = defaultdict(int)
    for t1, t2 in trusted_edges:
        degree[t1] += 1
        degree[t2] += 1

    print(f"\n=== Graph-theory structural checks on the trusted graph ({len(trusted_edges)} unique edges) ===")

    branch_points = {t: d for t, d in degree.items() if d > 2}
    endpoints = [t for t, d in degree.items() if d == 1]
    odd_degree_nodes = {t: d for t, d in degree.items() if d % 2 == 1}
    n_components_for_verdict = None
    if tile_cluster is not None:
        cluster_sizes_for_verdict = defaultdict(int)
        for t, c in tile_cluster.items():
            if t in registered:
                cluster_sizes_for_verdict[c] += 1
        n_components_for_verdict = len(cluster_sizes_for_verdict)

    # --- TOP-LINE VERDICT: can this be traced as one continuous stroke,
    # without lifting the pen, using ONLY positions and which links exist?
    # This is the classic EULERIAN TRAIL condition: a node can have any
    # number of connections and still be fine, AS LONG AS that number is
    # even (you pass through it, one edge in, one edge out, each time).
    # Only ODD-degree nodes actually break traceability -- and there must
    # be exactly 0 of them (closed loop) or exactly 2 (open path, which
    # become your start and end points).
    reasons = []
    if n_components_for_verdict is not None and n_components_for_verdict > 1:
        reasons.append(f"the structure is split into {n_components_for_verdict} disconnected pieces")
    if n_components_for_verdict is None or n_components_for_verdict <= 1:
        if len(odd_degree_nodes) not in (0, 2):
            reasons.append(f"{len(odd_degree_nodes)} tile(s) have an odd number of connections (expected exactly 0 for a closed loop or 2 for an open path): {list(odd_degree_nodes.keys())[:10]}{'...' if len(odd_degree_nodes) > 10 else ''}")

    # NOTE: this log-correlation-based verdict is auxiliary info only.
    # The DEFINITIVE traceability answer comes from check_index_guided_trace()
    # further down (filename-order trace, sanity-checked against positions) --
    # that one is what's actually reported as the final YES/NO. This log-based
    # version is kept only for its diagnostic detail (disconnected components,
    # degree info), since a tightly-folded strip can legitimately show many
    # real extra correlation-based links at each fold without that being a
    # genuine problem.
    if not reasons:
        shape_word = "CLOSED LOOP" if len(odd_degree_nodes) == 0 else "OPEN PATH"
        print(f"\n(Log-correlation-based check, informational only -- see the definitive index-guided trace verdict below): consistent with a {shape_word}.")
        if branch_points:
            print(f"    (Note: {len(branch_points)} tile(s) have more than 2 correlation-based links -- common at a tight fold in the path, not necessarily a problem.)")
        traceable_verdict = True
    else:
        print(f"\n(Log-correlation-based check, informational only -- see the definitive index-guided trace verdict below): {'; '.join(reasons)}.")
        traceable_verdict = False

    if n_components_for_verdict is not None:
        if n_components_for_verdict > 1:
            print(f"\nRULE VIOLATION -- {n_components_for_verdict} DISCONNECTED components found. "
                  f"The whole structure should be one single connected piece -- these clusters have "
                  f"no trusted link tying them together, so their relative position is unverified.")
        else:
            print("\nConnectivity rule: OK -- everything is one single connected piece, no isolated islands.")

    if branch_points:
        even_branch = {t: d for t, d in branch_points.items() if d % 2 == 0}
        odd_branch = {t: d for t, d in branch_points.items() if d % 2 == 1}
        if even_branch:
            print(f"\nInfo -- {len(even_branch)} tile(s) have more than 2 connections but EVEN degree (fine, just a revisited point):")
            for t, d in sorted(even_branch.items(), key=lambda kv: -kv[1]):
                print(f"    {t}: degree {d}")
        if odd_branch:
            print(f"\nRULE VIOLATION -- {len(odd_branch)} tile(s) have more than 2 connections AND odd degree (this breaks traceability):")
            for t, d in sorted(odd_branch.items(), key=lambda kv: -kv[1]):
                print(f"    {t}: degree {d}")
    else:
        print("\nDegree rule: OK -- no tile has more than 2 trusted links.")

    if len(odd_degree_nodes) == 0 and degree:
        print("Shape: looks like a CLOSED LOOP or fully-revisitable structure (0 odd-degree tiles).")
    elif len(odd_degree_nodes) == 2:
        print(f"Shape: looks like an OPEN PATH (2 odd-degree tiles = start/end: {list(odd_degree_nodes.keys())[0]}, {list(odd_degree_nodes.keys())[1]}).")
    elif degree:
        print(f"RULE VIOLATION -- {len(odd_degree_nodes)} odd-degree tile(s) found ({list(odd_degree_nodes.keys())[:10]}{'...' if len(odd_degree_nodes) > 10 else ''}); expected exactly 0 (loop) or 2 (path).")

    # Acquisition-order jump check
    all_indices = [extract_scan_index(t) for t in set(t for e in trusted_edges for t in e)]
    all_indices = [i for i in all_indices if i is not None]
    idx_range = (max(all_indices) - min(all_indices)) if all_indices else None

    jump_violations = []
    for t1, t2 in trusted_edges:
        i1, i2 = extract_scan_index(t1), extract_scan_index(t2)
        if i1 is None or i2 is None:
            continue
        raw_diff = abs(i1 - i2)
        if expect_loop and idx_range:
            jump = min(raw_diff, idx_range - raw_diff)
        else:
            jump = raw_diff
        if jump > max_scan_jump:
            jump_violations.append((t1, t2, jump))
    if jump_violations:
        mode_note = "circular distance (--loop mode)" if expect_loop else "raw numeric difference"

        # Corroboration check: an isolated long-jump link (no other nearby tile
        # pair showing a similarly-long jump between the same two neighborhoods)
        # is much more likely a one-off spurious match than genuine physical
        # proximity (e.g. a 'roller coaster' curve that swoops close without
        # crossing). A genuine close-approach should show up as a small
        # NEIGHBORHOOD of consistent links, not a single lone connection.
        all_jump_edges = [(t1, t2, extract_scan_index(t1), extract_scan_index(t2)) for t1, t2, _ in jump_violations]

        isolated = []
        corroborated = []
        for t1, t2, jump in jump_violations:
            i1, i2 = extract_scan_index(t1), extract_scan_index(t2)
            has_support = False
            for ot1, ot2, oi1, oi2 in all_jump_edges:
                if (ot1, ot2) == (t1, t2):
                    continue
                if abs(oi1 - i1) <= max_scan_jump and abs(oi2 - i2) <= max_scan_jump:
                    has_support = True
                    break
            if has_support:
                corroborated.append((t1, t2, jump))
            else:
                isolated.append((t1, t2, jump))

        if isolated:
            print(f"\nRULE VIOLATION -- {len(isolated)} ISOLATED trusted link(s) jump more than {max_scan_jump} in acquisition order ({mode_note}), with NO other nearby tile pair supporting the same connection -- likely spurious matches, not real proximity:")
            for t1, t2, jump in sorted(isolated, key=lambda x: -x[2])[:20]:
                print(f"    {t1} <-> {t2}  (index jump = {jump})")
        if corroborated:
            print(f"\n{len(corroborated)} long-jump link(s) ARE corroborated by nearby tile pairs making similar connections -- consistent with a genuine physical close-approach (e.g. a curve that swoops near itself), not flagged as violations, but still worth a glance:")
            for t1, t2, jump in sorted(corroborated, key=lambda x: -x[2])[:20]:
                print(f"    {t1} <-> {t2}  (index jump = {jump})")
        if not isolated and not corroborated:
            print(f"\nAcquisition-order rule: OK -- no trusted link jumps more than {max_scan_jump} in scan order.")
    else:
        isolated, corroborated = [], []
        print(f"\nAcquisition-order rule: OK -- no trusted link jumps more than {max_scan_jump} in scan order.")

    # Crossing check (only meaningful if positions are available)
    crossings = []
    edge_list = [(t1, t2) for t1, t2 in trusted_edges if t1 in registered and t2 in registered]
    for i in range(len(edge_list)):
        a1, a2 = edge_list[i]
        for j in range(i + 1, len(edge_list)):
            b1, b2 = edge_list[j]
            if len({a1, a2, b1, b2}) < 4:
                continue  # shares an endpoint, not a real crossing
            if segments_intersect(registered[a1], registered[a2], registered[b1], registered[b2]):
                crossings.append(((a1, a2), (b1, b2)))
    if crossings:
        print(f"\nRULE VIOLATION -- {len(crossings)} pair(s) of trusted edges geometrically CROSS each other (shouldn't happen in a simple non-self-intersecting strip):")
        for (a1, a2), (b1, b2) in crossings[:20]:
            print(f"    [{a1} <-> {a2}]  crosses  [{b1} <-> {b2}]")
    else:
        print("\nCrossing rule: OK -- no trusted edges cross each other.")

    return {
        "branch_points": branch_points,
        "endpoints": endpoints,
        "jump_violations": isolated,
        "corroborated_proximity": corroborated,
        "crossings": crossings,
        "n_components": n_components_for_verdict,
        "traceable": traceable_verdict,
        "traceable_reasons": reasons,
    }


def plot_topology_graph(pair_records, registered, out_html, shift_discrepancy_threshold=50.0, tile_cluster=None, rule_violations=None):
    if not registered:
        print("No registered positions available to draw a topology graph.")
        return

    edge_traces = {"trusted": {"x": [], "y": []}, "rejected": {"x": [], "y": []}, "suspicious": {"x": [], "y": []}}

    for rec in pair_records:
        t1, t2 = rec["file1"], rec["file2"]
        pos1, pos2 = registered.get(t1), registered.get(t2)
        if not pos1 or not pos2:
            continue

        geometrically_consistent = True
        dx, dy = rec.get("dx"), rec.get("dy")
        if dx is not None and dy is not None:
            actual_dx = pos2[0] - pos1[0]
            actual_dy = pos2[1] - pos1[1]
            vector_error = ((actual_dx - dx) ** 2 + (actual_dy - dy) ** 2) ** 0.5
            if vector_error > shift_discrepancy_threshold:
                geometrically_consistent = False

        if rec["valid"] and geometrically_consistent:
            key = "trusted"
        elif rec["valid"] and not geometrically_consistent:
            key = "suspicious"
        else:
            key = "rejected"

        edge_traces[key]["x"] += [pos1[0], pos2[0], None]
        edge_traces[key]["y"] += [pos1[1], pos2[1], None]

    fig = go.Figure()

    style = {
        "rejected": dict(color="lightgray", width=1, dash="dot", label="Rejected link"),
        "suspicious": dict(color="orange", width=2, dash="dash", label="High-corr but geometrically contradicted"),
        "trusted": dict(color="seagreen", width=2, dash="solid", label="Trusted & consistent link"),
    }
    for key in ["rejected", "suspicious", "trusted"]:
        et = edge_traces[key]
        if et["x"]:
            fig.add_trace(
                go.Scatter(
                    x=et["x"], y=et["y"],
                    mode="lines",
                    line=dict(color=style[key]["color"], width=style[key]["width"], dash=style[key]["dash"]),
                    name=style[key]["label"],
                    hoverinfo="skip",
                    visible="legendonly",
                )
            )

    tiles = sorted(registered)
    xs = [registered[t][0] for t in tiles]
    ys = [registered[t][1] for t in tiles]

    marker_kwargs = dict(size=12, line=dict(width=1, color="black"))
    if tile_cluster:
        marker_kwargs["color"] = [tile_cluster.get(t, -1) for t in tiles]
        marker_kwargs["colorscale"] = "Turbo"
        marker_kwargs["showscale"] = False
    else:
        marker_kwargs["color"] = "steelblue"

    # Highlight tiles actually involved in the DEFINITIVE verdict (the
    # index-guided trace: sequential_jump / trace_crossing), not the old
    # log-correlation branch-point count. A tightly-folded strip will
    # legitimately show many "extra" trusted links at each fold (real,
    # high-correlation overlaps between tiles close in space but far apart
    # in scan order) -- that's not a mistake, so it shouldn't be flagged.
    trace_problem_tiles = set()
    for t1, t2, _ in (rule_violations or {}).get("trace_jumps", []):
        trace_problem_tiles.add(t1)
        trace_problem_tiles.add(t2)
    for (a1, a2), (b1, b2) in (rule_violations or {}).get("trace_crossings", []):
        trace_problem_tiles.update([a1, a2, b1, b2])

    if trace_problem_tiles:
        marker_kwargs["size"] = [22 if t in trace_problem_tiles else 12 for t in tiles]
        marker_kwargs["line"] = dict(width=[3 if t in trace_problem_tiles else 1 for t in tiles], color=["red" if t in trace_problem_tiles else "black" for t in tiles])

    fig.add_trace(
        go.Scatter(
            x=xs, y=ys,
            mode="markers+text",
            text=tiles,
            textposition="top center",
            textfont=dict(size=7),
            marker=marker_kwargs,
            name="Tile",
            hovertemplate="%{text}<br>(%{x:.1f}, %{y:.1f})<extra></extra>",
            visible="legendonly",
        )
    )

    jump_violations = (rule_violations or {}).get("jump_violations", [])
    if jump_violations:
        jx, jy = [], []
        for t1, t2, _ in jump_violations:
            if t1 in registered and t2 in registered:
                jx += [registered[t1][0], registered[t2][0], None]
                jy += [registered[t1][1], registered[t2][1], None]
        if jx:
            fig.add_trace(go.Scatter(x=jx, y=jy, mode="lines", line=dict(color="red", width=3, dash="dashdot"), name="Isolated jump (likely spurious)", hoverinfo="skip", visible="legendonly"))

    corroborated_proximity = (rule_violations or {}).get("corroborated_proximity", [])
    if corroborated_proximity:
        px, py = [], []
        for t1, t2, _ in corroborated_proximity:
            if t1 in registered and t2 in registered:
                px += [registered[t1][0], registered[t2][0], None]
                py += [registered[t1][1], registered[t2][1], None]
        if px:
            fig.add_trace(go.Scatter(x=px, y=py, mode="lines", line=dict(color="gold", width=2, dash="dashdot"), name="Corroborated jump (likely genuine proximity)", hoverinfo="skip", visible="legendonly"))

    spatial_isolated = (rule_violations or {}).get("spatial_isolated", [])
    if spatial_isolated:
        sx, sy = [], []
        for t1, t2, _, _ in spatial_isolated:
            if t1 in registered and t2 in registered:
                sx += [registered[t1][0], registered[t2][0], None]
                sy += [registered[t1][1], registered[t2][1], None]
        if sx:
            fig.add_trace(go.Scatter(x=sx, y=sy, mode="lines", line=dict(color="purple", width=3, dash="solid"), name="Unexpected spatial neighbor (pure position check)", hoverinfo="skip", visible="legendonly"))

    spatial_corroborated = (rule_violations or {}).get("spatial_corroborated", [])
    if spatial_corroborated:
        cpx, cpy = [], []
        for t1, t2, _, _ in spatial_corroborated:
            if t1 in registered and t2 in registered:
                cpx += [registered[t1][0], registered[t2][0], None]
                cpy += [registered[t1][1], registered[t2][1], None]
        if cpx:
            fig.add_trace(go.Scatter(x=cpx, y=cpy, mode="lines", line=dict(color="lightpink", width=2, dash="dot"), name="Spatial proximity (corroborated, likely genuine)", hoverinfo="skip", visible="legendonly"))

    crossings = (rule_violations or {}).get("crossings", [])
    if crossings:
        cx, cy = [], []
        for (a1, a2), (b1, b2) in crossings:
            for (p, q) in [(a1, a2), (b1, b2)]:
                if p in registered and q in registered:
                    cx += [registered[p][0], registered[q][0], None]
                    cy += [registered[p][1], registered[q][1], None]
        if cx:
            fig.add_trace(go.Scatter(x=cx, y=cy, mode="lines", line=dict(color="magenta", width=4, dash="solid"), name="Crossing edge (rule violation)", hoverinfo="skip", visible="legendonly"))

    fig.update_layout(
        title="Tile connectivity graph at actual registered positions<br><sup>Green = trusted, orange = high-corr but contradicted, gray dotted = rejected. Red ring = tile involved in a trace jump/crossing (the actual verdict). Purple = unexpected spatial neighbor.</sup>",
        xaxis_title="X position (registered)",
        yaxis_title="Y position (registered)",
        yaxis=dict(autorange="reversed", scaleanchor="x", scaleratio=1),
        width=1100,
        height=950,
    )

    fig.write_html(out_html)
    print(f"\nSaved topology graph to {out_html}")


def write_issues_csv(rule_violations, out_csv):
    """Flatten every flagged issue for one sequence into rows of a CSV: category, tile1, tile2, detail."""
    rows = []

    for t1, t2, d in rule_violations.get("trace_jumps", []):
        rows.append({"category": "sequential_jump", "tile1": t1, "tile2": t2, "tile3": "", "tile4": "", "detail": f"distance={d:.1f}px"})

    for (a1, a2), (b1, b2) in rule_violations.get("trace_crossings", []):
        rows.append({"category": "trace_crossing", "tile1": a1, "tile2": a2, "tile3": b1, "tile4": b2, "detail": "these two steps cross"})

    for t, d in rule_violations.get("branch_points", {}).items():
        rows.append({"category": "branch_point", "tile1": t, "tile2": "", "tile3": "", "tile4": "", "detail": f"degree={d}"})

    for t1, t2, idx_diff, d in rule_violations.get("spatial_isolated", []):
        rows.append({"category": "unexpected_spatial_neighbor", "tile1": t1, "tile2": t2, "tile3": "", "tile4": "", "detail": f"scan_gap={idx_diff}, distance={d:.1f}px"})

    for t1, t2, jump in rule_violations.get("jump_violations", []):
        rows.append({"category": "log_link_acquisition_jump", "tile1": t1, "tile2": t2, "tile3": "", "tile4": "", "detail": f"index_jump={jump}"})

    for (a1, a2), (b1, b2) in rule_violations.get("crossings", []):
        rows.append({"category": "trusted_link_crossing", "tile1": a1, "tile2": a2, "tile3": b1, "tile4": b2, "detail": "these two trusted links cross"})

    n_components = rule_violations.get("n_components")
    if n_components and n_components > 1:
        rows.append({"category": "disconnected_components", "tile1": "", "tile2": "", "tile3": "", "tile4": "", "detail": f"{n_components} separate pieces"})

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["category", "tile1", "tile2", "tile3", "tile4", "detail"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {len(rows)} flagged issue(s) to {out_csv}")
    return len(rows)


def write_summary_csv(summary_rows, out_csv):
    """One row per sequence: name, traceable YES/NO, reasons, and issue count."""
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sequence", "traceable", "reasons", "n_issues"])
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"\nSaved summary CSV to {out_csv}")


def process_one_sequence(log_files, registered_path, out_html, shift_discrepancy_threshold=50.0, use_cluster_colors=True, max_scan_jump=3, expect_loop=False, spatial_threshold_multiplier=1.5):
    print(f"Log file(s): {[os.path.basename(l) for l in log_files]}")
    print(f"Registered config: {registered_path}")

    all_records = []
    for lf in log_files:
        recs = parse_log_file(lf)
        print(f"  {os.path.basename(lf)}: {len(recs)} pairwise comparisons")
        all_records.extend(recs)

    registered = parse_tile_configuration(registered_path)

    tile_cluster = analyze_clusters(all_records, registered, shift_discrepancy_threshold) if use_cluster_colors else None
    rule_violations = analyze_graph_rules(all_records, registered, tile_cluster, max_scan_jump=max_scan_jump, expect_loop=expect_loop)
    spatial_violations = analyze_spatial_proximity_violations(registered, max_scan_jump=max_scan_jump, spatial_threshold_multiplier=spatial_threshold_multiplier, expect_loop=expect_loop)
    rule_violations["spatial_isolated"] = spatial_violations["isolated"]
    rule_violations["spatial_corroborated"] = spatial_violations["corroborated"]

    # Definitive traceability verdict: trace built from FILENAMES (scan
    # order), sanity-checked against positions (no unexpected jumps, no
    # self-crossing). This is what's used for the top-line YES/NO and the
    # batch summary.
    trace_verdict = check_index_guided_trace(registered, expect_loop=expect_loop)
    rule_violations["traceable"] = trace_verdict["traceable"]
    reasons = []
    if trace_verdict.get("jumps"):
        reasons.append(f"{len(trace_verdict['jumps'])} unexpected jump(s) between sequential tiles")
    if trace_verdict.get("crossings"):
        reasons.append(f"{len(trace_verdict['crossings'])} self-crossing(s) in the trace")
    rule_violations["traceable_reasons"] = reasons
    rule_violations["trace_jumps"] = trace_verdict.get("jumps", [])
    rule_violations["trace_crossings"] = trace_verdict.get("crossings", [])

    plot_topology_graph(all_records, registered, out_html, shift_discrepancy_threshold, tile_cluster, rule_violations)

    out_csv = os.path.splitext(out_html)[0] + "_issues.csv"
    n_issues = write_issues_csv(rule_violations, out_csv)
    rule_violations["n_issues"] = n_issues

    return rule_violations


def run_batch(project_root, out_dir=None, shift_discrepancy_threshold=50.0, use_cluster_colors=True, max_scan_jump=3, expect_loop=False, spatial_threshold_multiplier=1.5):
    slices_dir = os.path.join(project_root, "aorta")
    log_dir = os.path.join(project_root, "log_dir")

    if not os.path.isdir(slices_dir):
        print(f"Expected a 'slices' folder inside {project_root}, none found.", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(log_dir):
        print(f"Expected a 'log' folder inside {project_root}, none found.", file=sys.stderr)
        sys.exit(1)

    seq_dirs = sorted(d for d in glob.glob(os.path.join(slices_dir, "*")) if os.path.isdir(d))
    if not seq_dirs:
        print(f"No sequence subfolders found inside {slices_dir}", file=sys.stderr)
        sys.exit(1)

    all_log_files = sorted(glob.glob(os.path.join(log_dir, "*.txt")) + glob.glob(os.path.join(log_dir, "*.log")))

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    summary = []

    for seq_dir in seq_dirs:
        seq_name = os.path.basename(seq_dir.rstrip("/"))
        print(f"\n=== Sequence: {seq_name} ===")

        matched_logs = [lf for lf in all_log_files if seq_name.lower() in os.path.basename(lf).lower()]
        if not matched_logs:
            print(f"  No log files matched '{seq_name}' in {log_dir} -- skipping.")
            summary.append({"sequence": seq_name, "traceable": "SKIPPED", "reasons": "no matching logs", "n_issues": 0})
            continue

        registered_path = find_registered_file(seq_dir)
        if not registered_path:
            print(f"  No TileConfiguration.registered.txt found under {seq_dir} -- skipping.")
            summary.append({"sequence": seq_name, "traceable": "SKIPPED", "reasons": "no registered file", "n_issues": 0})
            continue

        out_html = os.path.join(out_dir, f"{seq_name}_topology.html") if out_dir else os.path.join(seq_dir, "topology.html")
        rule_violations = process_one_sequence(matched_logs, registered_path, out_html, shift_discrepancy_threshold, use_cluster_colors, max_scan_jump=max_scan_jump, expect_loop=expect_loop, spatial_threshold_multiplier=spatial_threshold_multiplier)
        traceable = rule_violations.get("traceable")
        reasons = "; ".join(rule_violations.get("traceable_reasons", [])) if not traceable else ""
        summary.append({"sequence": seq_name, "traceable": "YES" if traceable else "NO", "reasons": reasons, "n_issues": rule_violations.get("n_issues", 0)})

    print("\n\n=== SUMMARY: traceable without lifting the pen? ===")
    name_width = max((len(s["sequence"]) for s in summary), default=10) + 2
    for s in summary:
        verdict = f"{s['traceable']}" + (f" ({s['reasons']})" if s["reasons"] else "")
        print(f"  {s['sequence']:<{name_width}} {verdict}")

    summary_csv = os.path.join(out_dir, "summary.csv") if out_dir else os.path.join(project_root, "summary.csv")
    write_summary_csv(summary, summary_csv)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--tile-path", help="Full path to ONE tile image file inside the sequence folder, e.g. '/path/to/project/slices/MySeq/Snap-1.czi'. Auto-derives the sequence folder and looks for logs matching that sequence name inside a sibling 'log/' folder. Processes just that one sequence.")
    mode.add_argument("--log-dir", help="Manual mode: directory containing log file(s)")
    mode.add_argument("--project-root", help="Batch mode: project folder containing 'log/' and 'slices/<seq>/' subfolders -- processes EVERY sequence automatically in one run.")

    parser.add_argument("--registered", help="Manual mode: path to TileConfiguration.registered.txt (required with --log-dir)")
    parser.add_argument("--out", default="topology.html", help="Output HTML path (single-sequence modes only)")
    parser.add_argument("--out-dir", default=None, help="Batch mode only: write all HTMLs into this one folder (named '<sequence>_topology.html') instead of inside each sequence's own folder")
    parser.add_argument("--shift-discrepancy-threshold", type=float, default=50.0, help="Pixel threshold for treating a high-correlation link as contradicted by the final registration (default 50px)")
    parser.add_argument("--no-cluster-colors", action="store_true", help="Don't color nodes by detected cluster")
    parser.add_argument("--max-scan-jump", type=int, default=3, help="A trusted link connecting tiles whose trailing filename number differs by more than this is flagged as an acquisition-order violation (default 3)")
    parser.add_argument("--loop", action="store_true", help="Tell the acquisition-order check that this structure is a CLOSED LOOP (e.g. a vessel cross-section) -- uses circular distance so the legitimate seam link (lowest-numbered tile <-> highest-numbered tile) isn't wrongly flagged.")
    parser.add_argument("--spatial-threshold-multiplier", type=float, default=1.5, help="For the pure-position proximity check: flag tiles closer than (typical tile spacing * this multiplier) if they are not expected to be neighbors. Default 1.5.")
    args = parser.parse_args()

    if args.project_root:
        run_batch(args.project_root, out_dir=args.out_dir, shift_discrepancy_threshold=args.shift_discrepancy_threshold, use_cluster_colors=not args.no_cluster_colors, max_scan_jump=args.max_scan_jump, expect_loop=args.loop, spatial_threshold_multiplier=args.spatial_threshold_multiplier)
        return

    if args.tile_path:
        tile_path = os.path.normpath(args.tile_path)
        sequence_dir = os.path.dirname(tile_path)
        slices_dir = os.path.dirname(sequence_dir)
        project_root = os.path.dirname(slices_dir)
        seq_name = os.path.basename(sequence_dir.rstrip("/"))
        log_dir = os.path.join(project_root, "log")

        print(f"Derived:\n  project root: {project_root}\n  sequence dir: {sequence_dir}\n  log dir:      {log_dir}\n")

        log_files = [
            lf for lf in sorted(glob.glob(os.path.join(log_dir, "*.txt")) + glob.glob(os.path.join(log_dir, "*.log")))
            if seq_name.lower() in os.path.basename(lf).lower()
        ]
        registered_path = find_registered_file(sequence_dir)
    else:
        if not args.registered:
            parser.error("--registered is required when using --log-dir")
        log_files = sorted(glob.glob(os.path.join(args.log_dir, "*_log.txt")))
        registered_path = args.registered

    if not log_files:
        print("No matching log files found.", file=sys.stderr)
        sys.exit(1)
    if not registered_path:
        print("No TileConfiguration.registered.txt found.", file=sys.stderr)
        sys.exit(1)

    process_one_sequence(log_files, registered_path, args.out, args.shift_discrepancy_threshold, not args.no_cluster_colors, max_scan_jump=args.max_scan_jump, expect_loop=args.loop, spatial_threshold_multiplier=args.spatial_threshold_multiplier)


if __name__ == "__main__":
    main()
