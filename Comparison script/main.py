import re
import math

def parse_fortran_nodes(filepath):
    nodes = {}
    edge = {}
    triangle = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.strip().startswith("Node"):
                match = re.match(
                    r"Node\s+(\d+)\s+-> x:\s+([-0-9.E+]+)\s+y:\s+([-0-9.E+]+)\s+z:\s+([-0-9.E+]+)\s+n_neighbor:\s+(\d+)\s+boundary:\s+([TF])\s+state:\s+(\d+)\s+voronoi_area:\s+([-0-9.E+]+)",
                    line.strip()
                )
                if match:
                    idx = int(match.group(1))
                    nodes[idx] = {
                        "x": float(match.group(2)),
                        "y": float(match.group(3)),
                        "z": float(match.group(4)),
                        "n_neighbor": int(match.group(5)),
                        "boundary": match.group(6) == "T",
                        "state": int(match.group(7)),
                        "voronoi_area": float(match.group(8)),
                    }
            elif line.strip().startswith("Triangle"):
                match = re.match(
                    # r"Triangle\s+(\d+)\s+-> x:\s+([-0-9.E+]+)\s+y:\s+([-0-9.E+]+)\s+z:\s+([-0-9.E+]+)\s+n_neighbor:\s+(\d+)\s+boundary:\s+([TF])\s+state:\s+(\d+)\s+voronoi_area:\s+([-0-9.E+]+)",
                    r"Triangle\s+(\d+)\s+-> center_x:\s+([-0-9.E+]+)\s+center_y:\s+([-0-9.E+]+)\s+radius:\s+([-0-9.E+]+)\s+area:\s+([-0-9.E+]+)\s+aspect_ratio:\s+([-0-9.E+]+)\s+",
                    line.strip()
                )
                # Triangle           1 -> center_x:  -49.5773201      center_y:  -29.4167919      radius:   3.48062813E-02  area:   8.68029019E-06  aspect_ratio:   3.50617655E-02
                if match:
                    idx = int(match.group(1))
                    triangle[idx] = {
                        "center_x": float(match.group(2)),
                        "center_y": float(match.group(3)),
                        "radius": float(match.group(4)),
                        "area": float(match.group(5)),
                        "aspect_ratio": float(match.group(6))
                    }
            elif line.strip().startswith("Edge"):
                match = re.match(
                    # r"Edge\s+(\d+)\s+-> x:\s+([-0-9.E+]+)\s+y:\s+([-0-9.E+]+)\s+z:\s+([-0-9.E+]+)\s+n_neighbor:\s+(\d+)\s+boundary:\s+([TF])\s+state:\s+(\d+)\s+voronoi_area:\s+([-0-9.E+]+)",
                    r"Edge\s+(\d+)\s+-> id:\s+(\d+)\s+start:\s+(\d+)\s+end:\s+(\d+)\s+voronoi_x:\s+([-0-9.E+]+)\s+voronoi_y:\s+([-0-9.E+]+)\s+length:\s+([-0-9.E+]+)\s+slope:\s+([-0-9.E+]+)\s+boundary:\s+([TF])\s",
                    line.strip()
                )
                #  Edge        2931 -> id:        1415  start:         500  end:         466  voronoi_x:  -49.4125214      voronoi_y:  -29.3069267      length:   5.89246675E-03  slope:   0.00000000      boundary: F
                if match:
                    idx = int(match.group(1))
                    edge[idx] = {
                        "id": int(match.group(2)),
                        "start": int(match.group(3)),
                        "end": int(match.group(4)),
                        "voronoi_x": float(match.group(5)),
                        "voronoi_y": float(match.group(6)),
                        "length": float(match.group(7)),
                        "slope": float(match.group(8)),
                        "boundary": match.group(9) == "T",
                    }
    return nodes, edge, triangle

def parse_rust_nodes(filepath):
    nodes = {}
    edge = {}
    triangle = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.strip().startswith("Node"):
                match = re.match(
                    r"Node\s+(\d+)\s+-> x:\s+([-0-9.E+]+), y:\s+([-0-9.E+]+), z:\s+([-0-9.E+]+), n_neighbor:\s+(\d+), boundary:\s+(true|false), state:\s+(\d+), voronoi_area:\s+([-0-9.E+]+)",
                    line.strip()
                )
                if match:
                    idx = int(match.group(1))
                    nodes[idx] = {
                        "x": float(match.group(2)),
                        "y": float(match.group(3)),
                        "z": float(match.group(4)),
                        "n_neighbor": int(match.group(5)),
                        "boundary": match.group(6) == "true",
                        "state": int(match.group(7)),
                        "voronoi_area": float(match.group(8)),
                    }
            elif line.strip().startswith("Triangle"):
                match = re.match(
                    r"Triangle\s+(\d+)\s+-> center_x:\s+([-0-9.E+]+), center_y:\s+([-0-9.E+]+), radius:\s+([-0-9.E+]+), area:\s+([-0-9.E+]+), aspect_ratio:\s+([-0-9.E+]+)",
                    line.strip()
                )
                if match:
                    idx = int(match.group(1))
                    triangle[idx] = {
                        "center_x": float(match.group(2)),
                        "center_y": float(match.group(3)),
                        "radius": float(match.group(4)),
                        "area": float(match.group(5)),
                        "aspect_ratio": float(match.group(6)),
                    }
            elif line.strip().startswith("Edge"):
                match = re.match(
                    r"Edge\s+(\d+)\s+-> id:\s+(\d+), start:\s+(\d+), end:\s+(\d+), voronoi_x:\s+([-0-9.E+]+), voronoi_y:\s+([-0-9.E+]+), length:\s+([-0-9.E+]+), slope:\s+(\d+), boundary:\s+(true|false)",
                    line.strip()
                )
                if match:
                    idx = int(match.group(1))
                    edge[idx] = {
                        "id": int(match.group(2)),
                        "start": int(match.group(3)),
                        "end": int(match.group(4)),
                        "voronoi_x": float(match.group(5)),
                        "voronoi_y": float(match.group(6)),
                        "length": float(match.group(7)),
                        "slope": float(match.group(8)),
                        "boundary": match.group(9) == "true"
                    }
    return nodes, edge, triangle

def compare_nodes(fortran_nodes, rust_nodes, tol=1e-6):
    discrepancies = []
    for f_idx, f_data in fortran_nodes.items():
        r_idx = f_idx - 1  # ajustar índice
        r_data = rust_nodes.get(r_idx)
        if not r_data:
            discrepancies.append((f_idx, "missing in rust"))
            continue
        for key in ["x", "y", "z"]:
            if not math.isclose(f_data[key], r_data[key], rel_tol=tol, abs_tol=tol):
                discrepancies.append((f_idx, f"{key} mismatch: {f_data[key]} vs {r_data[key]}"))
        if f_data["n_neighbor"] != r_data["n_neighbor"]:
            discrepancies.append((f_idx, f"n_neighbor mismatch: {f_data['n_neighbor']} vs {r_data['n_neighbor']}"))
        if f_data["boundary"] != r_data["boundary"]:
            discrepancies.append((f_idx, f"boundary mismatch: {f_data['boundary']} vs {r_data['boundary']}"))
        # state e voronoi_area talvez não sejam relevantes, mas podemos verificar:
        if f_data["state"] != r_data["state"]:
            discrepancies.append((f_idx, f"state mismatch: {f_data['state']} vs {r_data['state']}"))
        if not math.isclose(f_data["voronoi_area"], r_data["voronoi_area"], rel_tol=tol, abs_tol=tol):
            discrepancies.append((f_idx, f"voronoi_area mismatch: {f_data['voronoi_area']} vs {r_data['voronoi_area']}"))
    return discrepancies

def compare_edges(fortran_edges, rust_edges, tol=1e-6):
    discrepancies = []
    #  edge[idx] = {
    #                     "id": int(match.group(2)),
    #                     "start": int(match.group(3)),
    #                     "end": int(match.group(4)),
    #                     "voronoi_x": float(match.group(5)),
    #                     "voronoi_y": float(match.group(6)),
    #                     "length": float(match.group(7)),
    #                     "slope": int(match.group(8)),
    #                     "boundary": match.group(9) == "true"
    #                 }
    for f_idx, f_data in fortran_edges.items():
        r_idx = f_idx - 1  # ajustar índice
        r_data = rust_edges.get(r_idx)
        if not r_data:
            discrepancies.append((f_idx, "missing in rust"))
            continue
        for key in ["x", "y", "z"]:
            if not math.isclose(f_data[key], r_data[key], rel_tol=tol, abs_tol=tol):
                discrepancies.append((f_idx, f"{key} mismatch: {f_data[key]} vs {r_data[key]}"))
        if f_data["n_neighbor"] != r_data["n_neighbor"]:
            discrepancies.append((f_idx, f"n_neighbor mismatch: {f_data['n_neighbor']} vs {r_data['n_neighbor']}"))
        if f_data["boundary"] != r_data["boundary"]:
            discrepancies.append((f_idx, f"boundary mismatch: {f_data['boundary']} vs {r_data['boundary']}"))
        # state e voronoi_area talvez não sejam relevantes, mas podemos verificar:
        if f_data["state"] != r_data["state"]:
            discrepancies.append((f_idx, f"state mismatch: {f_data['state']} vs {r_data['state']}"))
        if not math.isclose(f_data["voronoi_area"], r_data["voronoi_area"], rel_tol=tol, abs_tol=tol):
            discrepancies.append((f_idx, f"voronoi_area mismatch: {f_data['voronoi_area']} vs {r_data['voronoi_area']}"))
    return discrepancies

def compare_triangles(fortran_nodes, rust_nodes, tol=1e-6):
    discrepancies = []
    for f_idx, f_data in fortran_nodes.items():
        r_idx = f_idx - 1  # ajustar índice
        r_data = rust_nodes.get(r_idx)
        if not r_data:
            discrepancies.append((f_idx, "missing in rust"))
            continue
        for key in ["x", "y", "z"]:
            if not math.isclose(f_data[key], r_data[key], rel_tol=tol, abs_tol=tol):
                discrepancies.append((f_idx, f"{key} mismatch: {f_data[key]} vs {r_data[key]}"))
        if f_data["n_neighbor"] != r_data["n_neighbor"]:
            discrepancies.append((f_idx, f"n_neighbor mismatch: {f_data['n_neighbor']} vs {r_data['n_neighbor']}"))
        if f_data["boundary"] != r_data["boundary"]:
            discrepancies.append((f_idx, f"boundary mismatch: {f_data['boundary']} vs {r_data['boundary']}"))
        # state e voronoi_area talvez não sejam relevantes, mas podemos verificar:
        if f_data["state"] != r_data["state"]:
            discrepancies.append((f_idx, f"state mismatch: {f_data['state']} vs {r_data['state']}"))
        if not math.isclose(f_data["voronoi_area"], r_data["voronoi_area"], rel_tol=tol, abs_tol=tol):
            discrepancies.append((f_idx, f"voronoi_area mismatch: {f_data['voronoi_area']} vs {r_data['voronoi_area']}"))
    return discrepancies


if __name__ == "__main__":
    fortran_nodes, fortran_edges, fortran_triangles = parse_fortran_nodes("../Data/Output/Fortran/inicial_1_000.txt")
    rust_nodes, rust_edges, rust_triangles = parse_rust_nodes("../Data/Output/Rust_Sequencial/inicial_1_000.txt")
    node_diffs = compare_nodes(fortran_nodes, rust_nodes)
    edge_diffs = compare_nodes(fortran_edges, rust_edges)
    triangle_diffs = compare_nodes(fortran_triangles, rust_triangles)

    if not node_diffs:
        print("Arquivos são consistentes dentro da tolerância definida.")
    else:
        print("Discrepâncias encontradas:")
        for d in node_diffs:
            print(d)


# ---------------------------------------------------------------------------------


#!/usr/bin/env python3
"""
Comparador de saídas Fortran x Rust (Nodes, Edges, Triangles).

Uso:
    python compare_fortran_rust.py fortran_inicial_1_000.txt rust_inicial_1_000.txt

Saída:
    Resumo de quantos elementos foram comparados e lista (limitada) de mismatches.
"""
# import re
# import math
# import sys
# import argparse
# from collections import defaultdict

# # ---------------------------
# # Helpers numéricos / parser
# # ---------------------------
# def ffloat(s):
#     try:
#         return float(s)
#     except Exception:
#         s2 = s.replace('D', 'E').replace('d', 'e').replace(',', '.')
#         return float(s2)

# def is_close(a, b, rtol=1e-6, atol=1e-12):
#     if a is None or b is None:
#         return a == b
#     return abs(a - b) <= atol + rtol * abs(b)

# # ---------------------------
# # Regex patterns (fortran & rust tolerant)
# # ---------------------------
# # Node lines: handle "x: -49.4, y: -29.4, ..." or "x:  -49.4354172      y:  -29.3645840  ..."
# NODE_PATTERNS = [
#     re.compile(r'Node\s+(\d+)\s*->\s*x:\s*([-\dEe+.]+)\s*,?\s*y:\s*([-\dEe+.]+)\s*,?\s*z:\s*([-\dEe+.]+)\s*,?\s*n_neighbor:\s*(\d+)\s*,?\s*boundary:\s*(T|F|true|false)\s*,?\s*state:\s*([-\d]+)\s*,?\s*voronoi_area:\s*([-\dEe+.]+)', re.IGNORECASE),
#     re.compile(r'Node\s+(\d+)\s*->\s*x:\s*([-\dEe+.]+)\s*y:\s*([-\dEe+.]+)\s*z:\s*([-\dEe+.]+)\s*n_neighbor:\s*(\d+)\s+boundary:\s*(T|F|true|false)\s+state:\s*([-\d]+)\s+voronoi_area:\s*([-\dEe+.]+)', re.IGNORECASE),
# ]

# # Edge lines:
# EDGE_PATTERNS = [
#     re.compile(r'Edge\s+(\d+)\s*->\s*id:\s*(\d+)\s*,?\s*start:\s*(\d+)\s*,?\s*end:\s*(\d+)\s*,?\s*voronoi_x:\s*([-\dEe+.]+)\s*,?\s*voronoi_y:\s*([-\dEe+.]+)\s*,?\s*length:\s*([-\dEe+.]+)\s*,?\s*slope:\s*([-\dEe+.]+)\s*,?\s*boundary:\s*(T|F|true|false)', re.IGNORECASE),
#     re.compile(r'Edge\s+(\d+)\s*->\s*id:\s*(\d+)\s*start:\s*(\d+)\s*end:\s*(\d+)\s*voronoi_x:\s*([-\dEe+.]+)\s*voronoi_y:\s*([-\dEe+.]+)\s*length:\s*([-\dEe+.]+)\s*slope:\s*([-\dEe+.]+)\s*boundary:\s*(T|F|true|false)', re.IGNORECASE),
# ]

# # Triangle lines:
# TRI_PATTERNS = [
#     re.compile(r'Triangle\s+(\d+)\s*->\s*center_x:\s*([-\dEe+.]+)\s*,?\s*center_y:\s*([-\dEe+.]+)\s*,?\s*radius:\s*([-\dEe+.]+)\s*,?\s*area:\s*([-\dEe+.]+)\s*,?\s*aspect_ratio:\s*([-\dEe+.]+)', re.IGNORECASE),
#     re.compile(r'Triangle\s+(\d+)\s*->\s*center_x:\s*([-\dEe+.]+)\s*center_y:\s*([-\dEe+.]+)\s*radius:\s*([-\dEe+.]+)\s*area:\s*([-\dEe+.]+)\s*aspect_ratio:\s*([-\dEe+.]+)', re.IGNORECASE),
# ]

# # ---------------------------
# # Parsers
# # ---------------------------
# def parse_file(path):
#     """
#     Retorna um dict com chaves 'nodes','edges','triangles', cada um é dict index->data
#     """
#     nodes = {}
#     edges = {}
#     triangles = {}

#     with open(path, 'r', encoding='utf-8', errors='ignore') as f:
#         for raw in f:
#             line = raw.rstrip('\n')
#             # try node patterns
#             for p in NODE_PATTERNS:
#                 m = p.search(line)
#                 if m:
#                     idx = int(m.group(1))
#                     x = ffloat(m.group(2))
#                     y = ffloat(m.group(3))
#                     z = ffloat(m.group(4))
#                     n_neighbor = int(m.group(5))
#                     boundary = m.group(6).strip().lower()
#                     boundary_bool = boundary in ('t', 'true')
#                     state = int(m.group(7))
#                     voronoi_area = ffloat(m.group(8))
#                     nodes[idx] = {
#                         'index': idx, 'x': x, 'y': y, 'z': z,
#                         'n_neighbor': n_neighbor, 'boundary': boundary_bool,
#                         'state': state, 'voronoi_area': voronoi_area,
#                         'raw': line
#                     }
#                     break
#             else:
#                 # try edges
#                 for p in EDGE_PATTERNS:
#                     m = p.search(line)
#                     if m:
#                         idx = int(m.group(1))
#                         eid = int(m.group(2))
#                         start = int(m.group(3))
#                         end = int(m.group(4))
#                         vx = ffloat(m.group(5))
#                         vy = ffloat(m.group(6))
#                         length = ffloat(m.group(7))
#                         slope = ffloat(m.group(8))
#                         boundary = m.group(9).strip().lower()
#                         boundary_bool = boundary in ('t', 'true')
#                         edges[idx] = {
#                             'index': idx, 'id': eid, 'start': start, 'end': end,
#                             'voronoi_x': vx, 'voronoi_y': vy, 'length': length,
#                             'slope': slope, 'boundary': boundary_bool,
#                             'raw': line
#                         }
#                         break
#                 else:
#                     # try triangle patterns
#                     for p in TRI_PATTERNS:
#                         m = p.search(line)
#                         if m:
#                             idx = int(m.group(1))
#                             cx = ffloat(m.group(2))
#                             cy = ffloat(m.group(3))
#                             radius = ffloat(m.group(4))
#                             area = ffloat(m.group(5))
#                             asp = ffloat(m.group(6))
#                             triangles[idx] = {
#                                 'index': idx, 'center_x': cx, 'center_y': cy,
#                                 'radius': radius, 'area': area, 'aspect_ratio': asp,
#                                 'raw': line
#                             }
#                             break
#     return {
#         'nodes': nodes,
#         'edges': edges,
#         'triangles': triangles
#     }

# # ---------------------------
# # Matching nodes by coordinates
# # ---------------------------
# def map_nodes_by_coords(for_nodes, rust_nodes, tol=1e-6):
#     """
#     Retorna dicionário map_for_to_rust: for_idx -> rust_idx
#     Procura por correspondência pelo par (x,y) mínimo; exige distância <= tol.
#     Se não encontrar tolerância, dá match pelo mais próximo e registra aviso.
#     """
#     rust_list = [(idx, d['x'], d['y']) for idx, d in rust_nodes.items()]
#     mapping = {}
#     warnings = []

#     for fidx, f in for_nodes.items():
#         fx, fy = f['x'], f['y']
#         best = None
#         bestd2 = None
#         for ridx, rx, ry in rust_list:
#             dx = fx - rx
#             dy = fy - ry
#             d2 = dx*dx + dy*dy
#             if best is None or d2 < bestd2:
#                 best = ridx
#                 bestd2 = d2
#         if best is None:
#             warnings.append(f'Node {fidx}: nenhum node rust encontrado')
#             continue
#         if bestd2 <= tol*tol:
#             mapping[fidx] = best
#         else:
#             # fallback: still map to nearest but note distance
#             mapping[fidx] = best
#             warnings.append(f'Node {fidx} mapped to Rust {best} (dist={math.sqrt(bestd2):.3e})')
#     return mapping, warnings

# # ---------------------------
# # Compare functions
# # ---------------------------
# def compare_nodes(for_nodes, rust_nodes, mapping, rtol=1e-6, atol=1e-12, max_show=20):
#     mismatches = []
#     total = 0
#     for fidx, f in for_nodes.items():
#         total += 1
#         ridx = mapping.get(fidx)
#         if ridx is None:
#             mismatches.append((fidx, 'no_rust_match', f, None))
#             continue
#         r = rust_nodes.get(ridx)
#         if r is None:
#             mismatches.append((fidx, 'rust_missing_index', f, None))
#             continue
#         # compare fields x,y,z,n_neighbor,boundary,voronoi_area
#         diffs = {}
#         for k in ('x','y','z','n_neighbor','boundary','voronoi_area'):
#             a = f.get(k)
#             b = r.get(k)
#             if isinstance(a, (float, int)) and isinstance(b, (float, int)):
#                 ok = is_close(float(a), float(b), rtol=rtol, atol=atol)
#             else:
#                 ok = (a == b)
#             if not ok:
#                 diffs[k] = (a, b)
#         if diffs:
#             mismatches.append((fidx, diffs, f, r))
#     return {'total': total, 'mismatches': mismatches, 'count_mismatch': len(mismatches)}

# def compare_edges(for_edges, rust_edges, node_map_for_to_rust, rtol=1e-6, atol=1e-12, max_show=20):
#     # build rust lookup by unordered endpoints and by id
#     rust_by_ends = {}
#     rust_by_id = {}
#     for ridx, e in rust_edges.items():
#         s, t = e['start'], e['end']
#         key = tuple(sorted((s, t)))
#         rust_by_ends[key] = e
#         rust_by_id[e['id']] = e

#     mismatches = []
#     total = 0
#     for fidx, fe in for_edges.items():
#         total += 1
#         # map endpoints using node_map
#         fs, fe_end = fe['start'], fe['end']
#         rs = node_map_for_to_rust.get(fs, fs)   # if no mapping, use original - fallback
#         rt = node_map_for_to_rust.get(fe_end, fe_end)
#         key = tuple(sorted((rs, rt)))
#         re_edge = rust_by_ends.get(key)
#         reason = None
#         if re_edge is None:
#             # try matching by id
#             re_edge = rust_by_id.get(fe.get('id'))
#             if re_edge is None:
#                 reason = f'no_rust_edge_with_ends_or_id (ends mapped {rs}-{rt}, id {fe.get("id")})'
#                 mismatches.append((fidx, reason, fe, None))
#                 continue
#         # compare numeric fields
#         diffs = {}
#         for k in ('voronoi_x','voronoi_y','length','slope'):
#             a = fe.get(k)
#             b = re_edge.get(k)
#             ok = is_close(float(a), float(b), rtol=rtol, atol=atol)
#             if not ok:
#                 diffs[k] = (a, b)
#         # boundary bool compare
#         if fe.get('boundary') != re_edge.get('boundary'):
#             diffs['boundary'] = (fe.get('boundary'), re_edge.get('boundary'))
#         if diffs:
#             mismatches.append((fidx, diffs, fe, re_edge))
#     return {'total': total, 'mismatches': mismatches, 'count_mismatch': len(mismatches)}

# def compare_triangles(for_tris, rust_tris, rtol=1e-6, atol=1e-12, max_show=20):
#     mismatches = []
#     total = 0
#     for fidx, ft in for_tris.items():
#         total += 1
#         rt = rust_tris.get(fidx)
#         if rt is None:
#             mismatches.append((fidx, 'no_rust_triangle_index', ft, None))
#             continue
#         diffs = {}
#         for k in ('center_x','center_y','radius','area','aspect_ratio'):
#             a = ft.get(k); b = rt.get(k)
#             if not is_close(float(a), float(b), rtol=rtol, atol=atol):
#                 diffs[k] = (a, b)
#         if diffs:
#             mismatches.append((fidx, diffs, ft, rt))
#     return {'total': total, 'mismatches': mismatches, 'count_mismatch': len(mismatches)}

# # ---------------------------
# # Orquestrador
# # ---------------------------
# def compare_files(fortran_path, rust_path, tol_node_match=1e-6, rtol=1e-6, atol=1e-12, max_show=20):
#     print(f'Parsing {fortran_path} (Fortran) and {rust_path} (Rust)...')
#     fdata = parse_file(fortran_path)
#     rdata = parse_file(rust_path)

#     print('Nodes: Fortran=%d  Rust=%d' % (len(fdata['nodes']), len(rdata['nodes'])))
#     print('Edges: Fortran=%d  Rust=%d' % (len(fdata['edges']), len(rdata['edges'])))
#     print('Triangles: Fortran=%d  Rust=%d' % (len(fdata['triangles']), len(rdata['triangles'])))

#     # map nodes (fortran -> rust) by coords
#     node_map, node_warnings = map_nodes_by_coords(fdata['nodes'], rdata['nodes'], tol=tol_node_match)
#     if node_warnings:
#         print('Node mapping WARNINGS (first 10):')
#         for w in node_warnings[:10]:
#             print('  -', w)

#     # compare nodes
#     cn = compare_nodes(fdata['nodes'], rdata['nodes'], node_map, rtol=rtol, atol=atol, max_show=max_show)
#     print('\nNodes compared:', cn['total'], 'mismatches:', cn['count_mismatch'])
#     if cn['count_mismatch']:
#         print('First node mismatches (up to %d):' % max_show)
#         for item in cn['mismatches'][:max_show]:
#             fidx, diffs, f, r = item
#             print(' Node', fidx, 'differences:')
#             if isinstance(diffs, str):
#                 print('   ', diffs)
#             else:
#                 for k, (a,b) in diffs.items():
#                     print('   ', k, 'Fortran=', a, 'Rust=', b)

#     # compare edges
#     ce = compare_edges(fdata['edges'], rdata['edges'], node_map, rtol=rtol, atol=atol, max_show=max_show)
#     print('\nEdges compared:', ce['total'], 'mismatches:', ce['count_mismatch'])
#     if ce['count_mismatch']:
#         print('First edge mismatches (up to %d):' % max_show)
#         for item in ce['mismatches'][:max_show]:
#             fidx, diffs, fe, re = item
#             print(' Edge', fidx, '->', diffs if isinstance(diffs,str) else '')
#             if isinstance(diffs, dict):
#                 for k,(a,b) in diffs.items():
#                     print('   ', k, 'Fortran=', a, 'Rust=', b)
#             if re is None:
#                 print('   Rust edge not found for this fortran edge.')

#     # compare triangles
#     ct = compare_triangles(fdata['triangles'], rdata['triangles'], rtol=rtol, atol=atol, max_show=max_show)
#     print('\nTriangles compared:', ct['total'], 'mismatches:', ct['count_mismatch'])
#     if ct['count_mismatch']:
#         print('First triangle mismatches (up to %d):' % max_show)
#         for item in ct['mismatches'][:max_show]:
#             fidx, diffs, ft, rt = item
#             print(' Triangle', fidx, 'differences:')
#             if isinstance(diffs, str):
#                 print('   ', diffs)
#             else:
#                 for k,(a,b) in diffs.items():
#                     print('   ', k, 'Fortran=', a, 'Rust=', b)

#     return {'nodes': cn, 'edges': ce, 'triangles': ct, 'node_map': node_map}


# # ---------------------------
# # CLI
# # ---------------------------
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Compare Fortran and Rust text outputs (Nodes, Edges, Triangles)')
#     parser.add_argument('fortran_file', help='Arquivo Fortran (.txt)')
#     parser.add_argument('rust_file', help='Arquivo Rust (.txt)')
#     parser.add_argument('--tol-node', type=float, default=1e-6, help='Tolerância para mapear nodes por coord (default 1e-6)')
#     parser.add_argument('--rtol', type=float, default=1e-6, help='rtol para comparações numéricas')
#     parser.add_argument('--atol', type=float, default=1e-12, help='atol para comparações numéricas')
#     args = parser.parse_args()

#     compare_files(args.fortran_file, args.rust_file, tol_node_match=args.tol_node, rtol=args.rtol, atol=args.atol)
