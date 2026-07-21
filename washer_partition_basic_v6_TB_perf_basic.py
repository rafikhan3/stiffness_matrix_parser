# -*- coding: utf-8 -*-
"""
hex_tubesheet_washer_partition.py
---------------------------------
Washer partitions for the HEXAGONAL TUBESHEET (~740 through-holes, r = 0.16").

Adapted from washer_partition.py for the tubesheet's non-equilateral
triangular hole array. Key differences (see the .md for the reasoning):

  * HOLE_RADIUS   = 0.16"                (Ø0.32" holes)
  * WASHER_RADIUS = 0.25"                capped by the tight B-C pitch 0.505"
                                         (half-pitch 0.2525" = overlap limit)
  * ~740 holes  -> the per-hole datum planes / spider splits from the original
                   script are dropped and the work is BATCHED into ~2 partition
                   features (one imprint, one extrude) so CAE stays responsive.
  * Edge holes  -> any hole whose full washer would cross the hexagon's free
                   edge is SKIPPED (reported), so no partitions fail.
  * Seed        = 12 elements around every partitioned hole.

Assumes the plate lies in the global XY plane with thickness along Z (holes
run along Z). If the run reports "found 0 hole(s)", the thickness axis is not
Z -- reorient the part or tell me and I'll switch the axis.

Run in the Abaqus PDE with the target set to KERNEL (or GUI). Not LOCAL.
"""

import math
from abaqus import *
from abaqusConstants import *
import part

# ----------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------
MODEL_NAME    = 'Model-1'   # name of the model in the mdb
PART_NAME     = 'Tubesheet' # name of the 3D solid tubesheet part
HOLE_RADIUS   = 0.16        # hole radius (Ø0.32")
WASHER_RADIUS = 0.25        # washer outer radius (< 0.2525" = half B-C pitch)
MIN_PITCH     = 0.505       # governing (smallest) hole pitch, B-C
N_AROUND      = 12          # elements around each hole circumference

RTOL          = 0.02        # tolerance when matching the hole radius
CTOL          = 1.0e-3      # tolerance for de-duplicating hole centres
BOUNDARY_SAMPLES = 12       # points sampled on each washer circle for the
                            # inside-the-part test (raise for safety, lower
                            # for speed)

AUTO_DETECT   = True
MANUAL_CENTERS = [
    # (x, y),
]
# ----------------------------------------------------------------------

m = mdb.models[MODEL_NAME]
p = m.parts[PART_NAME]
EPS = 1.0e-4


# ======================================================================
# Helpers reused from washer_partition.py
# ======================================================================
def _flatten_pt(pt):
    """getCentroid() may return (x,y,z) or ((x,y,z),) -> normalise."""
    if len(pt) == 1 and hasattr(pt[0], '__len__'):
        return pt[0]
    return pt


def find_hole_centers(prt, radius, rtol, ctol):
    """Hole centre = centroid of the cylindrical WALL face (on the hole
    axis). Wall faces have their centroid at mid-thickness, so the flat
    top/bottom faces that share the same rim edges are skipped."""
    zvals = [v.pointOn[0][2] for v in prt.vertices]
    zmin, zmax = min(zvals), max(zvals)
    sanity = 0.5 * radius
    centers = []
    for f in prt.faces:
        rim = None
        for eid in f.getEdges():
            e = prt.edges[eid]
            try:
                r = e.getRadius()
            except Exception:
                continue
            if abs(r - radius) <= rtol:
                rim = e
                break
        if rim is None:
            continue
        cx, cy, cz = _flatten_pt(f.getCentroid())
        if not (zmin + ctol < cz < zmax - ctol):
            continue
        px, py = rim.pointOn[0][0], rim.pointOn[0][1]
        if abs(((px - cx) ** 2 + (py - cy) ** 2) ** 0.5 - radius) > sanity:
            continue  # centroid not on the hole axis (split cylinder) -> skip
        if not any(abs(qx - cx) < ctol and abs(qy - cy) < ctol
                   for (qx, qy) in centers):
            centers.append((cx, cy))
    return centers


# ======================================================================
# Tubesheet-specific helpers
# ======================================================================
def on_material(prt, x, y, z):
    """True if (x, y, z) lies on a face of the part (i.e. on solid top
    surface, not over a hole void or off the part edge)."""
    try:
        return len(prt.faces.findAt(((x, y, z),))) > 0
    except Exception:
        return False


def washer_inside(prt, cx, cy, z, R, nsamp):
    """True if the whole washer circle of radius R around (cx, cy) stays on
    the (convex) top face -> the washer does not cross the hexagon edge."""
    for k in range(nsamp):
        a = 2.0 * math.pi * k / nsamp
        if not on_material(prt, cx + R * math.cos(a), cy + R * math.sin(a), z):
            return False
    return True


# ======================================================================
# Main
# ======================================================================
zvals = [v.pointOn[0][2] for v in p.vertices]
zmin, zmax = min(zvals), max(zvals)

if WASHER_RADIUS >= 0.5 * MIN_PITCH:
    print('WARNING: WASHER_RADIUS %.4f >= half of MIN_PITCH %.4f -> washers '
          'will OVERLAP on the tight pair. Reduce it.'
          % (WASHER_RADIUS, MIN_PITCH))

centers = find_hole_centers(p, HOLE_RADIUS, RTOL, CTOL) if AUTO_DETECT \
    else MANUAL_CENTERS
print('Tubesheet washer: found %d hole(s).' % len(centers))
if not centers:
    raise ValueError('No holes of radius ~%g found. If the part is fine, its '
                     'thickness axis is probably not Z -- reorient or adjust.'
                     % HOLE_RADIUS)

# --- split interior holes (washer fits) from edge holes (skip) --------
interior, skipped = [], []
for (cx, cy) in centers:
    if washer_inside(p, cx, cy, zmax, WASHER_RADIUS, BOUNDARY_SAMPLES):
        interior.append((cx, cy))
    else:
        skipped.append((cx, cy))
print('  %d interior hole(s) will get a washer; %d edge hole(s) skipped '
      '(washer would cross the boundary).' % (len(interior), len(skipped)))
if not interior:
    raise ValueError('Every hole was classified as an edge hole -- check '
                     'WASHER_RADIUS / BOUNDARY_SAMPLES.')

# --- datum plane on top + orientation/extrude axes --------------------
topdp = p.datums[p.DatumPlaneByPrincipalPlane(
    principalPlane=XYPLANE, offset=zmax).id]
yaxis = p.datums[p.DatumAxisByPrincipalAxis(principalAxis=YAXIS).id]
zaxis = p.datums[p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS).id]

# --- 1. ONE sketch with every interior washer circle -> ONE imprint ----
ext = max(max(abs(c[0]) for c in interior), max(abs(c[1]) for c in interior))
sheet = 2.0 * ext + 4.0 * WASHER_RADIUS + 1.0
t = p.MakeSketchTransform(sketchPlane=topdp, sketchUpEdge=yaxis,
                          sketchPlaneSide=SIDE1, origin=(0.0, 0.0, zmax))
s = m.ConstrainedSketch(name='__washers__', sheetSize=sheet, transform=t)
for (cx, cy) in interior:
    s.CircleByCenterPerimeter(center=(cx, cy), point1=(cx + WASHER_RADIUS, cy))
p.PartitionCellBySketch(sketchPlane=topdp, sketchUpEdge=yaxis,
                        cells=p.cells, sketch=s)
del m.sketches['__washers__']
print('  washer circles imprinted on the top face.')

# --- 2. collect every ring edge, extrude them all through the thickness -
rings, missing = [], 0
for (cx, cy) in interior:
    r = p.edges.findAt(((cx + WASHER_RADIUS, cy, zmax),))
    if len(r):
        rings.append(r[0])
    else:
        missing += 1
if missing > 0.2 * len(interior):
    print('  WARNING: %d of %d ring edges not found -- the batched sketch may '
          'be mirrored (local vs global axes). Check placement in the '
          'viewport.' % (missing, len(interior)))

ok = False
for sns in (REVERSE, FORWARD):
    try:
        p.PartitionCellByExtrudeEdge(line=zaxis, cells=p.cells,
                                     edges=rings, sense=sns)
        ok = True
        break
    except Exception:
        pass
if not ok:
    print('  Batched extrude failed; falling back to per-hole extrude...')
    fails = 0
    for (cx, cy) in interior:
        r = p.edges.findAt(((cx + WASHER_RADIUS, cy, zmax),))
        if not len(r):
            fails += 1
            continue
        done = False
        for sns in (REVERSE, FORWARD):
            try:
                p.PartitionCellByExtrudeEdge(line=zaxis, cells=p.cells,
                                             edges=(r[0],), sense=sns)
                done = True
                break
            except Exception:
                pass
        if not done:
            fails += 1
    print('  per-hole extrude finished, %d failure(s).' % fails)
print('Tubesheet washer: %d washer ring(s) cut through the thickness.'
      % len(interior))

# --- 3. seed 12 elements around each partitioned hole ------------------
seed_edges = []
for (cx, cy) in interior:
    for (rr, zz) in ((HOLE_RADIUS, zmax), (HOLE_RADIUS, zmin),
                     (WASHER_RADIUS, zmax), (WASHER_RADIUS, zmin)):
        r = p.edges.findAt(((cx + rr, cy, zz),))
        if len(r):
            seed_edges.append(r[0])
if seed_edges:
    p.seedEdgeByNumber(edges=seed_edges, number=N_AROUND, constraint=FINER)
print('Seeding done: %d element(s) around each of %d hole(s).'
      % (N_AROUND, len(interior)))

if skipped:
    print('Skipped edge holes (no washer), first few centres:')
    for (cx, cy) in skipped[:10]:
        print('    (%.4f, %.4f)' % (cx, cy))
    if len(skipped) > 10:
        print('    ... and %d more.' % (len(skipped) - 10))