# -*- coding: utf-8 -*-
"""
washer_partition.py
-------------------
Around every circular THROUGH-HOLE of radius ~2 in a 3D solid plate:

  1. cut a concentric "washer" ring from the hole edge (r = 2) out to r = 3,
     THROUGH the full thickness (top face to bottom face);
  2. add radial "spider" splits so the ring is hex/structured mesh-able and
     has radial + through-thickness edges to seed against;
  3. seed the mesh: >= N_AROUND elements around each hole, N_RADIAL across
     the washer gap, N_THICK through the thickness.

How to run  (Abaqus PDE):  set the run target to KERNEL (or GUI), then Run.
    - KERNEL / GUI  -> connected to the model database (mdb). USE THESE.
    - LOCAL         -> standalone Python, no mdb -> this script will fail.
"""

import math
from abaqus import *
from abaqusConstants import *
import part

# ----------------------------------------------------------------------
# USER SETTINGS  -- edit these to match your model
# ----------------------------------------------------------------------
MODEL_NAME    = 'Model-1'   # name of the model in the mdb
PART_NAME     = 'Plate'     # name of the 3D solid part to partition
HOLE_RADIUS   = 2.0         # radius of the through-holes to look for
WASHER_RADIUS = 3.0         # OUTER radius of the washer ring
RTOL          = 0.05        # tolerance when matching the hole radius
CTOL          = 1.0e-4      # geometric tolerance

# Mesh options
ADD_RADIAL_SPLITS = True    # quadrant cuts -> structured mesh + radial edges
DO_SEEDING        = True    # apply the edge seeds below
N_AROUND          = 12      # >= this many elements AROUND each hole (total)
N_RADIAL          = 4       # elements ACROSS the washer gap (r=2 -> r=3)
N_THICK           = 2       # elements THROUGH the plate thickness

# Auto-detect holes; if it ever mis-fires set AUTO_DETECT=False and list
# the centres by hand -- that path is 100% reliable.
AUTO_DETECT    = True
MANUAL_CENTERS = [
    # (x, y),
]
# ----------------------------------------------------------------------

m = mdb.models[MODEL_NAME]
p = m.parts[PART_NAME]
EPS = 1.0e-3


def _flatten_pt(pt):
    """getCentroid() may return (x,y,z) or ((x,y,z),) -> normalise."""
    if len(pt) == 1 and hasattr(pt[0], '__len__'):
        return pt[0]
    return pt


def find_hole_centers(prt, radius, rtol, ctol):
    """Hole centre = centroid of the cylindrical WALL face (it lies on the
    hole axis). Wall faces have their centroid at mid-thickness, so we skip
    the flat top/bottom faces that share the same rim edges."""
    zvals = [v.pointOn[0][2] for v in prt.vertices]
    zmin, zmax = min(zvals), max(zvals)
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
        d = ((px - cx) ** 2 + (py - cy) ** 2) ** 0.5
        if abs(d - radius) > max(10 * rtol, 0.1):
            print('  WARNING: hole near (%.3f, %.3f) may be a split cylinder;'
                  ' centre uncertain -- use MANUAL_CENTERS if off.' % (cx, cy))
        if not any(abs(qx - cx) < ctol and abs(qy - cy) < ctol
                   for (qx, qy) in centers):
            centers.append((cx, cy))
    return centers


def washer_cells(prt, cx, cy, zmin, zmax, R):
    """The washer sub-cells = every cell fully inside the r=R cylinder around
    the hole (the surrounding plate cell extends past R and is excluded)."""
    return prt.cells.getByBoundingCylinder(center1=(cx, cy, zmax + EPS),
                                           center2=(cx, cy, zmin - EPS),
                                           radius=R + EPS)


def seed_hole(prt, cx, cy, zmin, zmax, Rw, n_around, n_radial, n_thick, ctol):
    """Classify the edges inside the washer cylinder and seed them:
       circular arcs -> around the hole, radial straights -> across the gap,
       vertical straights -> through the thickness."""
    eds = prt.edges.getByBoundingCylinder(center1=(cx, cy, zmax + EPS),
                                          center2=(cx, cy, zmin - EPS),
                                          radius=Rw + EPS)
    circ, radial, thick = [], [], []
    for e in eds:
        try:
            e.getRadius()
            circ.append(e)                 # circular arc / rim
            continue
        except Exception:
            pass
        vids = e.getVertices()
        if len(vids) < 2:
            continue
        p0 = prt.vertices[vids[0]].pointOn[0]
        p1 = prt.vertices[vids[-1]].pointOn[0]
        if abs(p1[2] - p0[2]) > ctol:
            thick.append(e)                # vertical -> through thickness
        else:
            radial.append(e)               # in-plane -> radial gap
    if circ:
        prt.seedEdgeByNumber(edges=circ, number=n_around, constraint=FINER)
    if radial:
        prt.seedEdgeByNumber(edges=radial, number=n_radial, constraint=FINER)
    if thick:
        prt.seedEdgeByNumber(edges=thick, number=n_thick, constraint=FINER)
    return len(circ), len(radial), len(thick)


# --- locate holes and plate extents -----------------------------------
zvals = [v.pointOn[0][2] for v in p.vertices]
zmin, zmax = min(zvals), max(zvals)
rmid = 0.5 * (HOLE_RADIUS + WASHER_RADIUS)

centers = find_hole_centers(p, HOLE_RADIUS, RTOL, CTOL) if AUTO_DETECT \
    else MANUAL_CENTERS
print('Washer partition: found %d hole(s).' % len(centers))
if not centers:
    raise ValueError('No holes of radius ~%g found -- check HOLE_RADIUS/RTOL '
                     'or use MANUAL_CENTERS.' % HOLE_RADIUS)

# datum plane on the top face (offset=zmax) + axes for sketch orientation
# (Y) and the through extrude (Z)
topdp = p.datums[p.DatumPlaneByPrincipalPlane(
    principalPlane=XYPLANE, offset=zmax).id]
yaxis = p.datums[p.DatumAxisByPrincipalAxis(principalAxis=YAXIS).id]
zaxis = p.datums[p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS).id]

# --- build the washer partitions --------------------------------------
for (cx, cy) in centers:
    # 1. imprint the washer circle on the top face. Sketching on a DATUM
    #    plane + PartitionCellBySketch is the step that ran reliably in your
    #    build; it creates the r=3 ring edge on the top face.
    t = p.MakeSketchTransform(sketchPlane=topdp, sketchUpEdge=yaxis,
                              sketchPlaneSide=SIDE1, origin=(cx, cy, zmax))
    s = m.ConstrainedSketch(name='__washer__', sheetSize=200.0, transform=t)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(WASHER_RADIUS, 0.0))
    p.PartitionCellBySketch(sketchPlane=topdp, sketchUpEdge=yaxis,
                            cells=p.cells, sketch=s)
    del m.sketches['__washer__']

    # 2. extrude that ring DOWN through the thickness (guarantees a real
    #    through-cut, top face to bottom face). Try both senses so it works
    #    regardless of how the part's normal is oriented.
    ring = p.edges.findAt(((cx + WASHER_RADIUS, cy, zmax),))
    try:
        p.PartitionCellByExtrudeEdge(line=zaxis, cells=p.cells,
                                     edges=ring, sense=REVERSE)
    except Exception:
        p.PartitionCellByExtrudeEdge(line=zaxis, cells=p.cells,
                                     edges=ring, sense=FORWARD)

    # 3. radial "spider" splits (two planes through the hole axis -> quads)
    if ADD_RADIAL_SPLITS:
        dpx = p.datums[p.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=cx).id]
        p.PartitionCellByDatumPlane(
            datumPlane=dpx,
            cells=washer_cells(p, cx, cy, zmin, zmax, WASHER_RADIUS))
        dpy = p.datums[p.DatumPlaneByPrincipalPlane(
            principalPlane=XZPLANE, offset=cy).id]
        p.PartitionCellByDatumPlane(
            datumPlane=dpy,
            cells=washer_cells(p, cx, cy, zmin, zmax, WASHER_RADIUS))

print('Washer partition: %d washer ring(s) cut through the thickness.'
      % len(centers))

# --- mesh seeding ------------------------------------------------------
if DO_SEEDING:
    # if the ring was split into 4 quadrant arcs, seed each arc by
    # ceil(N_AROUND/4) so the full circle still gets at least N_AROUND.
    n_arc = int(math.ceil(N_AROUND / 4.0)) if ADD_RADIAL_SPLITS else N_AROUND
    total_around = 4 * n_arc if ADD_RADIAL_SPLITS else n_arc
    for (cx, cy) in centers:
        try:
            nc, nr, nt = seed_hole(p, cx, cy, zmin, zmax, WASHER_RADIUS,
                                   n_arc, N_RADIAL, N_THICK, CTOL)
            print('  hole (%.3f, %.3f): ~%d elems around, %d radial edge(s), '
                  '%d thickness edge(s) seeded.'
                  % (cx, cy, total_around, nr, nt))
        except Exception as exc:
            print('  WARNING: seeding failed at (%.3f, %.3f): %s'
                  % (cx, cy, exc))
    print('Seeding done (>= %d elements around each hole).' % total_around)