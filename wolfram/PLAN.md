# Plan: Polyhedra in Wolfram Language + WLJS

## TL;DR

Reimplement the polyhedra project using **Wolfram Language** running in **WLJS Notebook** (free, open-source). Wolfram provides rich native polyhedra support (`PolyhedronData`, `DualPolyhedron`, `TruncatedPolyhedron`, etc.) and WLJS gives interactive THREE.js-based 3D visualization with HTML export. Custom stellation and discovery will be implemented where Wolfram's built-ins fall short.

## Architecture

```
wolfram/
├── polyhedra.wln          ← Main WLJS notebook (exploration + visualization)
├── Polyhedra.wl           ← Reusable package: transformations, stellation, discovery
├── exports/               ← Generated HTML exports
└── README.md              ← Setup instructions + usage guide
```

- **Computation**: Wolfram Engine (free Community Edition)
- **Interface**: WLJS Notebook (Electron app, `.wln` format, git-friendly)
- **3D Rendering**: THREE.js via WLJS `Graphics3D` (interactive rotation/zoom/pan)
- **Output**: Self-contained interactive HTML files (no server needed to view)

## Phases

### Phase 1: Environment Setup

1. Install Wolfram Engine (`brew install --cask wolfram-engine` or download from wolfram.com)
2. Install WLJS Notebook (`brew install --cask wljs-notebook`)
3. Verify connection: open WLJS, evaluate `$Version` and `PolyhedronData["Icosahedron"]`
4. Create `wolfram/README.md` with setup instructions

### Phase 2: Core Package — Native Operations (*parallel with Phase 3*)

5. Create `wolfram/Polyhedra.wl` package with:
   - `polyhedronFromData[name]` — wrapper around `PolyhedronData` returning a standard internal representation: `<|"coords" -> ..., "faces" -> ..., "name" -> ...|>`
   - Leverage built-in transforms that map directly:
     - `dual[p]` → wraps `DualPolyhedron`
     - `truncate[p]` → wraps `TruncatedPolyhedron`
     - `augment[p]` → wraps `AugmentedPolyhedron`
     - `bevel[p]` → wraps `BeveledPolyhedron`
   - Custom transforms (port from R logic):
     - `quasi[p]` — rectification (edge midpoints become vertices)
     - `rhombic[p]` — cantellation (edges → squares)
     - `compose[p1, p2]` — compound of two polyhedra
   - `description[p]` — Schläfli symbol generation

6. Internal representation design:
   - Use `Association`: `<|"coords" -> {{x,y,z},...}, "faces" -> {{v1,v2,...},...}, "name" -> "...", "bodies" -> {...}|>`
   - Topology computed lazily via `topology[p]` (edge-to-face, vertex figures, etc.)

### Phase 3: Visualization Layer (*parallel with Phase 2*)

7. Create `drawPoly[p]` function that renders using `Graphics3D` + `GraphicsComplex`:
   - Per-face coloring (using `FaceForm` + color palette)
   - Edge rendering via `EdgeForm[{Thin, Black}]`
   - Labeling (polyhedron name + Schläfli description below)
8. Create `drawMultiple[{p1, p2, ...}]` — grid layout of multiple polyhedra (using `GraphicsRow` or `Grid`)
9. Verify interactive 3D works in WLJS (rotation, zoom, pan)
10. Test HTML export of a visualization cell

### Phase 4: Stellation & Discovery

11. Implement `stellate[p]` — general stellation:
    - For named stellations: look up in `PolyhedronData` first
    - For novel stellations: implement face-plane stellation algorithm (extend face planes, find intersections, construct new faces)
    - Fallback: try `ResourceFunction["PolyhedronStellate"]` if available
12. Implement `buildRegularPoly[coords, polygonSize, vertexSize, exampleEdge]` — fitting faces onto coordinates (port from R)
13. Implement `discover[p]` — grid search over edge lengths × polygon sizes × vertex sizes

### Phase 5: Extended Transformations (new capabilities beyond R)

14. Add `snub[p]` — snub transformation (chiral, not available in R)
15. Add `chamfer[p]` — chamfer operation (expand edges to hexagons)
16. Add Conway operator notation: `conway["tD"]` → truncated dual, etc.
17. Add `net[p]` — 2D unfolding/net (Wolfram has `PolyhedronData[..., "Net"]`)

### Phase 6: Main Notebook

18. Create `wolfram/polyhedra.wln` notebook with sections mirroring R notebook:
    - Platonic solids (all 5, from `PolyhedronData`)
    - Kepler-Poinsot solids (from coordinates or `PolyhedronData`)
    - Transformations: dual, truncate, quasi, rhombic, bevel, augment
    - Compounds: 5 tetrahedra, 5 cubes, 10 tetrahedra, etc.
    - Stellations (using custom `stellate`)
    - Discovery mode
    - Extended: snub, chamfer, Conway operators
    - Icosahedron stellations family
    - Catalan solids (duals of Archimedean)
19. Export complete notebook as self-contained interactive HTML

## Key Decisions

- **Pure WLJS** — no Python/Jupyter layer. Wolfram Language only.
- **Stellation strategy**: Use `PolyhedronData` for named stellations, implement general face-plane stellation algorithm for novel cases.
- **Internal representation**: Wolfram `Association` with coords/faces/name/bodies (analogous to R list structure).
- **Built-in transforms where possible**: `DualPolyhedron`, `TruncatedPolyhedron`, `AugmentedPolyhedron`, `BeveledPolyhedron` are native — wrap them rather than reimplementing.
- **Custom transforms for**: `quasi` (rectification), `rhombic` (cantellation), `compose`, `stellate`, `buildRegularPoly`, `discover`.
- **Scope excluded (for now)**: 2D paper layouts/nets (may add later via `PolyhedronData[..., "Net"]`), 4D polytopes.

## Verification Criteria

1. `PolyhedronData["Icosahedron"]` renders interactively in WLJS
2. `dual[icosahedron]` produces dodecahedron (compare vertex count: 20 → 12)
3. `truncate[cube]` produces truncated cube (14 faces: 8 triangles + 6 octagons)
4. Compound of 5 tetrahedra matches R output visually
5. HTML export of notebook opens in browser with interactive 3D intact
6. Custom `stellate[icosahedron]` produces small stellated dodecahedron
7. `discover[cube]` finds compound of two tetrahedra (as R version does)

## Further Considerations

- **Package loading in WLJS**: The `.wl` package can be loaded via `Get["Polyhedra.wl"]` or `Needs["Polyhedra`"]`. Need to verify WLJS handles package loading cleanly.
- **Performance for discovery**: Wolfram's compiled functions (`Compile`, `FunctionCompile`) can accelerate the grid search if needed — but try interpreted first.
- **Conway notation**: Could use `ResourceFunction["ConwayNotation"]` if it exists, otherwise implement a simple parser mapping single-letter operators to transform chains.

## Reference: R Transformations to Port

| R Function | Wolfram Equivalent | Notes |
|------------|-------------------|-------|
| `dual(p)` | `DualPolyhedron` | Native |
| `truncate(p)` | `TruncatedPolyhedron` | Native |
| `quasi(p)` | — | Custom: edge midpoints → vertices |
| `rhombic(p)` | — | Custom: cantellation |
| `compose(p1, p2)` | — | Custom: merge coords + faces |
| `buildRegularPoly(...)` | — | Custom: face-fitting algorithm |
| `discover(p)` | — | Custom: grid search |
| `description(p)` | — | Custom: Schläfli symbol |
| — | `AugmentedPolyhedron` | New (faces → pyramids) |
| — | `BeveledPolyhedron` | New |
| — | `snub` | New (chiral) |
| — | `chamfer` | New |
| — | `stellate` | New (face-plane stellation) |
