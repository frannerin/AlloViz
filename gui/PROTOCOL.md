# AlloViz GUI Protocol

This document describes the interface between the Python GUI and the Tcl/VMD plugin used by AlloViz.

## Components and Roles

- Python GUI (`alloviz_gui.py`)
  - PyQt5 application that computes/filters/analyzes data using the `AlloViz` Python package.
  - Communicates with VMD via a TCP socket to drive visualizations.

- Tcl/VMD plugin (`alloviz_vmd.tcl`)
  - Runs inside VMD, starts a TCP server, exposes functions under the `::alloviz::` namespace.
  - Executes visualization actions (atoms/edges) and can dump trajectories.
  - Spawns the Python GUI and sets environment variables the GUI relies on.

## Environment and Paths

- The Tcl plugin sets `::env(ALLOVIZ_GUI_DIR)` to the directory of the GUI files and appends it to `auto_path`.
- The Python GUI imports UI modules using `os.environ["ALLOVIZ_GUI_DIR"]`.

## Network Channel

- Host: `localhost`
- Port: `9990` (started by `::alloviz::start 9990`)
- Transport: TCP; one request per connection.
- Framing: single line terminated by `\n` per request and single line reply per response.

Client flow (Python):
1. Open a socket to `localhost:9990`.
2. Send one line with the Tcl command to execute.
3. Read the reply string.
4. Close the socket.

Server flow (Tcl/VMD):
1. Accept a connection.
2. Read a line with `gets`.
3. Evaluate the command and capture its return value.
4. Write the return value and close the connection.

## Command Framing and Argument Passing

- Python builds commands via `doVMDcall(fcn, *args)` which generates:
  - `::alloviz::jsonwrap <fcn> {<json-arg-1>} {<json-arg-2>} ...`
  - Each Python argument is JSON-serialized and wrapped in Tcl braces to keep it as a single token.
- The Tcl side converts each JSON argument into Tcl data structures using `json::json2dict` and then calls the target function.

## Supported Functions and Semantics

All functions live under the `::alloviz::` namespace and are invoked through `jsonwrap` unless otherwise noted.

### dump_trajectory

- Signature: `::alloviz::dump_trajectory asel`
- Args:
  - `asel` (string): VMD atom selection to dump.
- Behavior: Writes PDB/PSF/DCD to `/var/tmp/alloviz_<pid>.{pdb,psf,dcd}` using the selection.
- Returns: Base path string (without extension), e.g., `/var/tmp/alloviz_12345`.

### check_vmd_topology_conformity

- Signature: `::alloviz::check_vmd_topology_conformity asel rlist`
- Args:
  - `asel` (string): VMD atom selection.
  - `rlist` (list): List of `[resname, resid]` pairs to check.
- Behavior: For each pair, selects `($asel) and name CA and resid <resid>` and ensures it matches exactly one atom.
- Returns: `1` if all unique, `0` otherwise.

### visualize_nodes

- Signature: `::alloviz::visualize_nodes asel rnl rvl`
- Args:
  - `asel` (string): VMD atom selection.
  - `rnl` (list): Residue IDs (integers/strings acceptable to `resid`).
  - `rvl` (list): Values (floats) associated with each residue.
- Behavior:
  - Sets the `beta` field for residues in `rnl` to the corresponding values in `rvl`.
  - Removes the previous visualization, then creates a VDW representation on CA atoms colored by `Beta`.
- Returns: Empty string.

### visualize_edges

- Signature: `::alloviz::visualize_edges asel r1l r2l rvl`
- Args:
  - `asel` (string): VMD atom selection.
  - `r1l` (list): First residues of each edge.
  - `r2l` (list): Second residues of each edge.
  - `rvl` (list): Edge weights (floats).
- Behavior:
  - Builds a new molecule named `AlloViz`.
  - For all unique residues in `r1l ∪ r2l`, fetches CA coordinates from the current top molecule.
  - Draws cylinders between residue pairs with color and thickness proportional to the edge weight.
  - Restores the original top molecule and resets the view.
- Returns: Empty string.

## Python Client Usage Pattern

High-level sequence used by the GUI (`AlloVizWindow.runAnalysis`):
1. Dump trajectory via `::alloviz::dump_trajectory {<asel>}` and receive the base path.
2. Load `AlloViz.Protein` with the dumped `pdb/psf/dcd` and compute the selected `method`.
3. Optionally compute `GetContacts` if requested by the user.
4. Filter results with selected filters and analyze for chosen `elements` (`nodes` or `edges`) and `metrics` (`btw`, `cfb`, or `raw`).
5. Optionally validate selection consistency via `::alloviz::check_vmd_topology_conformity`.
6. Select a subset above a threshold and visualize:
   - Nodes: `::alloviz::visualize_nodes {asel} {rnl} {rvl}`.
   - Edges: `::alloviz::visualize_edges {asel} {r1l} {r2l} {rvl}`.

## Tcl/VMD Plugin Behavior

- Server startup: `::alloviz::start 9990` initializes the TCP server.
- GUI menu: `::alloviz::register_menu` installs a VMD menu entry under `Analysis/AlloViz GUI` that launches the GUI.
- GUI launcher: `::alloviz::alloviz_gui_start` executes `run_standalone.sh` in the GUI directory.
- Visualization lifecycle: `::alloviz::delete_current_viz` ensures only one active visualization (representation or molecule) at a time.

## How to Run

- From VMD: source the plugin and start the server (menu entry is also installed):
  - `source gui/alloviz_vmd.tcl`
  - or run `vmd -e gui/run_vmd.tcl` to also load sample data and launch the GUI automatically.
- The GUI launcher (`run_standalone.sh`) sets `DYLD_LIBRARY_PATH` for macOS (Homebrew sqlite3) and runs `python alloviz_gui.py`.

## Data Flow Summary

GUI (Python) → Server (Tcl/VMD):
- Requests: single-line commands invoking `::alloviz::jsonwrap` with JSON-encoded args.
- Typical: dump trajectory, check topology, visualize nodes/edges.

Server (Tcl/VMD) → GUI (Python):
- Replies: a single line string per request (e.g., base path, `1`/`0`, or empty string).

