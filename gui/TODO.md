# TODO and Bugs

This file tracks known issues and improvement ideas for the AlloViz GUI protocol and integration between Python and Tcl/VMD.

## Bugs

- `alloviz_gui.py`: In `sendVMDCommand`, the warning message references `self.HOST` and `self.PORT`, but the module-level constants are `_HOST` and `_PORT`. Use the correct names in the message.

## Improvements

- Safer command dispatch on the Tcl side: replace `eval $cmd` with a dispatcher that only invokes whitelisted `::alloviz::*` procedures, parsing arguments explicitly.

- Consistent return encoding: have `jsonwrap` JSON-encode return values so Python can reliably parse complex results (beyond short scalars/strings).

- Robust reply handling: support multi-chunk responses on the Python side (loop on `recv`) and/or adopt a length-prefixed or delimiter-based protocol for larger replies.

- Configurable port and better error reporting: allow configuring the TCP port and detect/report if the port is already in use; surface clear errors to the GUI.

- Retries and timeouts: add client-side timeouts and optional retry/backoff for transient connection errors.

- Argument quoting robustness: using Tcl braces to group JSON works for typical data, but JSON string values containing unmatched braces can interfere with Tcl bracing. Consider alternative encodings (e.g., base64) or escaping to make argument passing fully robust.

