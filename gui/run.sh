#!/bin/bash

# Some dependency requires Homebrew's sqlite3, while brew does not link it 
# by default. Add it to the DLL path.

DIRECTORY=$(cd `dirname $0` && pwd)

# Looks like anaconda python's sqlite3 requires _sqlite3_enable_load_extension,
# which is only provided by the brew version of sqlite3 (not the system's).
# Errors is ImportError: dlopen(/.../lib/python3.10/lib-dynload/_sqlite3.cpython-310-darwin.so, 0x0002): Symbol not found: (_sqlite3_enable_load_extension)
if [[ $(uname) == "Darwin" ]]; then
    export DYLD_LIBRARY_PATH="/opt/homebrew/opt/sqlite3/lib" 
fi

python $DIRECTORY/alloviz_gui.py
