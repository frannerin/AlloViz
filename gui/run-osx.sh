#!/bin/bash

# Some dependency requires Homebrew's sqlite3, while brew does not link it 
# by default. Add it to the DLL path.

DIRECTORY=$(cd `dirname $0` && pwd)

export DYLD_LIBRARY_PATH="/opt/homebrew/opt/sqlite3/lib" 
python $DIRECTORY/alloviz_gui.py
