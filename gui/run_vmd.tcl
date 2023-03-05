source gui/alloviz_vmd.tcl

# Load sample data if present
catch {
    mol new 117/11160_dyn_117.psf
    mol addfile 117/11156_trj_117_r.dcd waitfor all
}

# Simulate the menu action
alloviz::alloviz_tk

