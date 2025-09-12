
# https://www.tcl.tk/man/tcl8.4/TclCmd/vwait.html
# https://wiki.tcl-lang.org/page/The+simplest+possible+socket+demonstration

set alloviz_gui_dir [file dirname [file normalize [info script]]]
lappend auto_path $alloviz_gui_dir
set ::env(ALLOVIZ_GUI_DIR) $alloviz_gui_dir


package provide alloviz 1.0
package require json

namespace eval alloviz {
    variable already_registered 0
    variable current_viz none
    variable gui_chan

    # Anonymous pipe handling: Python GUI writes requests to stdout; we read them
    # here and write replies back to its stdin over the same channel.
    proc on_gui_read {} {
        variable gui_chan
        if {[eof $gui_chan]} {
            puts "AlloViz GUI pipe closed"
            fileevent $gui_chan readable {}
            set rc [catch { close $gui_chan } msg]
            if {$rc} {
                puts stderr "Error closing AlloViz GUI pipe: $msg"
                if {[info exists ::errorInfo]} {
                    puts stderr $::errorInfo
                }
            }
            return
        }
        if {[gets $gui_chan cmd] < 0} {
            return
        }
        # puts "alloviz::pipe <-- $cmd"
        set out [ eval $cmd ]
        # puts "alloviz::pipe --> $out"
        puts $gui_chan $out
        flush $gui_chan
    }
    
    proc alloviz_gui_start {} {
        global alloviz_gui_dir
        variable gui_chan
        # Indicate to Python that pipe transport should be used.
        set ::env(ALLOVIZ_COMM) pipe
        puts "Starting AlloViz GUI Python component $alloviz_gui_dir/run_standalone.sh"
        # Open bidirectional pipeline to the GUI process
        set gui_chan [open "|$alloviz_gui_dir/run_standalone.sh" r+]
        fconfigure $gui_chan -buffering line -translation lf -blocking 0
        fileevent $gui_chan readable ::alloviz::on_gui_read
    }

    proc alloviz_tk {} {
        alloviz_gui_start
        return ""
    }

    proc register_menu {} {
        variable already_registered
        if {$already_registered==0} {
            incr already_registered
            vmd_install_extension alloviz \
                ::alloviz::alloviz_tk \
                "Analysis/AlloViz GUI"
        }
    }

    proc dump_trajectory {sel} {
        set tmpbase /var/tmp/alloviz_[pid]
        set s [atomselect top $sel]
        $s writepdb $tmpbase.pdb
        $s writepsf $tmpbase.psf
        animate write dcd $tmpbase.dcd sel $s
        puts "::alloviz::dump_trajectory writing to $tmpbase.{pdb,psf,dcd}"
        return $tmpbase
    }

    proc jsonwrap {fcn args} {
        # Call TCL function fcn with args after un-jsoning them
        set tclargs {}
        foreach ja $args {
            lappend tclargs [::json::json2dict $ja]
        }
        set r [$fcn {*}$tclargs]
        # TODO: json result
        return $r
    }

    proc check_vmd_topology_conformity {asel rlist} {
        set isok 1
        foreach rp $rlist {
            lassign $rp rn ri
            set as [atomselect top "($asel) and name CA and resid $ri"]
            set nm [$as num]
            $as delete
            if {$nm != 1} {
                set isok 0
                puts "Residue $ri resname $rn matched $nm times"
            }
        }
        return $isok
    }

    proc delete_current_viz {} {
        variable current_viz
        if {[lindex $current_viz 0] == "rep"} {
            mol delrep [lindex $current_viz 1] [lindex $current_viz 2]
        } elseif {[lindex $current_viz 0] == "mol"} {
            mol delete [lindex $current_viz 1]
        }
        set current_viz none
    }

# molinfo top 
# molinfo top get numreps
# Info) mol delrep 1 0
# Info) mol color Name
# Info) mol representation NewCartoon 0.300000 10.000000 4.100000 0
# Info) mol selection protein
# Info) mol material Opaque
# Info) mol addrep 0
# Info) mol modselect 1 0 all 

# Info) mol color Beta
# Info) mol representation VDW 1.000000 12.000000
# Info) mol selection protein and name CA
# Info) mol material Opaque
# Info) mol addrep 0

    proc visualize_nodes {asel rnl rvl} {
        variable current_viz
        foreach rn $rnl rv $rvl {
            set as [atomselect top "($asel) and resid $rn"]
            $as set beta $rv
            $as delete
        }
        delete_current_viz
        set mn [molinfo top]
        set ri [molinfo top get numreps]
        mol color Beta
        mol representation VDW 1 12
        mol selection "($asel) and name CA and resid $rnl"
        mol material Opaque
        mol addrep $mn
        set current_viz "rep $ri $mn"
    }

    proc visualize_edges {asel r1l r2l rvl} {
        variable current_viz

        set mtop [molinfo top]
        set mn [mol new]
        mol rename $mn AlloViz

        set rl [lsort -unique -integer [concat $r1l $r2l]]
        set coors [dict create]
        foreach r $rl {
            set as [atomselect $mtop "($asel) and resid $r and name CA"]
            set x [$as get {x y z}]
            dict set coors $r $x
            $as delete
        }

        set rvls [lsort -real $rvl]
        set vmin [lindex $rvls 0]
        set vmax [lindex $rvls end]
        set vrange [expr $vmax-$vmin]

        # colorscalebar.tcl
        set mincolorid [colorinfo num]
        set maxcolorid [expr [colorinfo max] - 1]
        set numscaleids [expr $maxcolorid - $mincolorid]

        set minthickness 0.0
        set maxthickness 1.0

        delete_current_viz
        # graphics $mn materials off
        # graphics $mn material Transparent
        foreach r1 $r1l r2 $r2l rv $rvl {
            #graphics $mn cylinder $x1 $x2 radius 1 filled yes
            if {$rv == 0.0} continue
            set rv_norm [expr ($rv-$vmin)/$vrange]
            set colid [expr int($rv_norm*$numscaleids+$mincolorid)]
            set thk [expr $rv_norm*($maxthickness-$minthickness)+$minthickness]
            graphics $mn color $colid
            # graphics $mn line {*}[dict get $coors $r1] {*}[dict get $coors $r2] width 1
            graphics $mn cylinder {*}[dict get $coors $r1] {*}[dict get $coors $r2] radius $thk filled yes
        }
        mol top $mtop
        display resetview ;# Why is it needed?

        set current_viz "mol $mn"
    }

}

::alloviz::register_menu
