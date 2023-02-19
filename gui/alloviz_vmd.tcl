
# https://www.tcl.tk/man/tcl8.4/TclCmd/vwait.html
# https://wiki.tcl-lang.org/page/The+simplest+possible+socket+demonstration

set alloviz_gui_dir [file dirname [file normalize [info script]]]
lappend auto_path $alloviz_gui_dir

# molinfo top 
# molinfo top get numreps
# Info) mol delrep 1 0
# Info) mol color Name
# Info) mol representation NewCartoon 0.300000 10.000000 4.100000 0
# Info) mol selection protein
# Info) mol material Opaque
# Info) mol addrep 0
# Info) mol modselect 1 0 all 

package provide alloviz 1.0
package require json

namespace eval alloviz {
    variable already_registered 0
    variable current_viz none

    proc accept {channel clientaddr clientport} {
        set cmd [ gets $channel ]
        puts "alloviz::accept <-- $cmd"
        set out [ eval  $cmd]
        puts "alloviz::accept --> $out"
        puts $channel $out
        close $channel
    }

    proc start port {
        puts "Starting server on port $port"
        socket -server ::alloviz::accept $port
    }
    
    proc alloviz_gui_start {} {
        global alloviz_gui_dir
        puts "Starting AlloViz GUI Python component $alloviz_gui_dir/run.sh"
        exec $alloviz_gui_dir/run.sh &
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

    proc visualize_nodes {asel rnl rvl} {
        foreach rn $rnl rv $rvl {
            set as [atomselect top "($asel) and resid $rn"]
            $as set beta $rv
            $as delete
        }
        set mn [molinfo top]
        set ri [molinfo top get numreps]
    }

}

::alloviz::register_menu
::alloviz::start 9990

