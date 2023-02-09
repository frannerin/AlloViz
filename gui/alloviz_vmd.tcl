
# https://www.tcl.tk/man/tcl8.4/TclCmd/vwait.html
# https://wiki.tcl-lang.org/page/The+simplest+possible+socket+demonstration

set alloviz_gui_dir [file dirname [file normalize [info script]]]
lappend auto_path $alloviz_gui_dir


package provide alloviz 1.0
package require json

namespace eval alloviz {
    variable already_registered 0

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

    proc check_vmd_topology_conformity {jlist} {
        set rl [::json::json2dict $jlist]
        set isok 1
        foreach rp $rl {
            lassign $rp rn ri
            set as [atomselect top "name CA and resid $ri and resname $rn"]
            set nm [$as num]
            if {$nm != 1} {
                set isok 0
                puts "Residue $ri resname $rn matched $nm times"
            }
        }
        return $isok
    }

}

::alloviz::register_menu
::alloviz::start 9990

