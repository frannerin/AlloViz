
# https://www.tcl.tk/man/tcl8.4/TclCmd/vwait.html
# https://wiki.tcl-lang.org/page/The+simplest+possible+socket+demonstration

set here [file dirname [file normalize [info script]]]
lappend auto_path $here


package provide alloviz 1.0

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
        global here
        puts "Starting AlloViz GUI Python component $here/run.sh"
        exec $here/run.sh &
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
        $s writepsf $tmpbase.psf
        animate write dcd $tmpbase.dcd sel $s
        puts "::alloviz::dump_trajectory writing to $tmpbase.{psf,dcd}"
        return $tmpbase
    }

}

::alloviz::register_menu
::alloviz::start 9990

