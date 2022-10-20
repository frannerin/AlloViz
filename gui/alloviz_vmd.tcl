
# https://www.tcl.tk/man/tcl8.4/TclCmd/vwait.html
# https://wiki.tcl-lang.org/page/The+simplest+possible+socket+demonstration

set here [file dirname [file normalize [info script]]]
lappend auto_path $here


package provide alloviz 1.0

namespace eval alloviz {
    variable already_registered 0

    proc accept {channel clientaddr clientport} {
        set cmd [ gets $channel ]
        puts "Executing $cmd"
        set out [ eval  $cmd]
        puts "Returning $out"
        puts $channel $out
        close $channel
    }

    proc start port {
        puts "Starting server on port $port"
        socket -server accept $port
    }
    

    proc alloviz_gui_start {} {
        global here
        puts "Starting AlloViz GUI Python component $here/run.sh"
        exec $here/run.sh &
    }

    proc alloviz_tk {} {
        alloviz_gui_start
    }

    proc register_menu {} {
        variable already_registered
        if {$already_registered==0} {
            incr already_registered
            vmd_install_extension alloviz \
                alloviz::alloviz_tk \
                "Analysis/AlloViz GUI"
        }
    }
}

alloviz::register_menu
alloviz::start 9990

