
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

    proc readable {args} {
        variable guipipe
        puts "Readable with $args"
        set data [chan gets $guipipe]
        puts "Got $data"
        puts "Sending 2"
        chan puts $guipipe 2
    }

    proc writable {args} {
     puts "Writable with $args"
    }

    proc start port {
        puts "Starting server on port $port"
        socket -server accept $port
    }
    
    proc alloviz_tk {} {
        global here
        variable guipipe
        puts "Starting AlloViz GUI Python component $here/run.sh"
        # exec $here/run.sh &
        set guipipe [open "|$here/run.sh" r+]
        fconfigure $guipipe -blocking false -buffering line
        chan event $guipipe readable ::alloviz::readable
       # chan event $guipipe writable ::alloviz::writable
    }

    proc register_menu {} {
        variable already_registered
        if {$already_registered==0} {
            incr already_registered
            vmd_install_extension alloviz \
                alloviz::alloviz_tk \
                "Analysis/AlloViz Allostery Analysis GUI"
        }
    }
}

alloviz::register_menu
# alloviz::start 9990

