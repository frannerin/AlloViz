# vmdcollab
# 8-23-2000 by Justin Gullingsrud, based on echo server from Welch, 3rd ed.
#
# vmdcollab: let two or more VMD's broadcast their commands to each other.
# This allows, for example, mouse rotations in one VMD to be reflected in
# all the other VMD's.  
#
# Start the chat session with the command vmdcollab::start <host> <port>
# VMD will try to connect to a server on the given host and port.  Once
# connected, all commands issued by VMD that have text equivalents will be
# sent to the server, which is then expected to forward the text commands
# to the other VMD's (but _not_ back to the originator).  
# 
# Important trick: the chat client avoids endless bouncing of messages by
# turning off broadcasts while an incoming command is being evaluated.
# 
# The server in bounce.tcl can serve as the chat host, and can even run
# in the same Tcl interpreter as a chat client.

namespace eval vmdcollab {
  variable chatsock
  variable broadcast 

  proc start { host port } {
    variable chatsock 
    variable broadcast 

    set chatsock [socket $host $port]
    fconfigure $chatsock -buffering line
    fileevent $chatsock readable [list vmdcollab::recv $chatsock]
    set broadcast 1 

    # Assume that text will be placed in the global variable vmd_logfile
    uplevel #0 trace variable vmd_logfile w vmdcollab::send 
  }
 
  proc stop { } {
    variable chatsock
    puts "closing connection"
    close $chatsock
    trace vdelete vmd_logfile w vmdcollab::send
    return
  }
 
  proc recv { sock } {
    variable broadcast
    if { [eof $sock] || [catch {gets $sock line}]} {
      # end of file or abnormal connection drop
      puts "closing connection"
      close $sock
      trace vdelete vmd_logfile w vmdcollab::send
      return
    }
    # Turn off broadcast while evaluating, otherwise we would echo every
    # command we receive.
    set broadcast 0 
    eval $line 
    set broadcast 1
  }

  proc send { name1 name2 op } {
    variable broadcast
    variable chatsock

    if { $broadcast == 0 } { return }
    # Grab the text out of vmd_logfile and send it
    upvar #0 vmd_logfile line
    puts $chatsock $line
  }
}

  
