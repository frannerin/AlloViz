proc Server {channel clientaddr clientport} {
    set cmd [ gets $channel ]
    puts "Executing $cmd"
    set out [ eval  $cmd]
    puts "Returning $out"
    puts $channel $out
    close $channel
}

socket -server Server 9900
vwait forever
