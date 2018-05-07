#!/usr/bin/expect --
spawn ssh -l XuetongZhai 10.211.55.3
expect "password:"
send "1234\r"
expect ">"
send "net use X: /delete\r"
expect ">"
send "net use X: \"\\\\Mac\\Debug\"\r"
expect ">"
send "cd /D X:\r"
expect ">"
send "CORROC2.bat\r"
expect ">"
send "cd /D C:\r"
expect ">"
send "exit\r"
expect eof
