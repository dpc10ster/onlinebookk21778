#!/usr/bin/expect --
##spawn ssh -l XuetongZhai 10.211.55.3
spawn ssh -l dev 172.16.36.232
##expect "login:"
##send "dev\r"
expect "password:"
send "1234\r"
expect ">"
send "net use X: /delete\r"
expect ">"
send "net use X: \"\\\\vmware-host\\Shared Folders\\Debug\"\r"
expect ">"
send "cd /D X:\r"
expect ">"
send "CORROC2.bat\r"
expect ">"
send "cd /D C:\r"
expect ">"
send "exit\r"
expect eof
