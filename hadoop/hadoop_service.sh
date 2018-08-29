

[program:hadoop-dfs]
command = /opt/soft/hadoop/sbin/start-dfs.sh
startsecs = 15
autorestart = true
startretries = 3
user = root
logfile = /tmp/dfs.log
pidfile = /tmp/dfs.pid
loglevel = debug
logfile_maxbytes = 50MB
logfile_backups=10

[program:hadoop-yarn]
command = /opt/soft/hadoop/sbin/start-yarn.sh
startsecs = 15
autorestart = true
startretries = 3
user = root
logfile = /tmp/yarn.log
pidfile = /tmp/yarn.pid
loglevel = debug
logfile_maxbytes = 50MB
logfile_backups=10

[program:hadoop]
command = /usr/local/bin/start_hadoop.sh
startsecs = 15
autorestart = true
startretries = 3
user = root
logfile = /tmp/hadoop.log
pidfile = /tmp/hadoop.pid
loglevel = trace # critical, error, warn, info, debug, trace, blather
logfile_maxbytes = 50MB
logfile_backups=10
# identifier = hadoop
# directory = /tmp
# nocleanup = true
# childlogdir = /tmp
