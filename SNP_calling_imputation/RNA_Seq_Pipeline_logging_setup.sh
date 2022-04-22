#!/bin/bash

#######################

log_newstage(){
    echo [`date '+%y-%m-%d %H:%M:%S'`] INFO: Start job $1 > $DIR_log/$logfile
    # echo [`date '+%y-%m-%d %H:%M:%S'`] INFO: $* > $DIR_log/$logfile
}

log_info(){
    echo [`date '+%y-%m-%d %H:%M:%S'`] INFO: $1 >> $DIR_log/$logfile
}

log_success(){
    echo [`date '+%y-%m-%d %H:%M:%S'`] DONE: $1 >> $DIR_log/$logfile
}

log_error(){
    echo [`date '+%y-%m-%d %H:%M:%S'`] ERROR: $1 >> $DIR_log/$logfile
    echo [`date '+%y-%m-%d %H:%M:%S'`] Exit job. >> $DIR_log/$logfile
    exit 1
}
