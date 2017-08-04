#!/bin/bash
date2stamp () {
    date --utc --date "$1" +%s
}
stamp2date (){
    date --utc --date "1970-01-01 $1 sec" "+%Y-%m-%d %T"
}

# USAGE # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# convert a date into a UNIX timestamp
    stamp=$(date2stamp "2006-10-01 15:00")
    echo $stamp

# from timestamp to date
    stamp2date $stamp


