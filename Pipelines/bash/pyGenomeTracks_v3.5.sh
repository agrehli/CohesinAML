#!/usr/bin/env bash

# (c) C Kohler -- 2020/07 
#
# simple wrapper shell script for pyGenomeTracks, 
# automagically picking the correct build for the host OS.
# Exits, if build for host OS is not available!
# 
# Important: we need to distinguish between the virtual python environment
# as created on Debian 9 and on Debian 10 (as the virtual environment on
# Debian 9 does not run on Debian 10) !
#

DEBIAN_VERS=$(cat /etc/debian_version | awk -F. '{print $1}')

_desc=$(lsb_release -d | awk '{print $2"" $3" " $4" "$5}')

DIR_PGT="/misc/software/ngs/pyGenomeTracks"
VERS_PGT="v3.5"
EXEC_PGT="${DIR_PGT}/${VERS_PGT}/debian${DEBIAN_VERS}_build/py3Env/bin/pyGenomeTracks"

if [ ! -e ${EXEC_PGT} ]; then
    printf "\n"
    printf "No pyGenomeTracks installation available for: ${_desc}\n"
    printf "\n"
    exit;
fi

exec "${EXEC_PGT}" "$@"
