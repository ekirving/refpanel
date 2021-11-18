#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# get the command line arguments
args=("$@")

datetime=$(date +%Y-%m-%d-%H%M.%S)

# make a unique log file
logfile="nohup-${datetime}.out"

refpanel=${args[0]:='test'}

# print the server name and start time to the log file
echo "SERVER: $HOSTNAME" >>${logfile}
echo "DATE: ${datetime}" >>${logfile}
echo "REFPANEL: ${refpanel}" >>$logfile

# load the conda environment
eval "$(conda shell.bash hook)"
conda activate refpanel

MAX_FTP=15

if ! command -v free &> /dev/null; then
  # MacOS does not have the free command
  MAX_MEM=$(sysctl -a | awk '/^hw.memsize:/{print $2/(1024)^2}')
else
  # but linux does
  MAX_MEM=$(free -m | awk '/^Mem:/{print $2}')
fi

flags="--cores all "
flags+="--nolock "
flags+="--keep-going "
flags+="--printshellcmds "
flags+="--show-failed-logs "
flags+="--rerun-incomplete "
flags+="--reason "
flags+="--use-conda "
flags+="--resources mem_mb=${MAX_MEM} ebi_ftp=${MAX_FTP} sanger_ftp=${MAX_FTP} "

(
  set -x
  snakemake ${flags} --config refpanel=${refpanel} -- refpanel
) &>>${logfile}

echo "DONE!" >>${logfile}
