#!/bin/csh -f
#  subset-allen-to-basal-ganglia.cmd
#
#  UGE job for subset-allen-to-basal-ganglia built Tue Aug  4 15:57:42 PDT 2015
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog
#$ -o /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.R
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=8192M,h_rt=8:00:00
#
#  Name of application for log
#$ -v QQAPP=R
#  Email address to notify
#$ -M dpolioud@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "R serial"
  set qqidir    = /u/home/d/dpolioud/a9-vivek-lipton-collab/code
  set qqjob     = subset-allen-to-basal-ganglia
  set qqodir    = /u/home/d/dpolioud/a9-vivek-lipton-collab/code
  cd     /u/home/d/dpolioud/a9-vivek-lipton-collab/code
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for subset-allen-to-basal-ganglia built Tue Aug  4 15:57:42 PDT 2015"
  echo ""
  echo "  subset-allen-to-basal-ganglia directory:"
  echo "    "/u/home/d/dpolioud/a9-vivek-lipton-collab/code
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "subset-allen-to-basal-ganglia started on:   "` hostname -s `
  echo "subset-allen-to-basal-ganglia started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load R
  module li
#
  echo ""
  echo R CMD BATCH  subset-allen-to-basal-ganglia.R subset-allen-to-basal-ganglia.out.$JOB_ID
#
  /usr/bin/time -p \
  R CMD BATCH  /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.R  /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.out.$JOB_ID 
#
  echo ""
  echo "subset-allen-to-basal-ganglia finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog.$JOB_ID  --------" >> /u/local/apps/queue.logs/R.log.serial
  if (`wc -l /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
        head -50 /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog.$JOB_ID >> /u/local/apps/queue.logs/R.log.serial
        echo " "  >> /u/local/apps/queue.logs/R.log.serial
        tail -10 /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog.$JOB_ID >> /u/local/apps/queue.logs/R.log.serial
  else
        cat /u/home/d/dpolioud/a9-vivek-lipton-collab/code/subset-allen-to-basal-ganglia.joblog.$JOB_ID >> /u/local/apps/queue.logs/R.log.serial
  endif
  exit (0)
