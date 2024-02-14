#! /usr/bin/env python

"""\
This script is responsible to submit jobs to the grid
Given a directory, it looks for existing genie gst.root files
and computes the radiative corrections for each root file

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
"""
import os, optparse, glob, tarfile, re

op = optparse.OptionParser(usage=__doc__)
op.add_option("--git-location", dest="GIT_LOCATION", default="https://github.com/e4nu/e4nuanalysiscode.git", help="e4nu code location in github. Defaulted to %default")
op.add_option("--git-branch", dest="BRANCH", default="master", help="Branch name. Default: %default")
op.add_option("--directory", dest="JOBSTD", default=os.getenv('PWD'), help="Directory where the unradiated gst files are located (default: %default)")
op.add_option("--ebeam-energy", dest="EnergyBeam", default="2", help="Comma separated list of beam energy for electrons. Default %default GeV")
op.add_option("--model", dest="MODEL", default="simple", help="Rad corr model to use. Default %default")
op.add_option("--target", dest="TARGET", default=1000010010, help="Target used for calculation. Default %default")
op.add_option("--thickness", dest="THICKNESS", default="0", help="Thickness. Defaulted to CLAS6 values")
op.add_option("--Emax", dest="EMAX", default="0.2",help="Maximum energy for your histogram axis. Default %default")

opts, args = op.parse_args()

# Check jobstopdir exists
if not os.path.exists(opts.JOBSTD) :
    print ( "Jobs top dir path "+opts.JOBSTD+" doesn't exist. Abort...\n")
    exit()

#JobSub is made available through the UPS package jobsub_client
os.system("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup" ) 

if opts.BRANCH: 
    print( ' Cloning e4nucode ' + opts.BRANCH ) 

# Configure grid
E4NUCODE=os.getenv('E4NUANALYSIS')
e4nu_setup_file = E4NUCODE+'/e4nu_gpvm_env.sh'
e4nu_pnfs_setup = opts.JOBSTD+"/e4nu_gpvm_env.sh"
if os.path.exists(e4nu_pnfs_setup) : 
    os.remove(e4nu_pnfs_setup)

# Configure additional optional options
options = "" 
#if( opts.THICKNESS ) options_str = options + " --thickness "+str(opts.THICKNESS)
#if( opts.EMAX ) options_str = options + " --Emax "+str(opts.EMAX)

os.system('cp '+e4nu_setup_file+' '+opts.JOBSTD )

rad_dir = opts.JOBSTD+"/radcorr/"
if not os.path.exists(rad_dir) : 
    os.mkdir(rad_dir)

counter = 0 
gst_file_names = glob.glob(opts.JOBSTD+"/*.gst.root")
name_out_file = "rad_corr"

# xml script
if os.path.exists(rad_dir+"grid_submission.xml") : 
    os.remove(rad_dir+"grid_submission.xml")
grid = open( rad_dir+"grid_submission.xml", 'w' ) 
grid.write("<parallel>\n")

for x in gst_file_names:
    genie_file=os.path.basename(gst_file_names[counter])
    
    number = re.findall(r'\d+',gst_file_names[counter])[-1:][0]
    
    if os.path.exists(rad_dir+name_out_file+"_"+str(number)+".sh"):
        os.remove(rad_dir+name_out_file+"_"+str(number)+".sh")
    script = open( rad_dir+name_out_file+"_"+str(number)+".sh", 'w' ) 

    script.write("#!/bin/bash \n")
    script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups \n")
    script.write("setup ifdhc v2_6_6 \n")
    script.write("export IFDH_CP_MAXRETRIES=0 \n")
    script.write("setup pdfsets v5_9_1b \n")
    script.write("setup gdb v8_1 \n")
    script.write("cd $CONDOR_DIR_INPUT \n")
    script.write("ifdh cp -D "+gst_file_names[counter]+" $CONDOR_DIR_INPUT/ ;\n \n")
    script.write("git clone "+opts.GIT_LOCATION+" -b "+opts.BRANCH+" ;\n")
    script.write("cd e4nuanalysiscode ; source e4nu_gpvm_env.sh ; make ;\n")
    #write main command
    script.write("./process_radweights --input-gst-file $CONDOR_DIR_INPUT/"+genie_file+" --output-gst-file $CONDOR_DIR_INPUT/rad_corr_"+str(number)+".gst.root --true-EBeam "+str(opts.EnergyBeam)+" --target "+str(opts.TARGET)+" --rad-model "+opts.MODEL+options+" ; \n\n")
    script.write("ifdh cp -D $CONDOR_DIR_INPUT/rad_corr_"+str(number)+".gst.root "+rad_dir+" \n")

    grid.write("jobsub_submit  -n --memory=4GB --disk=4GB --expected-lifetime=4h  --OS=SL7 --mail_on_error file://"+rad_dir+name_out_file+"_"+str(number)+".sh \n")

    counter += 1

grid.write("</parallel>\n")

if os.path.exists(rad_dir+"fnal_dag_submit.fnal"):
    os.remove(rad_dir+"fnal_dag_submit.fnal")

fnal_script = open( rad_dir+"fnal_dag_submit.fnal", 'w' ) 
fnal_script.write("#!/bin/bash \n")
fnal_script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups ;\n")
fnal_script.write("setup fife_utils ;\n")
fnal_script.write("jobsub_submit -G genie --OS=SL7 --memory=10GB --disk=10GB --expected-lifetime=25h -N 1 --role=Analysis --dag file://"+rad_dir+"/grid_submission.xml;\n")
