#PBS -l ncpus=1
#PBS -q longrun@waugh.cchem.berkeley.edu
#PBS -j oe
cd /usr/people/magneto/code/blochlib-1.0/examples/qmsims/perms

rectrain config_files
