#!/bin/bash

ALLEN_WG=/allen/programs/celltypes/workgroups/
ALLEN_PROD=/allen/programs/celltypes/production/
CHROMADIR=/allen/programs/celltypes/workgroups/hct/SEA-AD/environments

function ChromA_count() {
echo Running Chroma on i=$1 bed=$2 and c=$3
sbatch --ntasks-per-node=1 --cpus-per-task=32 --mem=250G -t 02:00:00 --qos=normal -p celltypes -J $2 -e $2.err --mail-type=FAIL <<END
#!/bin/bash

mkdir \$TMPDIR/{var,tmp}
export OWD=`pwd`

singularity exec --bind="\$TMPDIR/var/:/var" --bind="\$TMPDIR/tmp:/tmp" --bind "/sys/fs/cgroup/:/sys/fs/cgroup/" --bind "$ALLEN_WG:$ALLEN_WG" --bind "$ALLEN_PROD:$ALLEN_PROD" \
--contain --workdir \$TMPDIR/tmp $CHROMADIR/ChromA_biggermem.sif $CHROMADIR/ChromA_count_boilerplate.sh "$1" $2 $3 \$OWD

exit 0

END
}

ChromA_count	L8AT_220105_01_F11-1152530383.tsv.gz	grps13-6.merged.bed	L8AT_220105_01_F11-1152530383.scATAC.cellbarcode.txt
ChromA_count	L8AT_211222_01_C01-1152618570.tsv.gz	grps13-6.merged.bed	L8AT_211222_01_C01-1152618570.scATAC.cellbarcode.txt
ChromA_count	L8AT_211222_01_D01-1152618571.tsv.gz	grps13-6.merged.bed	L8AT_211222_01_D01-1152618571.scATAC.cellbarcode.txt
ChromA_count	L8AT_211222_01_E01-1152618572.tsv.gz	grps13-6.merged.bed	L8AT_211222_01_E01-1152618572.scATAC.cellbarcode.txt
ChromA_count	L8AT_211222_01_B01-1152618568.tsv.gz	grps13-6.merged.bed	L8AT_211222_01_B01-1152618568.scATAC.cellbarcode.txt
ChromA_count	L8AT_220105_01_G11-1153896887.tsv.gz	grps13-6.merged.bed	L8AT_220105_01_G11-1153896887.scATAC.cellbarcode.txt
ChromA_count	L8AT_220105_01_H11-1153896888.tsv.gz	grps13-6.merged.bed	L8AT_220105_01_H11-1153896888.scATAC.cellbarcode.txt
ChromA_count	L8AT_220105_01_A12-1153896889.tsv.gz	grps13-6.merged.bed	L8AT_220105_01_A12-1153896889.scATAC.cellbarcode.txt
ChromA_count	L8AT_220112_01_F03-1153896349.tsv.gz	grps13-6.merged.bed	L8AT_220112_01_F03-1153896349.scATAC.cellbarcode.txt
ChromA_count	L8AT_220112_01_G03-1153896350.tsv.gz	grps13-6.merged.bed	L8AT_220112_01_G03-1153896350.scATAC.cellbarcode.txt
ChromA_count	L8AT_220112_01_H03-1153896351.tsv.gz	grps13-6.merged.bed	L8AT_220112_01_H03-1153896351.scATAC.cellbarcode.txt
ChromA_count	L8AT_220112_01_A04-1155512003.tsv.gz	grps13-6.merged.bed	L8AT_220112_01_A04-1155512003.scATAC.cellbarcode.txt
