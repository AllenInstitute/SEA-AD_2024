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

ChromA_count	L8AT_211011_01_B01-1136687632.tsv.gz	grps13-6.merged.bed	L8AT_211011_01_B01-1136687632.scATAC.cellbarcode.txt
ChromA_count	L8AT_211018_01_B03-1138433801.tsv.gz	grps13-6.merged.bed	L8AT_211018_01_B03-1138433801.scATAC.cellbarcode.txt
ChromA_count	L8AT_211018_01_C03-1138433802.tsv.gz	grps13-6.merged.bed	L8AT_211018_01_C03-1138433802.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_F05-1140028668.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_F05-1140028668.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_G05-1140028669.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_G05-1140028669.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_H05-1140028670.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_H05-1140028670.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_A06-1140028671.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_A06-1140028671.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_C07-1140030210.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_C07-1140030210.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_D07-1140030211.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_D07-1140030211.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_A07-1140030208.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_A07-1140030208.scATAC.cellbarcode.txt
ChromA_count	L8AT_211025_01_B07-1140030209.tsv.gz	grps13-6.merged.bed	L8AT_211025_01_B07-1140030209.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_B09-1140885932.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_B09-1140885932.scATAC.cellbarcode.txt
