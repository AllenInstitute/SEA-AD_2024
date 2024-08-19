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

ChromA_count	L8AT_210804_01_B06-1123939344.tsv.gz	grps13-6.merged.bed	L8AT_210804_01_B06-1123939344.scATAC.cellbarcode.txt
ChromA_count	L8AT_210809_01_D07-1126220030.tsv.gz	grps13-6.merged.bed	L8AT_210809_01_D07-1126220030.scATAC.cellbarcode.txt
ChromA_count	L8AT_210809_01_E07-1126220031.tsv.gz	grps13-6.merged.bed	L8AT_210809_01_E07-1126220031.scATAC.cellbarcode.txt
ChromA_count	L8AT_210816_01_B09-1126220032.tsv.gz	grps13-6.merged.bed	L8AT_210816_01_B09-1126220032.scATAC.cellbarcode.txt
ChromA_count	L8AT_210816_01_C09-1126220033.tsv.gz	grps13-6.merged.bed	L8AT_210816_01_C09-1126220033.scATAC.cellbarcode.txt
ChromA_count	L8AT_210824_01_D12-1129176225.tsv.gz	grps13-6.merged.bed	L8AT_210824_01_D12-1129176225.scATAC.cellbarcode.txt
ChromA_count	L8AT_210824_01_E12-1129176226.tsv.gz	grps13-6.merged.bed	L8AT_210824_01_E12-1129176226.scATAC.cellbarcode.txt
ChromA_count	L8AT_210830_01_H03-1131257166.tsv.gz	grps13-6.merged.bed	L8AT_210830_01_H03-1131257166.scATAC.cellbarcode.txt
ChromA_count	L8AT_210916_01_B08-1131409693.tsv.gz	grps13-6.merged.bed	L8AT_210916_01_B08-1131409693.scATAC.cellbarcode.txt
ChromA_count	L8AT_211004_01_G10-1135042880.tsv.gz	grps13-6.merged.bed	L8AT_211004_01_G10-1135042880.scATAC.cellbarcode.txt
ChromA_count	L8AT_211004_01_H10-1135042881.tsv.gz	grps13-6.merged.bed	L8AT_211004_01_H10-1135042881.scATAC.cellbarcode.txt
ChromA_count	L8AT_211011_01_A01-1136687631.tsv.gz	grps13-6.merged.bed	L8AT_211011_01_A01-1136687631.scATAC.cellbarcode.txt
