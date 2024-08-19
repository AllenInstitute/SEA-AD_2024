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

ChromA_count	L8AT_220112_01_B04-1155512005.tsv.gz	grps13-6.merged.bed	L8AT_220112_01_B04-1155512005.scATAC.cellbarcode.txt
ChromA_count	L8AT_220126_01_A09-1157000752.tsv.gz	grps13-6.merged.bed	L8AT_220126_01_A09-1157000752.scATAC.cellbarcode.txt
ChromA_count	L4AT_201012_01_B04-1064293374.tsv.gz	grps13-6.merged.bed	L4AT_201012_01_B04-1064293374.scATAC.cellbarcode.txt
ChromA_count	L4AT_201012_01_C04-1064293376.tsv.gz	grps13-6.merged.bed	L4AT_201012_01_C04-1064293376.scATAC.cellbarcode.txt
ChromA_count	L4AT_201012_01_D04-1064293377.tsv.gz	grps13-6.merged.bed	L4AT_201012_01_D04-1064293377.scATAC.cellbarcode.txt
ChromA_count	L4AT_201019_01_E04-1064293378.tsv.gz	grps13-6.merged.bed	L4AT_201019_01_E04-1064293378.scATAC.cellbarcode.txt
ChromA_count	L4AT_201019_01_F04-1064293379.tsv.gz	grps13-6.merged.bed	L4AT_201019_01_F04-1064293379.scATAC.cellbarcode.txt
ChromA_count	L8AT_201028_01_H04-1074937000.tsv.gz	grps13-6.merged.bed	L8AT_201028_01_H04-1074937000.scATAC.cellbarcode.txt
ChromA_count	L8AT_201028_01_A05-1074937004.tsv.gz	grps13-6.merged.bed	L8AT_201028_01_A05-1074937004.scATAC.cellbarcode.txt
ChromA_count	L8AT_201028_01_B05-1074937005.tsv.gz	grps13-6.merged.bed	L8AT_201028_01_B05-1074937005.scATAC.cellbarcode.txt
ChromA_count	L8AT_201028_01_C05-1074937006.tsv.gz	grps13-6.merged.bed	L8AT_201028_01_C05-1074937006.scATAC.cellbarcode.txt
ChromA_count	L4AT_201019_01_G04-1075757819.tsv.gz	grps13-6.merged.bed	L4AT_201019_01_G04-1075757819.scATAC.cellbarcode.txt
