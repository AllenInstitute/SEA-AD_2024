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

ChromA_count	L8XR_210812_01_F11-1124987484.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210812_01_F11-1124987484.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210812_01_E11-1126220036.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210812_01_E11-1126220036.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210826_01_G03-1129238386.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210826_01_G03-1129238386.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210902_02_B08-1131257169.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210902_02_B08-1131257169.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210902_02_D08-1131257170.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210902_02_D08-1131257170.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210916_02_B11-1131607954.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210916_02_B11-1131607954.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210930_02_C03-1135448412.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210930_02_C03-1135448412.multiome.atac.cellbarcode.txt
