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

ChromA_count	L8XR_210729_01_E09-1122543707.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210729_01_E09-1122543707.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210729_01_C09-1122543708.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210729_01_C09-1122543708.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210729_01_D09-1124629224.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210729_01_D09-1124629224.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210805_01_H09-1124629228.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210805_01_H09-1124629228.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210805_01_E10-1124629229.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210805_01_E10-1124629229.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210812_01_D11-1124987482.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210812_01_D11-1124987482.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210812_01_A12-1124987483.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210812_01_A12-1124987483.multiome.atac.cellbarcode.txt
