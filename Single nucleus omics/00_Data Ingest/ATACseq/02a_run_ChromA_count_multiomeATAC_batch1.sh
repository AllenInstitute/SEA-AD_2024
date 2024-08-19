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

ChromA_count	L8XR_210715_01_D12-1121747486.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210715_01_D12-1121747486.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210715_01_E12-1121747487.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210715_01_E12-1121747487.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210715_01_A12-1121939868.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210715_01_A12-1121939868.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210722_01_D08-1122543703.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210722_01_D08-1122543703.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210722_01_B08-1122543704.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210722_01_B08-1122543704.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210722_01_H07-1122543705.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210722_01_H07-1122543705.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_210729_01_B09-1122543706.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_210729_01_B09-1122543706.multiome.atac.cellbarcode.txt
