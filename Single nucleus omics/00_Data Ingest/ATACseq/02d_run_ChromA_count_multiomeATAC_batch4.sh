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

ChromA_count	L8XR_211007_02_F03-1135448413.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211007_02_F03-1135448413.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_211007_02_G03-1135704088.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211007_02_G03-1135704088.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_211007_02_C04-1135704089.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211007_02_C04-1135704089.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_211014_02_D05-1136687625.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211014_02_D05-1136687625.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_211014_02_H05-1136687627.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211014_02_H05-1136687627.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_211021_02_E08-1138433803.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211021_02_E08-1138433803.multiome.atac.cellbarcode.txt
ChromA_count	L8XR_211021_02_G06-1138433804.multiome.atac.tsv.gz	grps13-6.merged.bed	L8XR_211021_02_G06-1138433804.multiome.atac.cellbarcode.txt
