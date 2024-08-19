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

ChromA_count	L8AT_210427_01_A05-1104418817.tsv.gz	grps13-6.merged.bed	L8AT_210427_01_A05-1104418817.scATAC.cellbarcode.txt
ChromA_count	L8AT_210427_01_B05-1104418819.tsv.gz	grps13-6.merged.bed	L8AT_210427_01_B05-1104418819.scATAC.cellbarcode.txt
ChromA_count	L8AT_220126_01_B09-1158586576.tsv.gz	grps13-6.merged.bed	L8AT_220126_01_B09-1158586576.scATAC.cellbarcode.txt
ChromA_count	L8AT_220126_01_C09-1158586577.tsv.gz	grps13-6.merged.bed	L8AT_220126_01_C09-1158586577.scATAC.cellbarcode.txt
ChromA_count	L8AT_220126_01_D09-1158586578.tsv.gz	grps13-6.merged.bed	L8AT_220126_01_D09-1158586578.scATAC.cellbarcode.txt
ChromA_count	L8AT_210830_01_G03-1131409694.tsv.gz	grps13-6.merged.bed	L8AT_210830_01_G03-1131409694.scATAC.cellbarcode.txt
