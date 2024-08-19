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

ChromA_count	L8AT_211117_01_G07-1143360947.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_G07-1143360947.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_H07-1143360948.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_H07-1143360948.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_A08-1143360949.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_A08-1143360949.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_B08-1143360950.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_B08-1143360950.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_C08-1143360951.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_C08-1143360951.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_D08-1143360952.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_D08-1143360952.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_E08-1143360953.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_E08-1143360953.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_E01-1149689563.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_E01-1149689563.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_F01-1149689564.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_F01-1149689564.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_G01-1149689565.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_G01-1149689565.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_H01-1149689566.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_H01-1149689566.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_A02-1149689567.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_A02-1149689567.scATAC.cellbarcode.txt
