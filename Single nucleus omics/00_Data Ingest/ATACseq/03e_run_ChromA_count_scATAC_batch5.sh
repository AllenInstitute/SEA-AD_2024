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

ChromA_count	L8AT_211101_01_C09-1140885933.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_C09-1140885933.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_D09-1140885934.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_D09-1140885934.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_E09-1140885935.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_E09-1140885935.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_B10-1142534944.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_B10-1142534944.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_C10-1142534945.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_C10-1142534945.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_D10-1142534946.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_D10-1142534946.scATAC.cellbarcode.txt
ChromA_count	L8AT_211101_01_E10-1142534947.tsv.gz	grps13-6.merged.bed	L8AT_211101_01_E10-1142534947.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_C01-1142534948.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_C01-1142534948.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_D01-1142534949.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_D01-1142534949.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_E01-1142581752.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_E01-1142581752.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_F01-1142581753.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_F01-1142581753.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_G01-1142581754.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_G01-1142581754.scATAC.cellbarcode.txt
