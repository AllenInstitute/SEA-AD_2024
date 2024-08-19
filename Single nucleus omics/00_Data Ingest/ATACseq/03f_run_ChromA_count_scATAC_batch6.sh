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

ChromA_count	L8AT_211108_01_H01-1142581755.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_H01-1142581755.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_A02-1142581756.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_A02-1142581756.scATAC.cellbarcode.txt
ChromA_count	L8AT_211108_01_B02-1142581757.tsv.gz	grps13-6.merged.bed	L8AT_211108_01_B02-1142581757.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_F07-1144788667.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_F07-1144788667.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_E09-1145495811.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_E09-1145495811.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_F09-1145495812.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_F09-1145495812.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_G09-1145495813.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_G09-1145495813.scATAC.cellbarcode.txt
ChromA_count	L8AT_211117_01_H09-1145495814.tsv.gz	grps13-6.merged.bed	L8AT_211117_01_H09-1145495814.scATAC.cellbarcode.txt
ChromA_count	L8AT_211123_01_A12-1146036866.tsv.gz	grps13-6.merged.bed	L8AT_211123_01_A12-1146036866.scATAC.cellbarcode.txt
ChromA_count	L8AT_211123_01_B12-1146036867.tsv.gz	grps13-6.merged.bed	L8AT_211123_01_B12-1146036867.scATAC.cellbarcode.txt
ChromA_count	L8AT_211123_01_C12-1149386627.tsv.gz	grps13-6.merged.bed	L8AT_211123_01_C12-1149386627.scATAC.cellbarcode.txt
ChromA_count	L8AT_211123_01_D12-1149386628.tsv.gz	grps13-6.merged.bed	L8AT_211123_01_D12-1149386628.scATAC.cellbarcode.txt
