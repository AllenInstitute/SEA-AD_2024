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

ChromA_count	L8AT_211201_01_B02-1149689568.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_B02-1149689568.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_C02-1149689569.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_C02-1149689569.scATAC.cellbarcode.txt
ChromA_count	L8AT_211201_01_D02-1149689570.tsv.gz	grps13-6.merged.bed	L8AT_211201_01_D02-1149689570.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_E07-1149773119.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_E07-1149773119.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_F07-1149773121.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_F07-1149773121.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_G07-1149773124.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_G07-1149773124.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_H07-1149773127.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_H07-1149773127.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_E08-1149773128.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_E08-1149773128.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_A08-1149773129.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_A08-1149773129.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_B08-1149773130.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_B08-1149773130.scATAC.cellbarcode.txt
ChromA_count	L8AT_211215_01_C08-1149773131.tsv.gz	grps13-6.merged.bed	L8AT_211215_01_C08-1149773131.scATAC.cellbarcode.txt
ChromA_count	L8AT_210804_01_A06-1123939343.tsv.gz	grps13-6.merged.bed	L8AT_210804_01_A06-1123939343.scATAC.cellbarcode.txt
