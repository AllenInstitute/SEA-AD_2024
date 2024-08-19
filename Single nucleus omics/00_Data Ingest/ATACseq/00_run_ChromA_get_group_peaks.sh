#!/bin/bash

ALLEN_WG=/allen/programs/celltypes/workgroups/
ALLEN_PROD=/allen/programs/celltypes/production/
CHROMADIR=/allen/programs/celltypes/workgroups/hct/SEA-AD/environments

function ChromA() {
echo Running Chroma on file=$1 output=$2 and genome=$3
sbatch --ntasks-per-node=1 --cpus-per-task=32 --mem=500G -t 14-12:00:00 --qos=normal -p celltypes -J $2 -e $2.err --mail-type=FAIL <<END
#!/bin/bash

mkdir \$TMPDIR/{var,tmp}
export OWD=`pwd`

singularity exec --bind="\$TMPDIR/var/:/var" --bind="\$TMPDIR/tmp:/tmp" --bind "/sys/fs/cgroup/:/sys/fs/cgroup/" --bind "$ALLEN_WG:$ALLEN_WG" --bind "$ALLEN_PROD:$ALLEN_PROD" \
--contain $CHROMADIR/ChromA.sif $CHROMADIR/ChromA_boilerplate.sh "$1" $2 $3 \$OWD

exit 0

END
}

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1140030209/1140030209/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod3/ARC_Analysis_Run_1135448413/1135448413/outs/atac_fragments.tsv.gz
" grp_0.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1064293379/1064293379/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1122543704/1122543704/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1075757819/1075757819/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1122543706/1122543706/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1129176225/1129176225/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1129238386/1129238386/outs/atac_fragments.tsv.gz
" grp_1.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1064293374/1064293374/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1121747487/1121747487/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1074937000/1074937000/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1124629224/1124629224/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1074937004/1074937004/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1122543703/1122543703/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1074937006/1074937006/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1122543707/1122543707/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1131257166/1131257166/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1131257170/1131257170/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1131409694/1131409694/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1131257169/1131257169/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1136687632/1136687632/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1136687627/1136687627/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1138433801/1138433801/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1138433803/1138433803/outs/atac_fragments.tsv.gz
" grp_2.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1126220032/1126220032/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1124987484/1124987484/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1131409693/1131409693/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod3/ARC_Analysis_Run_1131607954/1131607954/outs/atac_fragments.tsv.gz
" grp_3.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1064293376/1064293376/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1121747486/1121747486/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1126220031/1126220031/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1124987483/1124987483/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1135042880/1135042880/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod3/ARC_Analysis_Run_1135704088/1135704088/outs/atac_fragments.tsv.gz
" grp_4.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1064293377/1064293377/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1122543705/1122543705/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1064293378/1064293378/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod3/ARC_Analysis_Run_1121939868/1121939868/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1074937005/1074937005/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1122543708/1122543708/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1126220030/1126220030/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1124987482/1124987482/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1135042881/1135042881/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod3/ARC_Analysis_Run_1135704089/1135704089/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1136687631/1136687631/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1136687625/1136687625/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1138433802/1138433802/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1138433804/1138433804/outs/atac_fragments.tsv.gz
" grp_5.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1123939343/1123939343/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1124629228/1124629228/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1123939344/1123939344/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1124629229/1124629229/outs/atac_fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1126220033/1126220033/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod0/ARC_Analysis_Run_1126220036/1126220036/outs/atac_fragments.tsv.gz
" grp_6.bed hg38

ChromA "/allen/programs/celltypes/production/mousecelltypes/prod0/ATAC_Analysis_Run_1140028669/1140028669/outs/fragments.tsv.gz /allen/programs/celltypes/production/mousecelltypes/prod3/ARC_Analysis_Run_1135448412/1135448412/outs/atac_fragments.tsv.gz
" grp_7.bed hg38


