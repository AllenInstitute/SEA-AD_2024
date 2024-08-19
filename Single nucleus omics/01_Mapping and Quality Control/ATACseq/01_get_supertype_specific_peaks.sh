#!/bin/bash

ALLEN_WG=/allen/programs/celltypes/workgroups/
ALLEN_PROD=/allen/programs/celltypes/production/
CHROMADIR=/allen/programs/celltypes/workgroups/hct/SEA-AD/environments

function ChromA() {
echo Running Chroma on file=$1 output=$2 and genome=$3
sbatch --ntasks-per-node=1 --cpus-per-task=32 --mem=500G -t 3-23:00:00 --qos=normal -p celltypes -J $2 -e $2.err --mail-type=FAIL <<END
#!/bin/bash

mkdir \$TMPDIR/{var,tmp}
export OWD=`pwd`

singularity exec --bind="\$TMPDIR/var/:/var" --bind="\$TMPDIR/tmp:/tmp" --bind "/sys/fs/cgroup/:/sys/fs/cgroup/" --bind "$ALLEN_WG:$ALLEN_WG" --bind "$ALLEN_PROD:$ALLEN_PROD" \
--contain $CHROMADIR/ChromA_biggermem.sif $CHROMADIR/ChromA_boilerplate.sh "$1" $2 $3 \$OWD

exit 0

END
}

ChromA "/allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter68Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter12Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter7Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter74Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter10Astro.tsv.gz
" Astro_ADNC0.bed hg38

ChromA "/allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter62Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter100Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter41Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter18Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter48Astro.tsv.gz
" Astro_ADNC1.bed hg38

ChromA "/allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter40Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter79Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter88Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter31Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter99Astro.tsv.gz
" Astro_ADNC2.bed hg38

ChromA "/allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter57Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter47Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter35Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter91Astro.tsv.gz /allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/subclasses_files/filter_files/Astro/fragments.tsv.gz_filter77Astro.tsv.gz
" Astro_ADNC3.bed hg38

