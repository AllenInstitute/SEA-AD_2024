cat grp_0.bed_allpeaks.bed grp_1.bed_allpeaks.bed grp_3.bed_allpeaks.bed grp_4.bed_allpeaks.bed grp_5.bed_allpeaks.bed grp_6.bed_allpeaks.bed > grps13-6.concat.bed

bedSort grps13-6.concat.bed grps13-6.concat.sorted.bed

bedtools merge -i grps13-6.concat.sorted.bed > grps13-6.merged.bed
