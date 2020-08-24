# Presence absence of gene models

# get useful information from Exonerate output
cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Exonerate/Output

ml R_packages/3.6.0
R

for EXON in $( ls EXONERATE*.output )
do
# parse Exonerate output: get only protein name, start end protein alignment, scaffold matched, start end alignment on genome, similarity
grep -B 1 'END OF\|^>' $EXON | grep -v '#\|\--' > $EXON.temp
# get table in proper shape
awk '!/>/ {print $1, $4, $5} $1 ~/>/ {print $0}' OFS="\t" $EXON.temp | paste -sd ' \n' | tr -s ' ' '\t' | sed 's/>//' > $EXON.table
Rscript --vanilla /proj/uppstore2018073/private/Valentina/lycPyr/Code/summaryGenePresence.R $EXON.table
done
