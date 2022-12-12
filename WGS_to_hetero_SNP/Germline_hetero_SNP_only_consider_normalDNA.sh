#John Joo Yiul Lee
#directories from Sickkids hpf cluster and/or John Lee's laptop

sampleID=
vcfdir=/hpf/largeprojects/mdtaylor/hsuzuki/Platypus/${sampleID}
inputvcf=`ls ${vcfdir} | grep ${sampleID} | grep _T_VariantCalls_merge.vcf.gz | grep -v tbi`
blood_thresh_lower=0.3
blood_thresh_upper=0.7
blood_depth_threshold=7
tumor_depth_threshold=7
outdir=/hpf/largeprojects/mdtaylor/jjylee/WGS/platypus_output_filter/blood_tumor_separate_SNP_call/${sampleID}
tempdir=${outdir}/temp

mkdir -p ${tempdir}
#create temporary vcf file with only PASS variants.
zcat ${vcfdir}/${inputvcf} | grep PASS > ${tempdir}/temp_vcf

#REF and ALT - if any are too long, don't include. Most stringent would be 1, single nucleotides. 4 would be more loose.
#7 to 10 for total number of reads from germline. Use 7 for this analysis.
#For identifying SNPs that may have been lost in the tumors due to deletion, tumor read depth might be a little bit lower.
#Use tumor depth threshold of 7 instead of 10.

#write an awk function to compare between numbers.
#if n1 is bigger, 1
#if n1 is smaller, 2

numCompare() {
   awk -v n1="$1" -v n2="$2" 'BEGIN {if(n1 >= n2){print 1} else {print 2}}'
}

while read line; do
 blood_TOT_count=`echo ${line} | awk '{print $10}' | cut -f5 -d':'`
 blood_VAR_count=`echo ${line} | awk '{print $10}' | cut -f6 -d':'`
 tumor_TOT_count=`echo ${line} | awk '{print $11}' | cut -f5 -d':'`
 tumor_VAR_count=`echo ${line} | awk '{print $11}' | cut -f6 -d':'`
 alt_comma=`echo ${line} | awk '{print $5}' | grep \, | wc -l`
 #alt_comma should be 0. Do not use multiallelic variants
 ref_length=`echo ${line} | awk '{print $4}' | wc -c`
 alt_length=`echo ${line} | awk '{print $5}' | wc -c`
 blood_ratio=`awk -v TOT_B="${blood_TOT_count}" -v VAR_B="${blood_VAR_count}" 'BEGIN {print VAR_B/TOT_B}'`
 #if blood ratio meets the threshold limit, pass the allele on
 if [ `numCompare ${blood_ratio} ${blood_thresh_lower}` -eq 1 ] && [ `numCompare ${blood_ratio} ${blood_thresh_upper}` -eq 2 ] && [ ${blood_TOT_count} -ge ${blood_depth_threshold} ] && [ ${alt_comma} -eq 0 ] && [ ${ref_length} -le 3 ] && [ ${alt_length} -le 3 ] && [ ${tumor_TOT_count} -ge ${tumor_depth_threshold} ] ; then
  echo ${line} | tr ' ' '\t' >> ${outdir}/${sampleID}_platypus_germline_SNP_filtered_220416.txt
 fi
done < ${tempdir}/temp_vcf
rm -r ${tempdir}
