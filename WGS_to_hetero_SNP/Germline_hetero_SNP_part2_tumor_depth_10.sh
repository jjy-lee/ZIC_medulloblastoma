#John Joo Yiul Lee
#directories from Sickkids hpf cluster and/or John Lee's laptop

sampleID=
vcfdir=/hpf/largeprojects/mdtaylor/jjylee/WGS/platypus_output_filter/${sampleID}
inputvcf=`ls ${vcfdir} | grep ${sampleID} | grep _platypus_germline_SNP_filtered_201001.txt | grep -v pos`
tumor_depth_threshold=10
outdir=/hpf/largeprojects/mdtaylor/jjylee/WGS/platypus_output_filter/${sampleID}
tempdir=${outdir}/temp

mkdir -p ${tempdir}
#create temporary vcf file with only PASS variants.
cat ${vcfdir}/${inputvcf} | grep PASS > ${tempdir}/temp_vcf

#further filter the germline SNPs using tumor WGS read depth.
while read line; do
 tumor_TOT_count=`echo ${line} | awk '{print $11}' | cut -f5 -d':'`
 #if tumor depth meets the threshold limit, pass the allele on
 if [ ${tumor_TOT_count} -ge ${tumor_depth_threshold} ]; then
  echo ${line} | tr ' ' '\t' >> ${outdir}/${sampleID}_platypus_germline_SNP_filtered_t_10_201030.txt
 fi
done < ${tempdir}/temp_vcf

rm -r ${tempdir}
