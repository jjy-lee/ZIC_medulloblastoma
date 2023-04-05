#!/bin/bash

##### saturation analysis for ChIP-Seq peaks
##### author: Zhen Y

##### usage: ./saturation_analysis.sh -f input_file.txt -i iteration -o output_dir
##### input_file.txt: text file contain path and peaks file name; each line represents each peak file (no header for peak file)
##### iteration: how many iteration you want, default is 10
##### output_dir: output directory. default is current dir; if output dir does not exist, will create one.

##### output results: itera_{1 .. iteration}.txt; files contain peaks number as sample cumulation at each iteration  



helpFunc()
{
	echo "
	Usage: ./saturation_analysis.sh 
	-f <files contain peaks bed files of all samples> 
	-i <iteration number, default: 10>
	-o <output directory> 
	-h <showing this message>
    "
}


iteration=10
DATE=`date +"%d-%m-%y-%T" | tr ":" "-"`
while getopts "f:i:o:h" opt; do
  case $opt in
    f ) files="$OPTARG"
    ;;
    i ) iteration="$OPTARG"
    ;;
    o ) output="$OPTARG"
    ;;
    h ) helpFunc ; exit 0
    ;;
    \? )
     echo "Invalid Option" 
     exit 1
    ;;
  esac
done


BEDTOOLS=`which bedtools || true`


if [[ -z "$BEDTOOLS" ]]
then
    echo -e "Error: Can not find your bedtools ... \n"
    exit 1
fi


if test -z "$files"; then
	echo "no input file provide, please prepare plain text file with peaks bed file from all samples; each sample per line ..."
	helpFunc
	exit 2
fi



if [[ ! -d ${output} ]]; then
  echo ""
  echo "output folder: ${output} doesn't exist"
  echo "
  Creating new output folder:
  ${output}

  "
  mkdir -p ${output}
fi

if test -z "$output"; then
	echo "no output dir provide, set to current directory ..."
	output="./"
fi



 
echo "

$DATE
Start saturation analysis ...

"

# step 1 create merge function and save result as a temp file
cmerg1() {
	FILE1=$1;
	cat temp3.txt "$FILE1" > x3.txt;
	mv x3.txt temp3.txt
	cat temp3.txt | sort -k1,1 -k2,2n | ${BEDTOOLS} merge -i stdin | wc -l 
}


cd ${output}
for i in $(seq 1 ${iteration});do
	echo -e "start iteration ${i} ... \n"
	files1=(`cat files | shuf`)  # shuffle order
	count=1
	wc -l ${files1[0]} >> itera_${i}.txt;
	cp ${files[0]} temp3.txt
	while [[ $count -le `expr ${#files1[@]} - 1` ]]; do
	echo ${count}
	cmerg1 ${files[${count}]} >> itera_${i}.txt;
	count=$((count + 1))
	done
	echo -e "iteration ${i} is done\n"
done





