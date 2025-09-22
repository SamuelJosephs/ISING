#! /usr/bin/bash 


output_dir="./output-dir"
output_file="concatenated_data.csv"

files_list=("${output_dir}"/*.csv)
files=$(ls ${output_dir}/)

header=$(head -n1 ${files_list[0]})

echo "$header" > $output_file



find ${output_dir}/ -type f -name "*.csv" | xargs -n 200  awk -F, '

	
	#FNR == 1 {next};

	$0 !~ /^[[:space:]]*$/ && ($1 + 0 == $1) {

		l = $0;
		l = gensub(/[[:space:]]+/,"","G");
		print l;
	}




' >> $output_file
