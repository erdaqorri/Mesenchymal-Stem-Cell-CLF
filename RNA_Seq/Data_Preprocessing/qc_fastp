# Adapter and low base quality removal with FASTP 

input_dir="path/fastq"
output_dir="path/fastp_out"

for infile in "$input_dir"/*_R1.fastq.gz;
do
	base=$(basename "$infile" _R1.fastq.gz)
	read1="$input_dir/$base"_R1.fastq.gz
	read2="$input_dir/$base"_R2.fastq.gz

output_prefix="$output_dir/$base"
output_r1="$output_prefix"_R1_fastp.fastq.gz
output_r2="$output_prefix"_R2_fastp.fastq.gz
report_file="$output_prefix"_fastp.html

fastp --thread 16 --in1 "$read1" --in2 "$read2" --out1 "$output_r1" --out2 "$output_r2" adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 30 --length_required 35 -3 --html "$report_file" 

done
