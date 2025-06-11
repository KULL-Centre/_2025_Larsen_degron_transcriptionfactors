
module load eautils

ls ../fastq_230817/*_R1_001.fastq.gz  |  parallel -P 16 ./call_zerotol_paired.sh {};
ls ../fastq_230905/*_R1_001.fastq.gz  |  parallel -P 16 ./call_zerotol_paired.sh {};

