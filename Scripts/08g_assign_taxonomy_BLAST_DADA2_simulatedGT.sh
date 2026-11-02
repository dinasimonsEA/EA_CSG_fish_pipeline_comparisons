# script to BLAST sequences (alignment-based)

#!/bin/bash
set -euo pipefail

test_data_name="Simulated"

base="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data"

processed_dir="${base}/Processed/${test_data_name}/"
query_dir="${processed_dir}/blast_queries_DADA2"

db_base="${base}/Databases/Simulated_GroundTruth"

for i in {1..10}
do

    sample="riaz_EA_${i}"

    query="${query_dir}/${sample}.fasta"

    db_name="${sample}_db"
    db_src="${db_base}/${db_name}"

    output="${processed_dir}/08_groudtruth_BLAST_DADA2/08_${sample}_blast.txt"

    echo "Running ${sample}"

    blastn \
        -query "${query}" \
        -task blastn \
        -db "${db_src}/${db_name}" \
        -evalue 10 \
        -perc_identity 95 \
        -outfmt "6 qseqid saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "${output}"

done

echo "DONE"