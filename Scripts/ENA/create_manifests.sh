while IFS=$'\t' read -r alias accession; do
  cat > manifest_${alias}.txt <<EOF
STUDY        PRJEB110793
SAMPLE       ${accession}
NAME         ${alias}
INSTRUMENT   NextSeq 2000
INSERT_SIZE  300
LIBRARY_SOURCE    TRANSCRIPTOMIC
LIBRARY_SELECTION unspecified
LIBRARY_STRATEGY  ssRNA-seq
FASTQ        ${alias}_R1_001_trimmoNohuman_1P.fq.gz
FASTQ        ${alias}_R1_001_trimmoNohuman_2P.fq.gz
EOF

done < link_sample_ID.txt
