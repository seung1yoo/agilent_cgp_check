bedtools coverage -a MANE.GRCh38.v1.2.refseq_genomic.exon.gtf -b A3416642_Covered.bed > MANE.GRCh38.v1.2.refseq_genomic.exon.CoverageWith_A3416642.bed
cat MANE.GRCh38.v1.2.refseq_genomic.exon.CoverageWith_A3416642.bed | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$11"\t"$12"\t"$13"\t"$9}' > MANE.GRCh38.v1.2.refseq_genomic.exon.CoverageWith_A3416642.select_column.tsv')))
}
