## XRNs-DegradomeAnalysis

#### GuitarPlotFast.R  
Degradome signals density among mRNAs are calculated by Modified Guitar. Some functions of Guitar are modified from source code v2.10.0

#### FLEP-Seq_polyAEnds.R  
To find mapping reads of FLEP-Seq containing poly A and then assign 3' ends of alignemnts to non-reference A site.  
Bam files in dataset/Bam_FLEP_seq only contain three genes("AT4G27440","AT1G12090","AT1G15690").

#### FLEP-Seq_PolyAsite-GenomeToTranscriptCoordinates.R  
To correct position of published polyA sites and convert them from genome coordinates to transcript coordinates

#### Reference files:  
TAIR10_representative_gene_models.txt is download from TAIR10 database  
TAIR10_GFF3_genes_transposons.gtf is coverted from TAIR10_GFF3_genes_transposons.gff by gffread.
13059_2021_2543_MOESM2_ESM.xlsx if from Table S1 of Mo et al., Genome Biol 22, 322 (2021).

2022.12.15 written by Bo-Han Hou
