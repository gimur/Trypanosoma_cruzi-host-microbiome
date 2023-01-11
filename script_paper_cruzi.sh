#Project: Parasite-Microbiota_Host interaction in a murine model of Trypanosoma cruzi infection (Tulahuen strain).
#Data obtained by shotgun metagenomics (4Gb per sample).
# BALBc and BL6 mice
#For each model there are shots at 0-3-7-10-13-16 days 
#Days 0 - 10 - 16 have a NI
#NI = Not Infected
#DPI: Days Post Infection
#Raw reads

#Several folders are created there:
#In manual analysis there are each of the analysis subfolders: #In manual analysis there are each of the analysis subfolders.
#In tools there are the tools that have been installed outside of the controlled environment.


#Preprocessing

#Load working environment in the cluster
module load conda
source activate /datacnmat01/biologia/gimur/Anaconda/Anaconda3/envs/gimur-tools

#A confirmatory quality control is performed by means of fastqc.
#Fastqc installation was performed with the following command:

conda install -c bioconda/label/broken fastqc

fastqc *.gz

# Use of trimmomatic
# The tool was installed with the command:

conda install -c bioconda/label/broken trimmomatic

# Reads are trimmed to confirm the removal of adapters,
# to maintain sequences of at least 150 bp
# to perform filtering by quality (Phred Score)
# Repeat the following command for each sampl

trimmomatic PE -threads 16 _1.fq.gz _2.fq.gz _1_paired_trim.fq.gz _1_unpaired_trim.fq.gz _2_paired_trim.fq.gz _2_unpaired_trim.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads MINLEN:150 AVGQUAL:20 TRAILING:20

#Mapping with host genome with Bowtie2. Previously written script was used and adjusted for the present analysis.
#Manual mapping is possible. For this case the adapted script bowtie2_mice is used:

----------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
# Bowtie2 wrapper script to filter out contaminations from metagenomic samples.
# Adapted by: Sergio Castañeda
# Last updated on: 2021-09-11

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Remove contaminant sequences with Bowtie2.
Usage: metabiome bowtie2 [options] -i <in_dir> -o <out_dir> -opts bowtie2_options

Required:
  -i in_dir         Input directory containing FASTQ files.
  -o out_dir        Directory in which results will be saved. This directory will
                    be created if it doesn't exist.

Options:
  -ho host          Host reference genome in FASTA format.
  -ph PhiX          PhiX-174 phage reference genome in FASTA format.
  -m  mus           Mus musculus reference genome in FASTA format.
  -t  NUM           Number of threads to use. (default=1)
  -opts OPTIONS     Bowtie2's options.
  -h, --help        Show this help.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -ho )       host=$(readlink -f "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -ph )       PhiX=$(readlink -f "$2"); shift 2 ;;
        -m  )       Mus=$(readlink -f "$2"); shift 2 ;;
        -opts )     shift; bowtie2_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Activate Conda environment
#activate_env metabiome-preprocessing

# Download Mice and PhiX reference genomes
# Highlight: Next code downloads PhiX and Mus musculus genome and checks if the
# downloads were successfull
cd "$out_dir"
if [[ ! -e "$Mus" ]]; then
    URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
    for i in {1..10}; do
        echo "Downloading Mus musculus reference genome"
        wget "$URL" -P "$out_dir"
        if [[ -s GCF_000001635.27_GRCm39_genomic.fna.gz ]]; then
            gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
            echo "Mus musculus genome was downloaded"
            Mus="$(readlink -f GCF_000001635.27_GRCm39_genomic.fna)"
            break
        else
            echo "Try downloading again the Mus musculus reference genome"
        fi
    done
fi

if [[ ! -e "$PhiX" ]]; then
    for i in {1..10}; do
        echo "Downloading PhiX reference genome"
        esearch -db nucleotide -query "NC_001422.1" | efetch -format fasta > PhiX_NC_001422.1.fasta
        if [[ -s PhiX_NC_001422.1.fasta ]]; then
            echo "PhiX genome was downloaded"
            PhiX=$(readlink -f PhiX_NC_001422.1.fasta)
            break
        else
            echo "Try downloading again the PhiX reference genome"
        fi
    done
fi

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Phix Genome: ${PhiX:?'PhiX genome not set'}"
echo "Mus musculus Genome: ${Mus:?'Mus musculus genome not set'}"
echo "Bowtie2 called with options: $bowtie2_opts"

# Generate the mixed fasta and its index
# Small index if the genome is small
if [[ ! -f "$host" ]] && [[ ! -f Mix.rev.2.bt2 ]]; then
    # Pool together the Mus and host genome
    cat "$PhiX" "$Mus" > Mixed.fasta

    # Generate a small index
    echo "Small index needs to be generated"
    bowtie2-build Mixed.fasta Mix --threads "$threads"

# Large index if the genome is large due to the presence of the host genome
elif [[ -f "$host" ]] && [[ ! -f Mix.rev.2.bt2l ]]; then
    # Pool together the host, phage and Mus musculus genome
    cat "$host" "$PhiX" "$Mus" > Mixed.fasta

    # Generate a large index
    echo "Large index needs to be generated"
    bowtie2-build --large-index Mixed.fasta Mix --threads "$threads"
fi

# Alignment
for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        bowtie2 -x Mix -1 "$forward_file" -2 $(forward_to_reverse "$forward_file") -p "$threads" \
            -q --un-conc-gz "$out_dir"/$(echo "$core_name" | sed 's/_trim/_bt2/').fq.gz \
            2> "$out_dir"/$(echo "$core_name" | sed 's/_trim/_bt2/')_summary.txt \
            $bowtie2_opts > /dev/null # Bowtie2 output to terminal is excesive and we do not need it in this case

    # Unpaired reads
    elif [[ ! "$file" ==  *_@(R1_|1.|R2_|2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        bowtie2 -x Mix -U "$unpaired_file" -q -p "$threads" \
            --un-gz "$out_dir"/$(basename -- "$unpaired_file" | sed 's/_trim/_bt2/') \
            2> "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_trim/_bt2/')_summary.txt \
            $bowtie2_opts > /dev/null

    # Files that do not match the required extension
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# Rename bowtie2 output files
rename "s/.fq.1.gz/_1.fq.gz/" *.fq.1.gz 2> /dev/null
rename "s/.fq.2.gz/_2.fq.gz/" *.fq.2.gz 2> /dev/null

# Set the correct file extension (.gz) for unpaired output files
rename -f "s/fq/fq.gz/" *.fq 2> /dev/null
rename -f "s/fastq/fastq.gz/"  *.fastq 2> /dev/null

echo "Done."
echo "You can now use clean reads to:"
echo "- Assemble genomes using metabiome metaspades or megahit"
echo "- Perform taxonomic binning with metabiome kraken2, kaiju or metaphlan3"
echo "- Perform functional profiling using metabiome humann"
echo "- Extract 16S sequences with metabiome bbduk"

------------------------------------------------------------------------------------------------------------------------------------------------------------------


#From Clean reads in Clean_reads folder: 

#Centrifuge tool is explored, which showed better resolution at species level.

conda create --yes -n centrifuge centrifuge
conda activate centrifuge

#The centrifuge database was generated.

curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz

# alternatively we can use wget

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz

# once the download is finished, we need to extract the archive content
# It will extract a few files from the archive and may take a moment to finish.

tar -xvzf p_compressed+h+v.tar.gz

#The Tool is run. In the -x argument you put the path where the previously generated base was left. The previous step generates p_compressed+h+v with numbers.
#Simply put it without any number. -1 and -2 correspond to the PE reads. Change the names of report.txt and results.txt.
#This is done for each sample

centrifuge -x p_compressed+h+v -1 example.1.fq -2 example.2.fq -U single.fq --report-file report.txt -S results.txt


#In order to build graphs with Pavian in RStudio you need a Kraken type format. Therefore you must transform the file results.txt
#to kreport.txt
#This is done for each sample

centrifuge-kreport -x /p_compressed+h+v RT1_results.txt > RT1_kreport.txt

#With this kreport you can perform a statistics and assignment visualization with Pavian in RStudio.
#This is done for each sample
#Open RStudio and run it:

if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

pavian::runApp(port=5000)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")


#This generates an interactive window where you can view various analyses.
#From there it is worth downloading the CSVs with the assignments for Bacteria-Viruses-Eukaryotes to be considered.

#An exploratory figure is generated with the 10 most abundant taxa.
----------------------------------------------------------------------------------

# library
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(viridis)
library(hrbrthemes)
library(forcats)
library(cowplot)



cruzi_pivot <- file %>% pivot_longer(

  cols = 2:11,

  names_to = "Species",

  values_to = "Abundance",

  values_drop_na = TRUE)

# Stacked + percent
ggplot(cruzi_pivot, aes(fill=Species, y=Abundance, x=Name)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_brewer(palette = "Spectral")

--------------------------------------------------------------------------------


#Parallel to this, a functional analysis is performed on the reads with Humann3 (https://github.com/biobakery/biobakery/wiki/humann3) (https://github.com/biobakery/humann).
#The process indicated for its installation and use is followed https://huttenhower.sph.harvard.edu/humann/
#Important to work the installation and use inside conda biobakery3 so that the dependencies and others work correctly.


conda create --name biobakery3 python=3.7
conda activate biobakery3
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
conda install humann -c biobakery
humann test

#Base de Datos Pangenoma
humann_databases --download chocophlan full  /databases_humann/ --update-config yes


#Base de datos proteínas
humann_databases --download uniref uniref90_diamond /databases_humann/ --update-config yes

#Base de datos anotaciones

humann_databases --download utility_mapping full /databases_humann/ --update-config yes


bowtie2 -x /databases_humann/mpa_v30_CHOCOPhlAn_201901 -U merge_sample.fq.gz


#To perform the functional analysis on the readings directly (recommended) first do a cat of the paired end and then run Humann

#cat to do the merge of the F and R paired readings in

cat BL60NI_R1.fq.gz BL60NI_R2.fq.gz > merge_BL60NI.fq.gz

#Humann for the reads. Nucleotide and protein database parameters are indicated.
humann --threads 18 -i /merge_sample.fq.gz -o /sample --nucleotide-database /databases_humann/chocophlan --protein-database /databases_humann/uniref


#With the following commands you can visualize and inspect the output which are 3 files plus a temp folder.

column -t -s $'\t' merge_sample_genefamilies.tsv | less -S


# To facilitate comparisons between samples with different sequencing depths, it is important to normalize the RPK values 
#predetermined HUMAnN values before performing statistical analyses. 
#Sum-normalization to relative abundance or "copies per million" (CPM) units are common approaches for this; these methods are implemented in humann_renorm_tablescript 
#(with special considerations for HUMAnN output formats, e.g., feature layering). Let's normalize our gene abundance to CPM units.

humann_renorm_table --input merge_sample_genefamilies.tsv --output merge_sample_genefamilies-cpm.tsv --units cpm --update-snames

#With the following commands you can visualize and inspect the output.

column -t -s $'\t' merge_sample_genefamilies-cpm.tsv | less -S

# The default "units" of HUMAnN microbial function are the UniRef gene families (which we use internally to calculate reaction abundances and, 
# from there, the pathway abundances). However, from the abundance of gene families, it is possible to reconstruct the abundance of other functional categories 
# in a microbiome using the humann_regroup_tablescript. Inspect the regrouping options using the humann_regroup_table -hcommand.

# We regrouped our CPM-normalized gene family abundance values to the MetaCyc reaction (RXN) abundances, which are included with the default HUMAnN installation.

humann_regroup_table --input merge_sample_genefamilies-cpm.tsv  --output merge_sample_genefamilies_uniref90rxn-cpm.tsv --groups uniref90_rxn


#It is generally useful to attach some human-readable descriptions of these IDs to facilitate biological interpretation.

humann_rename_table --input merge_sample_genefamilies_uniref90rxn-cpm.tsv --output merge_sample_genefamilies_uniref90rxn-cpm-named.tsv --names metacyc-rxn

column -t -s $'\t' merge_sample_genefamilies_uniref90rxn-cpm-named.tsv | less -S


# Take all the final files to the same folder to unite all the output -named.tsv. Run the following command

humann_join_tables -i . -o sample_Reads_Final_Humann.tsv --file_name uniref90rxn-cpm

#To plot correctly using Humann, a stratified table must be obtained in Humann3 with the following command:

humann_split_stratified_table --input BL6_Reads_Final_Humann.tsv --output BL6_Reads_Final_Humann_estratificada.tsv

# Download this layered file to work in Maaslin2 in R Studio (https://github.com/biobakery/biobakery/wiki/maaslin2#1-installing-r)
# Important to review the downloaded table, organize it according to the structure required by the tool. Do not receive such long column names
# It is necessary to adjust manually 
# In RStudio you work as follows:

--------------------------------------------------------------------------------------------------------------------------------

library(dplyr)

# Maaslin2 execution

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

library(Maaslin2)
library(readr)

# As input data you can put the taxonomic assignment results 
# as well as the results of the analysis by Humann3
# Here the consolidated .tsv obtained by Humann3 in command line is used as input_data.

df_input_data = read.table(file = "/sample_Final_Humann.tsv", header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]
df_input_data=as.data.frame(df_input_data)


# Input metadata correspondiente a las muestras evaluadas

df_input_metadata = read.table(file = "/metadata.tsv", header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]


# Ejecución del análisis

fit_func = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "demo_functional4", 
  fixed_effects = c("Condition"),
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  normalization = "NONE",
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

--------------------------------------------------------------------------------------------------------------------------------

# Once we know which roads are statistically significant, we search directly with the corresponding term in the stratified table. 
# To this table we must add the metadata by rows of what we want to consider 
# for each sample according to the figure to be constructed

grep 'IMP-DEHYDROG-RXN' sample_Reads_Final_Humann_stratified.tsv | less -S

#Figure is generated

humann_barplot --input sample_Reads_Final_Humann_stratified.tsv --focal-metadata Time --last-metadata Condition --output BALBc_RXN-14117_2.png --focal-feature 'RXN-14117' --sort sum metadata --scaling logstack

#Binning


bowtie2-build contigs.fasta final.contigs

bowtie2 -x final.contigs -1 *R1.fq.gz -2 *R2.fq.gz | samtools view -bS -o /*_to_sort.bam

samtools sort *_to_sort.bam -o *.bam

samtools index *.bam

conda install metabat2

runMetaBat.sh -m 1500 -t 18 contigs.fasta *.bam

conda install maxbin2

run_MaxBin.pl -thread 16 -contig contigs.fasta -reads *_R1.fq.gz -reads2 *_R2.fq.gz -out name

git clone https://github.com/BinPro/CONCOCT.git
cd CONCOCT
pip install -r requirements.txt
python setup.py install

cut_up_fasta.py contigs.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

concoct_coverage_table.py contigs_10K.bed *.bam > coverage_table.tsv

# binning

concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_*/

merge_cutup_clustering.py concoct_*/clustering_gt1000.csv > concoct_*/clustering_merged.csv

mkdir concoct_*/fasta_bins

extract_fasta_bins.py contigs.fasta concoct_*/clustering_merged.csv --output_path concoct_*/fasta_bins

#CheckM

conda create -n checkm python=3.9
conda activate checkm
conda install numpy matplotlib pysam
conda install hmmer prodigal pplacer
pip3 install checkm-genome

wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvzf

checkm data setRoot .


checkm lineage_wf -t 18 -x fasta bins ./

checkm qa lineage.ms . -o 2 -f qa


# DAS Tool

module load usearch

# Config input  DAS Tool by tool

src/Fasta_to_Scaffolds2Bin.sh -i /tool/BALBc16NI -e fasta > tool.scaffolds2bin.tsv

# CONCOCT

perl -pe "s/,/\tconcoct./g;" concoct_*/clustering_gt1000.csv > concoct.scaffolds2bin.tsv


# DAS Tool

DAS_Tool -i /maxbin.scaffolds2bin.tsv,metabat.scaffolds2bin.tsv,concoct.scaffolds2bin.tsv -l maxbin,metabat,concoct -c contigs.fasta -t 16  --write_bins 1 -o *_Result

#GTDB Tk

conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk
conda activate gtdbtk

download-db.sh

#Workflow gtdbtk classify wf

#identify, aligny y classify.

gtdbtk classify_wf --e fa --genome_dir ensamblaje --out_dir out/gtdbtk --pplacer_cpus 1 --scratch_dir .

#Anotation

conda install -c conda-forge -c bioconda -c defaults prokka
conda install -c bioconda perl-bioperl
prokka --version

for k in *.fna; do prokka $k --outdir "$k".prokka.output --prefix PROKKA_$k; echo $k; done

prokka *.fasta --outdir *


#Pangenoma Roary
ml load miniconda3/4.8.3
conda create --name NOMBRE
conda env remove --name NOMBRE
mamba install roary
mamba install -c conda-forge -c bioconda -c defaults panaroo
mamba create -n bactopia -c conda-forge -c bioconda bactopia

roary -f roary_* -e -n -v -r -qc *.gff
FastTree -nt -gtr core_gene_alignment.aln > my_tree.newick


#Pangenoma Panaroo

panaroo -i *.gff -o panaroo_out --clean-mode strict -a core --aligner mafft --core_threshold 0.99 -t 32
FastTree -nt -gtr core_gene_alignment.aln > my_tree.newick







