# L'analyse RNAseq - Tutoriel

Dans cette exemple d'analyse RNAseq, nous allons utiliser des données d'expression short-reads de cellule humaine HEK293 accessibles sous le BioProject [PRJNA753061](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA753061/).


## 1. Téléchargement des données d'intérêts : les lectures ARN

Pour pouvoir télécharger les lectures, nous allons utiliser [SRA Toolkit](https://github.com/ncbi/sra-tools). Mais en premier, pour faciliter "l'automatisation" du téléchargement des données, nous allons récupérer la liste des accessions des SRA. 

### 1.1. Récupération des noms SRA

- Aller à l'adresse du BioProject
- Cliquer sur *SRA* dans le menu *Related information* ou dans le tableau *Project Data* cliquer sur le *Number of Links* correspondant au nombre total de SRA
- Dans le menu déroulant *Send to* choisir *File* puis *Accession List* et valider avec *Create File*

Un fichier au nom *SraAccList.csv* a été créé.

### 1.2. Téléchargement des SRA

- Installer SRA toolkit
- Avec l'outil lancer la ligne de commande : 
`prefetch --option-file path/to/SraAccList.csv -O path/to/output/folder`, cela permet de récupérer les fichiers *.sra*

Warning, pour deux SRA (SRR15376509, SRR15376518) le message d'erreur suivant est renvoyé :
```
2024-08-19T10:10:23 prefetch.2.5.7 err: name incorrect while evaluating path within network system module - Scheme is 'https'
2024-08-19T10:10:23 prefetch.2.5.7 err: path not found while resolving tree within virtual file system module - 'SRR15376509' cannot be found.
```
Pour palier à cela, il faut se rendre directement sur la page du SRA via le navigateur de recherche du [SRA Run Browser](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&display=metadata). Puis récupérer dans la section *Data access* et la sous-section *SRA archive data* l'adresse http afin de faire un `wget`.

### 1.3. Obtention des FASTQ liés aux SRA

- Dans le dossier où il y a tous les dossiers SRA, rentrer la ligne de commande : `for dir in SRR*; do fastq-dump --split-3 --skip-technical --gzip "$dir"/"$dir".sra; done`, cela permet d'avoir les fichiers FASTQ paired-end \_1 et \_2 et unpaired compressés (voir la figure ci-dessous provenant du HowTo de [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) qui est un fastq-dump plus rapide mais sans la possibilité de compresser directement les fichiers)
- Supprimer les dossiers SRA pour alléger le volume données

Attention, après cette étape les *\_1* et *\_2* des noms des fichiers FASTQ ont été renommés en *_R1* et *_R2*.

## 2. Téléchargement des fichiers d'intérêts : les références

### 2.1. Le transcriptome de référence

Nous allons télécharger le transcriptome humain de référence depuis le genome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) (étant à ce jour la dernière version).

- Aller sur la page du génome assemblé
- Cliquer sur l'onglet *FTP*
- Une liste d'*Index of /genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14* sera proposée (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)
- Choisir *GCF_000001405.40_GRCh38.p14_rna.fna.gz*

### 2.2. Le génome de référence

Tout comme le transcriptome de référence, nous allons également prendre en possession le génome humain de référence.

- Aller sur la page du génome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Sélectionner l'onglet *FTP*
- Une liste d'Index of /genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14 sera proposée (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)
- Choisir *GCF_000001405.40_GRCh38.p14_genomic.fna.gz*

### 2.3. Le GTF

Nous allons récupérer les annotations génomiques du génome humain.

- Aller sur la page du génome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Sélectionner l'onglet *FTP*
- Une liste d'Index of /genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14 sera proposée (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)
- Choisir *GCF_000001405.40_GRCh38.p14_genomic.gtf.gz*

## 3. Préparation des lectures

### 3.1. Inférence du strandedness

Nous allons déterminer les spécifités des librairies avec Salmon (v1.10.1) et le transcriptome de référence pour pouvoir renseigner correctement les paramètres d'outils.
Puisque les données ont été séquencées au même moment par le même séquenceur avec la même préparation de librairie, nous n'avons pas besoin d'effectuer cette manipulation pour chaque SRA.
Voici les lignes de commande : 
```
WORKDIR=path/of/work/directory
gunzip "$WORKDIR"/GCF_000001405.40_GRCh38.p14_rna.fna.gz

# Index human transcriptome
salmon index -t "$WORKDIR"/GCF_000001405.40_GRCh38.p14_rna.fna -i "$WORKDIR"/GCF_000001405.40_GRCh38.p14_rna.index -k 31

# Run quant for having librarie type
salmon quant -i "$WORKDIR"/GCF_000001405.40_GRCh38.p14_rna.index -l A -1 "$WORKDIR"/SRR15376509_R1.fastq.gz -2 "$WORKDIR"/SRR15376509_R2.fastq.gz --validateMappings -o "$WORKDIR"/salmon_for_libtype 

gzip "$WORKDIR"/GCF_000001405.40_GRCh38.p14_rna.fna
```

Ici, l'option :

- `-l` permet de renseigner le type de librairie. En informant `A`, la détection se fait automatiquement !

Ce qui est important lors du lancement de salmon quant, c'est de bien regarder les sorties terminales (ou le log si vous les avez redirigées) notamment la ligne où il est indiquée :  [date] `[jointLog] [info] Automatically detected most likely library type as` [libtype].
Dans notre cas, voici la ligne obtenue : `[2024-08-22 11:45:41.747] [jointLog] [info] Automatically detected most likely library type as ISR`.

Ainsi, nos librairies sont `ISR`. De ce fait, nous pouvons dire que nos lectures sont **fr-firststrand**.

### 3.2. Contrôle qualité des lectures avant trimming

Avant de débuter les analyses des lectures, nous allons réaliser un contrôle qualité des lectures.
Pour ce faire, [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9) et [MultiQC](https://multiqc.info/) (v1.9) sont utilisés. FastQC permet de faire un rapport contrôle qualité sur chacun des fichiers FASTQ là où MultiQC permet de regrouper tous les rapports en un.
Voici les lignes de commande :
```
mkdir "$WORKDIR"/data_qc
mkdir "$WORKDIR"/data_qc/qc_before_trimming
fastqc "$WORKDIR"/*fastq* -t 16 -o "$WORKDIR"/data_qc/qc_before_trimming
multiqc -i "Human RNAseq" \
        -b "Quality control report before trimming step" \
        -n RNAseq_wo_trimming_report \
        -o "$WORKDIR"/data_qc/qc_before_trimming \
        "$WORKDIR"/data_qc/qc_before_trimming/*fastq* 
```
Ici, l'option :

- `-t` permet de déterminer le nombre de threads
- `-i` de caractériser le titre du raport MultiQC
- `-b` le sous-titre du raport MultiQC
- `-n` le nom du fichier contenant le raport MultiQC

### 3.3. Trimming des lectures

Dans le but d'avoir des lectures viables pour les analyses, nous effectuons un trimming. Cette étape est essentielle car elle permet l'élimination des adaptateurs et les portions de lectures de mauvaise qualité. Nous utilisons [fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file) (v0.20.1) comme outil qui permet de faire le trimming des adaptateurs, le filtre de qualité mais aussi pas mal de choses de façon très rapide. 
```
# Declare a array with all samples names
declare -a SAMPLES
SAMPLES=(
    'Name_sample_1'
    'Name_sample_2'
    'Name_sample_3'
)

THREADS=16

mkdir "$WORKDIR"/data_trimmed
mkdir "$WORKDIR"/data_trimmed/fastp_report

# For each samples names in the array
for i in ${!SAMPLES[@]}
do
    echo -e "*---------- TRIMMING STEP WITH FASTP : ${SAMPLES[$i]}"
    # Make fastp
    fastp -i "$WORKDIR"/"${SAMPLES[$i]}"_R1.fastq.gz \
          -o "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R1.trim.fastq.gz \
          -I "$WORKDIR"/"${SAMPLES[$i]}"_R2.fastq.gz \
          -O "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R2.trim.fastq.gz \
          --unpaired1 "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R1_unpaired.trim.fastq.gz \
          --unpaired2 "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R2_unpaired.trim.fastq.gz \
          --detect_adapter_for_pe \
          -j "$WORKDIR"/data_trimmed/fastp_report/"${SAMPLES[$i]}"_fastp.json \
          -h "$WORKDIR"/data_trimmed/fastp_report/"${SAMPLES[$i]}".html \
          -R "fastp report for ${SAMPLES[$i]}" \
          -w "$THREADS" \
          --cut_right
done
```

Ici, l'option :

- `-i` / `-I` permet de renseigner les fichiers FASTQ d'entrée
- `-o` / `-O` permet de renseigner les fichiers FASTQ de sortie. NB : si un `.gz` est ajouté à la fin du nom du fichier, il est alors automatiquement compressé
- `--unpaired1` / `--unpaired2` permet de renseigner les fichiers FASTQ de sortie unpaired.
- `--detect_adapter_for_pe` pour activer la détection automatique de la séquence d'adaptateurs
- `-j` pour renommer le rapport JSON 
- `-h` pour renommer le rapport HTML
- `-R` pour redéfinir le titre du rapport  
- `-w` caractérise le nombre de threads à utiliser
- `--cut_right` permet de déplacez une fenêtre glissante d'avant en arrière et si une fenêtre avec une qualité moyenne < threshold est rencontrée, la lecture est coupée à partir de la fenêtre. La taille de la fenêtre et le threshold de qualité sont conservés par défaut, soit 4 et 20. Cette option est similaire à `SLIDINGWINDOW:4:20` de [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

### 3.4. Contrôle qualité des lectures après trimming 

Pour voir si le trimming s'est bien passé, nous réalisons un contrôle qualité post-trimming. De la même façon que précedemment, FastQC et MultiQC sont utilisés.
```
WORKDIR=path/of/work/directory
mkdir "$WORKDIR"/data_qc/qc_after_trimming
fastqc "$WORKDIR"/data_trimmed/*fastq* -t 16 -o "$WORKDIR"/data_qc/qc_after_trimming
multiqc -i "Human RNAseq" \
        -b "Quality control report after trimming step" \
        -n RNAseq_trimming_report \
        -o "$WORKDIR"/data_qc/qc_after_trimming \
        "$WORKDIR"/data_qc/qc_after_trimming/*fastq* 
```

Egalement, pour rassembler tous les logs de fastp en un seul, nous lançons MultiQC sur les fichiers de sortie JSON : 
```
multiqc -i "Human RNAseq" \
        -b "Trimming report" \
        -n RNAseq_fastp_report \
        -o "$WORKDIR"/data_trimmed/fastp_report \
        "$WORKDIR"/data_trimmed/fastp_report/*.json
```
Attention, il faut bien que le pattern *fastp* apparaisse dans le nom des fichiers JSON de fastp pour que MultiQC marche.

## 4. Alignement des lectures viables

### 4.1. Indexation du génome de référence

Avant d'aligner, il est essentiel d'indexer le génome de référence. Nous allons utiliser dans ce cas, mais également pour l'alignement plus tard, [HiSat2](http://daehwankimlab.github.io/hisat2/manual/) (v2.2.1).

HiSat2 index le génome avec un algorithme qui permet de faire un compromis entre le temps d'exécution et l'allocation de la mémoire. Le gros avantage est qu'il recherche automatiquement les paramètres pour aboutir à cet objectif. 

HiSat2 peut prendre en entrée le génome ou le transcriptome de référence. Mais, nous utiliserons le génome de référence pour éviter des alignements discordants. En effet, avec le transcriptome de référence, et notamment avec la présence des isoformes, il peut y avoir des multi-alignements et donc créer des discordances dans la quantification (sous-estimer ou sur-estimer un transcript) et dans l'assemblage de novo.

```
gunzip "$WORKDIR"/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mkdir "$WORKDIR"/hisat2_index

# Build genome index
hisat2-build -p 16 -f "$WORKDIR"/GCF_000001405.40_GRCh38.p14_genomic.fna "$WORKDIR"/hisat2_index/GCF_000001405.40_GRCh38.p14_genomic

gzip "$WORKDIR"/GCF_000001405.40_GRCh38.p14_genomic.fna
```

Ici, l'option :

- `-p` indique le nombre de threads à utiliser 
- `-f` indique que la référence à indexer est un fichier FASTA

### 4.2. Alignement contre le génome indexé

Ici toujours, nous allons utiliser HiSat2 et le génome indexé préalablement pour pouvoir réaliser l'alignement des lectures trimmées contre ce dernier. Une fois le fichier d'alignement généré, nous trions les alignements dans l'ordre de position nucléotidique et compressons le fichier (de [SAM](https://fr.wikipedia.org/wiki/SAM_(format_de_fichier)) vers [BAM](https://fr.wikipedia.org/wiki/Binary_Alignment_Map)) avec [SAMtools](http://www.htslib.org/) (v1.18)
Le fait de trier nos alignements permet de poursuivre les analyses. Pour ce qui est de la compression, cela donne la possibilité de gagner en espace mémoire.

```
# Declare a array with all samples names
declare -a SAMPLES
SAMPLES=(
    'Name_sample_1'
    'Name_sample_2'
    'Name_sample_3'
)

THREADS=16

mkdir "$WORKDIR"/data_aligned
mkdir "$WORKDIR"/data_aligned/hisat2_report
for i in ${!SAMPLES[@]}
do
    echo -e "*---------- ALIGN WITH HISAT2 : ${SAMPLES[$i]}"
    hisat2 -x "$WORKDIR"/hisat2_index/GCF_000001405.40_GRCh38.p14_genomic \\
           -1 "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R1.trim.fastq.gz \\
           -2 "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R2.trim.fastq.gz \\
           -U "$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R1_unpaired.trim.fastq.gz,"$WORKDIR"/data_trimmed/"${SAMPLES[$i]}"_R2_unpaired.trim.fastq.gz \\
           -S "$WORKDIR"/data_aligned/"${SAMPLES[$i]}".sam \\
           --rna-strandness RF \\
           --downstream-transcriptome-assembly \\
           --summary-file "$WORKDIR"/data_aligned/hisat2_report/"${SAMPLES[$i]}".hisat2.txt \\
           -p "$THREADS"

    echo -e "*---------- SORT AND COMPRESS WITH SAMTOOLS : ${SAMPLES[$i]}"
    samtools sort "$WORKDIR"/data_aligned/"${SAMPLES[$i]}".sam -@ "$THREADS" -o "$WORKDIR"/data_aligned/"${SAMPLES[$i]}".sort.bam -O bam

    rm "$WORKDIR"/data_aligned/"${SAMPLES[$i]}".sam
done
```


Ici, l'option :

- `-x` permet de renseigner les noms des fichiers index  
- `-1` et `-2` indique les fichiers FASTQ R1 et R2 des lectures pairées
- `-U` indique les fichiers FASTQ R1 et R2 des lectures non pairées
- `-S` désigne le fichier SAM contenant les alignements à créer 
- `--rna-strandness` renseigne le type de librairie des lectures
- `--downstream-transcriptome-assembly` est une option utile pour l'assemblage car elle permet de conserver les lectures avec un long ancrage
- `-p` détermine le nombre de threads à utiliser

## 5. Assemblage des transcriptomes

### 5.1. Préparation du GTF/GFF

Lorsque nous allons assembler le transcriptome, nous allons avoir besoin du GTF pour guider l'assemblage. Mais avant de faire cela, il est important de vérifier notre fichier GTF et de le corriger si jamais. Pour cela, nous utilisons [AGAT](https://github.com/NBISweden/AGAT) (v1.4.0) un outil spécialisé dans les GTF/GFF.

```
agat_convert_sp_gxf2gxf.pl -g "$WORKDIR"/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -o "$WORKDIR"/GCF_000001405.40_GRCh38.p14_genomic.agat.gff
```
Par défaut, le fichier de sortie est au format GFF mais il est possible de changer les configs. Pour la suite, ce n'est pas génant que le fichier soit un GFF.

### 5.2. Assemblage guidé

Les alignements ainsi que le GFF sont prêts, nous allons donc réaliser les assemblages des transcriptomes. Pour le faire, nous utiliserons l'outil [StringTie](https://ccb.jhu.edu/software/stringtie/) (v2.2.3) qui est un assembleur rapide et efficace d'alignements RNAseq en transcripts potentiels.

```
# Declare a array with all samples names
declare -a SAMPLES
SAMPLES=(
    'Name_sample_1'
    'Name_sample_2'
    'Name_sample_3'
)

THREADS=16

mkdir -p "$WORKDIR"/data_assembly

echo -e "*---------- ASSEMBLY SAMPLES"
for i in ${!SAMPLES[@]}
do
    echo -e "*---------- ASSEMBLY WITH STRINGTIE : ${SAMPLES[$i]}"
    stringtie -G "$DIRREFS"/GCF_000001405.40_GRCh38.p14_genomic.agat.gff -o "$WORKDIR"/data_assembly/"${SAMPLES[$i]}".gtf -p "$THREADS" --rf "$WORKDIR"/data_aligned/"${SAMPLES[$i]}".sort.bam
done
```

L'option : 

- `-G` indique un fichier GFF (GFF3) contenant des annotations de gènes (comme les exons et les gènes) qui vont servir de guide pour l'assemblage. Cela permet à StringTie de mieux assembler les transcrits en se basant sur des annotations connues.
- `--rf` pour informer que les données sont fr-firststrand
- `-p` pour le nombre de threads


## 6. Quantification des transcrits et des gènes

### 6.1. Création du transcriptome global

Pour chacun des échantillons, nous avons assemblé les transcriptomes avec les lectures. Maintenant, nous allons construire le transcriptome global en le générant avec l'ensemble non redondant des transcrits observés dans les échantillons. Pour effectuer cela, nous allons encore une fois utilisé StringTie mais le mode `--merge`. Ce mode prends en entrée une liste de GTF et renvoie un seul et unique GTF.

```
stringtie --merge -G "$DIRREFS"/GCF_000001405.40_GRCh38.p14_genomic.agat.gff -o "$WORKDIR"/merged_transcriptomes.gtf "$WORKDIR"/data_assembly/*.gtf
```

### 6.2. Estimation de l'abondance des transcripts et des gènes

Toujours avec StringTie, nous allons estimer les abondances des transcripts et générer un tableau avec l'abondance des gènes. Ainsi, nous allons donner en entrée les alignements triées que nous avions produit avec la combinaison HiSat2 - SAMtools précédemment et le transcriptome que nous venons de créer afin d'avoir pour chacun des échantillons les TPM.

```
declare -a SAMPLES
SAMPLES=(
    'Name_sample_1'
    'Name_sample_2'
    'Name_sample_3'
)

THREADS=16

mkdir "$WORKDIR"/transcript_abundances

echo -e "*---------- TRANSCRIPT ABUNDANCES"
for i in ${!SAMPLES[@]}
do
    echo -e "*---------- TRANSCRIPT ABUNDANCES WITH STRINGTIE : ${SAMPLES[$i]}"
    mkdir "$WORKDIR"/transcript_abundances/sample_"${SAMPLES[$i]}"
    stringtie -e -A "$WORKDIR"/transcript_abundances/sample_"${SAMPLES[$i]}"/"${SAMPLES[$i]}".tsv -G "$WORKDIR"/merged_transcriptomes.gtf -o "$WORKDIR"/transcript_abundances/sample_"${SAMPLES[$i]}"/"${SAMPLES[$i]}".gtf -p "$THREADS" --rf "$WORKDIR"/data_aligned/"${SAMPLES[$i]}".sort.bam
done
```

Ici, l'option : 

- `-e` indique à l'outil d'estimer l'expression
- `-A` pour dire que nous voulons le tableau avec l'abondance des gènes

En sortie de l'option `-e`, nous obtenons un GTF, format que nous connaissons, avec tous les transcrits et leur TPM associé. Avec l'option `-A`, c'est un tableau tabulé avec l'abondance des gènes sous forme de TPM :

![image](https://gist.github.com/user-attachments/assets/21256cad-ee2d-488f-937a-39b4769f70ee)

Il faut noter que les TPM des gènes sont la somme des TPM des transcripts les constituants. Aussi, dans la colonne *Gene Name*, le nom du gène est reporté si connu. S'il n'est pas reporté, c'est que :
- le gène présente au moins un nouveau transcripts que va au delà du *start* et du *end* du gène connu
- la région est non annotée
- il y a absence de correspondance

Egalement, il peut y avoir une redondance de *Gene Name* car comme nous avons pris la version patchée du génome humain de référence le gène peut être présent sur un chromosome et sur un scaffold. Voici un exemple :

![image](https://gist.github.com/user-attachments/assets/360b5fcf-27bd-484c-80c0-242105e0fbc6)

## 7. Analyse de l'expression différentielle

### 7.1. Création des matrices de comptes

Nous avons précédemment estimé les transcripts. Nous allons par la suite utiliser le package R [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) pour réaliser l'expression différentielle.
Pour ce faire, nous devons d'abord transformer les sorties d'estimation d'abondance en entrée DESeq2 à l'aide du script python `prepDE.py` fourni par Stringtie.

```
python prepDE.py3 -i "$WORKDIR"/transcript_abundances -p sample -g "$WORKDIR"/transcript_abundances/prepDE_gene_count_matrix.csv -t "$WORKDIR"/transcript_abundances/prepDE_transcript_count_matrix.csv
```

L'option : 

- `-i` renseigne le dossier où il faut aller chercher les estimations
- `-p` indique le pattern des sous-dossier où se trouve les fichiers GTF avec les estimations
- `-g` et `-t` servent juste à donner des noms aux fichiers de sorties CSV

A la fin de cette ligne de commande, nous aurons deux fichiers CSV prêts à être utilisés par DESeq2 : l'un relatant les comptes pour les transcrits et l'un pour les gènes. Il est important de noter, comme nous l'indique Stringtie, pour le compte des transcrits : 

>prepDE.py derives hypothetical read counts for each transcript from the coverage values estimated by StringTie for each transcript, by using this simple formula: reads_per_transcript = coverage * transcript_len / read_len 

Les comptes des gènes, quant à eux, correspondent à la somme des `reads_per_transcript` les constituants. 


### 7.2. Préparation de la matrice d'informations

Les comptes sont prêts à être donnés à DESeq2, mais il manque encore une chose : la matrice des informations sur les échantillons. Celle-ci est indispensable dans notre analyse d'expression différencielle. Nous allons donc créer un fichier au format tabulé que nous appelerons *pheno_data.csv*.

Voici quelques lignes du fichier : 
```
sample	time	condition	selection	replica
sample_D2_mock_noSel_rep43	S0	mock	noSel	rep43
sample_D2_mock_noSel_rep41	S0	mock	noSel	rep41
sample_D2_mock_noSel_rep33	S0	mock	noSel	rep33
sample_D2_mock_noSel_rep32	S0	mock	noSel	rep32
```

### 7.3. L'analyse d'expression différencielle

Tous les fichiers sont prêts : 

- le fichier avec les comptes 
- le fichier avec les informations

Dorénavant, tout se fera sous R avec la librairie DESeq2. 

```
# Import libraries
library(DESeq2)
library("BiocParallel") # librarie to parallelize analysis
register(MulticoreParam(4)) # number of cpu

# Import gene count
countData <- read.csv("prepDE_gene_count_matrix.csv", row.names="gene_id")

# Remove last line with total values
countData <- countData[-nrow(countData),]

# Keep line where isn't only 0 values
countData <- countData[rowSums(countData != 0) > 0, ]

# Change the dataframe into a matrix
countData <- as.matrix(countData)

# Import informations about data
colData <- read.csv("pheno_data.csv", sep=",", row.names=1)

# Check
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

# Create dds object
ddsCond <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

# Filtering counts
smallestGroupSize <- 6 #choose the minimum redondance of a vector of informations
keep <- rowSums(counts(ddsCond) >= 10) >= smallestGroupSize
ddsCond <- ddsCond[keep,]

# Calculate the variance stabilizing transformation
vsdCond <- vst(ddsCond, blind=FALSE)

# PCA on vst results
pdf("pca_ddscond_vst.pdf")
plotPCA(vsdCond, intgroup = c("time"))
dev.off()

# Run the default analysis of DESeq2
ddsCond <- DESeq(ddsCond, parallel=TRUE)

# Save dds object
saveRDS(ddsCond, "DESeqDataSet.ddsTimeCond.rds")

# View names of results comparisons
resultsNames(ddsCond)
# [1] "Intercept"                  "condition_mock_vs_empty"   
# [3] "condition_shble1_vs_empty"  "condition_shble10_vs_empty"
# [5] "condition_shble2_vs_empty"  "condition_shble3_vs_empty" 
# [7] "condition_shble4_vs_empty"  "condition_shble5_vs_empty" 
# [9] "condition_shble6_vs_empty"

# Create results table of dds object
resCondMockEmpty <- results(ddsCond, name="condition_mock_vs_empty")
resCondSh1Empty <- results(ddsCond, name="condition_shble1_vs_empty")
resCondSh2Empty <- results(ddsCond, name="condition_shble2_vs_empty")
resCondSh3Empty <- results(ddsCond, name="condition_shble3_vs_empty")
resCondSh4Empty <- results(ddsCond, name="condition_shble4_vs_empty")
resCondSh5Empty <- results(ddsCond, name="condition_shble5_vs_empty")
resCondSh6Empty <- results(ddsCond, name="condition_shble6_vs_empty")
resCondSh10Empty <- results(ddsCond, name="condition_shble10_vs_empty")
```

Une fois les résultats créés, vous pouvez les enregistrer sous format .csv sur votre ordinateur.