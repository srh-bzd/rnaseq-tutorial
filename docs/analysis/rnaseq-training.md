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

![image](https://gist.github.com/user-attachments/assets/c48eff7c-6e0f-490d-8308-e59fcc3a7463)

Attention, après cette étape les *\_1* et *\_2* des noms des fichiers FASTQ ont été renommés en *_R1* et *_R2*.

## 2. Téléchargement des fichiers d'intérêts : les références

### 2.1. Le transcriptome de référence

Le transcriptome de référence humain est une représentation complète de tous les transcrits (ARN) exprimés dans les cellules humaines à un moment donné. Contrairement au génome, qui est fixe et identique dans toutes les cellules d'un individu, le transcriptome varie en fonction du type cellulaire, du stade de développement, des conditions environnementales, et des stimuli externes. Il constitue donc un instantané dynamique de l'expression génique dans une cellule ou un tissu.

Nous allons télécharger le transcriptome humain de référence depuis le genome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) (étant à ce jour la dernière version).

#### 2.1.1. Téléchargement du transcriptome de référence

- Aller sur la page du génome assemblé
- Cliquer sur l'onglet *FTP*
- Une liste d'*Index of /genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14* sera proposée (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)
- Choisir *GCF_000001405.40_GRCh38.p14_rna.fna.gz*

#### 2.1.2. Informations sur le transcriptome de référence

**Le nom du fichier**

Nous allons décortiquer le nom du fichier du transcriptome en premier : 
- *GRCh38* (Genome Reference Consortium Human Build 38) correspond à la version du génome humain assemblé
- *p14* (patch release 14) correspond à la version 14 du patch, soit à la 14ème correction
- *GCF_000001405.40* est l'accession RefSeq (l'équivalent en accession GenBank est GCA_000001405.29)
- *_rna.fna.gz* correspond au format FASTA comprennant les ARN produits à partir du génome

**Accessions number et status code**

Pour ce qui est de l'intérieur du fichier, rien de changant par rapport à un fichier FASTA classique.
En revanche, on peut noter 4 types de numéro d'accessions dans les ids (après *>* dans le fichier) : *NM_*, *XM_*, *NR_* et *XR_*.
Si l'on regarde la page RefSeq de l'annotation d'[Homo sapiens](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Homo_sapiens/GCF_000001405.40-RS_2023_10/), notamment le tableau *Feature counts* dans la section *Gene and feature statistics*, on peut retrouver les patterns d'ids dans *mRNAs* (pour *NM* et *XM*) et dans *non-coding RNAs* (pour *NR* et *XR*). 

Les ids débutant par un *X* sont des *Model RefSeq*. Ce sont des séquences produites de façon automatique par un pipeline.
>  RNA and protein products that are generated by the eukaryotic genome annotation pipeline. These records use accession prefixes XM_, XR_, and XP_. (https://www.ncbi.nlm.nih.gov/refseq/about/)

Les ids débutant par un *N* sont des *Known RefSeq*. Ce sont des séquences manuellement revues par des membres du NCBI ou des collaborateurs.
>  RNA and protein products that are mainly derived from GenBank cDNA and EST data and are supported by the RefSeq eukaryotic curation group. These records use accession prefixes NM_, NR_, and NP_. (https://www.ncbi.nlm.nih.gov/refseq/about/)

(Pour voir plus sur les préfixes d'accessions et les codes RefSeq : https://gist.github.com/srh-bzd/94557751eef39e38a4235d00c54b93ea)

**Quantification des ARNs présents**

Comme nous l'avons vu précédement, nous avons 4 types d'accessions number correspondant à des ARNm et à des ARNnc.
Voici le nombre de chacun d'entre eux dans le transcriptome de référence téléchargé. Les informations peuvent être retrouvées dans la table *Feature counts* de la page RefSeq de l'annotation d'[Homo sapiens](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Homo_sapiens/GCF_000001405.40-RS_2023_10/), notamment aux lignes *mRNAs*, *non-coding RNAs* et *pseudo transcripts*.

<table>
  <tr>
       <th colSpan="2">mRNAs</th>  
       <th colSpan="2">ncRNAs</th>
       <th colSpan="2">pseudo transcripts</th>
  </tr>
  <tr>
       <th>NM</th>
       <th>XM</th>
       <th>NR</th>
       <th>XR</th>
       <th>NR</th>
       <th>XR</th>
   </tr>
  <tr>
       <th>67116</th>
       <th>69065</th>
       <th>21487</th>
       <th>25697</th>
       <th>1593</th>
       <th>163</th>
   </tr>
  <tr>
       <th colSpan="2">136181</th>  
       <th colSpan="2">47184</th>
       <th colSpan="2">1756</th>
  </tr>
  <tr>
       <th colSpan="6">185121</th>  
  </tr>
</table>

**Les isoformes**

Dans le transcriptome de référence humain, les isoformes des ARN sont présents. 
Prenons le cas du gène *MELK maternal embryonic leucine zipper kinase*. 
Lorsque l'on effectue une recherche de ce gène dans le fichier FASTA du transcriptome, on retrouve 53 séquences de transcripts de mRNA prédits ou non et d'un ncRNA review 
```
NM_001256685.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 2, mRNA
NM_001256687.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 3, mRNA
NM_001256688.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 4, mRNA
NM_001256689.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 5, mRNA
NM_001256690.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 6, mRNA
NM_001256691.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 7, mRNA
NM_001256692.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 8, mRNA
NM_001256693.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 9, mRNA
NM_014791.4	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 1, mRNA
NR_046337.2	Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant 10, non-coding RNA
XM_011518076.3	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X2, mRNA
XM_011518077.2	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X3, mRNA
XM_011518080.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X4, mRNA
XM_011518081.3	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X8, mRNA
XM_011518082.3	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X7, mRNA
XM_011518085.2	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X27, mRNA
XM_047424166.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X1, mRNA
XM_047424168.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X5, mRNA
XM_047424169.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X6, mRNA
XM_047424170.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X9, mRNA
XM_047424171.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X10, mRNA
XM_047424172.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X11, mRNA
XM_047424173.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X12, mRNA
XM_047424174.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X13, mRNA
XM_047424176.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X14, mRNA
XM_047424177.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X15, mRNA
XM_047424178.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X16, mRNA
XM_047424179.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X17, mRNA
XM_047424180.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X18, mRNA
XM_047424181.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X19, mRNA
XM_047424182.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X20, mRNA
XM_047424183.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X21, mRNA
XM_047424184.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X22, mRNA
XM_047424185.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X23, mRNA
XM_047424186.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X24, mRNA
XM_047424187.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X25, mRNA
XM_047424188.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X26, mRNA
XM_047424189.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X28, mRNA
XM_047424191.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X29, mRNA
XM_047424192.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X30, mRNA
XM_047424193.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X31, mRNA
XM_047424194.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X32, mRNA
XM_047424195.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X33, mRNA
XM_047424196.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X34, mRNA
XM_047424197.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X35, mRNA
XM_047424198.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X36, mRNA
XM_047424199.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X37, mRNA
XM_047424200.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X38, mRNA
XM_047424202.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X39, mRNA
XM_047424203.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X40, mRNA
XM_047424204.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X41, mRNA
XM_047424205.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X42, mRNA
XM_047424206.1	PREDICTED: Homo sapiens maternal embryonic leucine zipper kinase (MELK), transcript variant X43, mRNA
```

Il est possible de visualiser cette liste de transcripts en allant dans la section 
*Genomic regions, transcripts, and products* de la page NCBI du gène [MELK](https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=9833). Attention, dans le menu déroulant *Genomic Sequence:* il est important de bien prendre la référence *GRCh38.p14*.

![image](https://gist.github.com/user-attachments/assets/42b1f4c9-12d9-408d-9f36-8a6cc455c133)

### 2.2. Le génome de référence

Le génome humain de référence est une séquence de l’ADN humain qui sert de modèle standard pour représenter l'ensemble du génome humain. Il s'agit d'une carte génomique unifiée qui combine les séquences d'ADN provenant de plusieurs individus, utilisée pour la recherche, l'analyse clinique et la comparaison génétique.

Tout comme le transcriptome de référence, nous allons également prendre en possession le génome humain de référence.

#### 2.2.1. Téléchargement du génome de référence

- Aller sur la page du génome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Sélectionner l'onglet *FTP*
- Une liste d'Index of /genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14 sera proposée (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)
- Choisir *GCF_000001405.40_GRCh38.p14_genomic.fna.gz*

#### 2.2.2. Informations sur le génome de référence

**Le nom du fichier**

La seule différence avec le fichier du trancriptome de référence c'est qu'il n'y a que la fin du nom de fichier qui change. On passe de *_rna.fna.gz* à *_genomic.fna.gz* qui indique que le format FASTA comprends le génome. 

**Accessions number**

Quand l'on regarde les identifiants des séquences FASTA du génome de référence, nous pouvons constater que les numéros accessions débutent par *NT*, *NC* et *NW*. Nous savons dores et déjà que les ids débutant par *N* sont des *Known RefSeq*. Ainsi, ces séquences sont revues manuellement par des membres du NCBI ou des collaborateur. 
Sinon, si l'on se réfère à la Table 1 du [Chapter 18, The Reference Sequence (RefSeq) Database](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly), nous pouvoir voir : 

> NC_ |	Genomic	| Complete genomic molecule, usually reference assembly

> NT_	| Genomic	| Contig or scaffold, clone-based or WGS

> NW_	| Genomic	| Contig or scaffold, primarily WGS

Prenons des exemples d'identifiants dans le fichier du génome de référence : 
```
>NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly
>NC_000009.12 Homo sapiens chromosome 9, GRCh38.p14 Primary Assembly

>NT_187361.1 Homo sapiens chromosome 1 unlocalized genomic scaffold, GRCh38.p14 Primary Assembly HSCHR1_CTG1_UNLOCALIZED
>NT_187372.1 Homo sapiens chromosome 9 unlocalized genomic scaffold, GRCh38.p14 Primary Assembly HSCHR9_UNLOCALIZED_CTG1

>NW_012132914.1 Homo sapiens chromosome 1 genomic patch of type FIX, GRCh38.p14 PATCHES HG1342_HG2282_PATCH
>NW_013171804.1 Homo sapiens chromosome 9 genomic patch of type NOVEL, GRCh38.p14 PATCHES HSCHR9_1_CTG6
```
Nous pouvons constater avec ces exemples que les identifiants débutant par *NC_* correspondent à des chromosomes

**Quantification des séquences présentes**

Nous pouvons voir que dans le fichier FASTA du génome de référence il y a 705 séquences dont :
- 25 *NC_* soit 24 chromosomes + mitochondrie
- 355 *NT_*, unlocalized genomic scaffold, assemblées mais leur emplacement précis n'est pas déterminé avec certitude
- 325 *NW_*, genomic scaffold, souvent utilisé comme révision ou affinage de région

Ce qui comprends : 
- environ 20000 à 21000 gènes codant pour des protéines
- environ 10000 à 20000 gènes non codants (ARN non codants, microARN, etc.) 

### 2.3. Le GTF

Le GTF humain (General Feature Format) est un fichier contenant des annotations génomiques spécifiques au génome humain. Ces fichiers fournissent des informations détaillées sur la localisation des gènes, des exons, des introns et d'autres éléments fonctionnels dans le génome humain. Le GTF est identique à la version 2 du GFF (General Feature Format). Si nous téléchargeons le GTF et non le GFF c'est parce que le GTF contient les types d'intérêts pour de la transcriptomique à contrario du GFF qui contient les types d'intérêts mais également d'autres non utiles ici.

#### 2.3.1. Téléchargement du GTF

- Aller sur la page du génome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Sélectionner l'onglet *FTP*
- Une liste d'Index of /genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14 sera proposée (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)
- Choisir *GCF_000001405.40_GRCh38.p14_genomic.gtf.gz*

#### 2.3.2. Informations sur le GTF

**La constitution du fichier**

Contrairement aux transcriptomes et génomes de référence, le fichier GTF n'est pas constitué d'identifiants et de séquences mais de colonnes :
- seqid, l'identifiant de la séquence où se trouve la fonctionnalité
- source, l'algorithme ou la procédure qui a généré la fonctionnalité. Il s'agit généralement du nom d'un logiciel ou d'une base de données
- type, le nom du type de fonctionnalité comme "gène" ou "exon"
- start, position génomique de départ de la fonctionnalité
- end, position génomique de fin 
- score, valeur numérique qui indique généralement la confiance de la source dans l'entité annotée. Une valeur de "." (un point) est utilisé pour définir une valeur nulle
- strand, indique le sens avec "+" (positif ou 5'->3'), "-", (négatif ou 3'->5'), "." (indéterminé), ou "?" pour inconnus
- phase, la phase du CDS qui peut être 0, 1, 2 ou "."
- attributes, attributs supplémentaires au format clé-valeur (ex. : gene_id, gene_name)

**Les types présents dans le GTF**

Pour voir tous les types présents avec leurs comptes, il faut lancer cette commande sur le fichier GTF dézippé :

```
awk '{a[$3]++}END{for ( i in a) print i" "a[i]}' GCF_000001405.40_GRCh38.p14_genomic.gtf
```

Voici les informations que nous obtenons :

```
Genomic 22045
 5
RefSeq 1
exon 2298755
CDS 1837053
gene 50329
start_codon 144663
Genomic%2Ccmsearch 1
stop_codon 144508
transcript 200311
```
A titre informatif, voilà ce que nous pouvons en dire :

- Genomic (22045) : Ce type fait généralement référence à des éléments génomiques sans spécification claire (par exemple, des régions non annotées ou des segments de séquences génomiques).
- RefSeq (1) : Ce type fait référence à des éléments spécifiquement annotés selon le système RefSeq, qui est une base de données d'annotations génomiques maintenue par le NCBI (National Center for Biotechnology Information).
- exon (2298755) : Les exons sont les segments de gènes qui sont exprimés et qui sont traduits en protéines. Cette annotation indique combien d'exons sont présents dans le fichier.
- CDS (1837053) : CDS signifie Coding Sequence. Cela représente les régions du gène qui sont traduites en protéine. Cette annotation indique combien de CDS sont présents dans le fichier.
- gene (50329) : Cette annotation indique le nombre de gènes identifiés dans le fichier GTF. Chaque gène est généralement associé à plusieurs exons et éventuellement à plusieurs transcrits.
- start_codon (144663) : Cette annotation représente le nombre de codons de départ (ATG) présents dans le fichier. Un codon de départ marque le début de la traduction d'un ARNm en protéine.
- stop_codon (144508) : Cette annotation représente le nombre de codons d'arrêt présents dans le fichier. Un codon d'arrêt marque la fin de la traduction.
- transcript (200311) : Cela indique le nombre de transcrits différents présents dans le fichier GTF. Un transcrit est le produit de la transcription d'un gène, qui peut inclure des variantes d'épissage.
- Genomic%2Ccmsearch (1) : Ce type semble être une erreur ou une annotation non standard. %2C est une représentation URL pour une virgule, ce qui peut indiquer un problème de formatage ou une entrée non intentionnelle.

## 3. Préparation des lectures

### 3.1. Inférence du strandness

Les librairies RNAseq sont préparées de plusieurs façon ce qui fait que les lectures de séquençage peuvent présenter diverses caractéristiques. Les lectures présentent ainsi :
- une orientation
  - I = inward, vers l'intérieur (fr)
  - O = outward, vers l'extérieur (rf) 
  - M = matching (ff)
- la spécificité de librairie
  - S = stranded
  - U = unstranded
- la spécificité de brin (seulement dans le cas d'une librairie stranded)
  - F = forward, dans le cas où les lectures du fichier 1 dérivent du brin ADN sens et sont donc les copies du brin anti-sens (secondstrand)
  - R = reverse, dans le cas où les lectures du fichier 1 dérivent du brin ADN anti-sens et sont donc les copies du brin sens (firststrand)

![image](https://gist.github.com/user-attachments/assets/5620a856-5538-42fe-afc8-7016566cfe1f)
(voir la page de documentation de [Salmon](https://salmon.readthedocs.io/en/latest/library_type.html))

Les outils RNAseq d'assemblage et de comptage intègrent ces informations via un paramètre. Une utilisation incorrecte de ce paramètre peut avoir un impact sur le résultat des analyses RNA-Seq, notamment la perte de lectures lors d'alignement contre une référence et l'obtention de faux positifs et négatifs dans les résultats d'expression différentielle.

Nous allons donc déterminer les spécifités des librairies avec Salmon (v1.10.1) et le transcriptome de référence pour pouvoir renseigner correctement les paramètres d'outils.
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

Ainsi, nos librairies sont `ISR`. Ce qui veut dire que :
- les lectures pairées vont vers l'intérieur l'une de l'autre, soit une polarité fr
- les lectures sont strand spécifiques
- la lecture 1 est reverse, elle dérive donc du brin ADN anti-sens et est donc la copie du brin sens. La spécifité de brin est firststrand

De ce fait, nous pouvons dire que nos lectures sont **fr-firststrand**.
Selon l'outil utilisé, voici donc l'option à renseigner (https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) : 

![image](https://gist.github.com/user-attachments/assets/9acdcb99-db6a-44ed-82ff-7d8d22087d39)


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
