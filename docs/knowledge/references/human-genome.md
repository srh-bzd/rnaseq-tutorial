# Le génome humain

Le génome humain de référence est une séquence de l’ADN humain qui sert de modèle standard pour représenter l'ensemble du génome humain. Il s'agit d'une carte génomique unifiée qui combine les séquences d'ADN provenant de plusieurs individus, utilisée pour la recherche, l'analyse clinique et la comparaison génétique.

La dernière version du génome de référence humain à ce jour est [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). Il est possible de la télécharger directement à l'adresse [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/) sous le nom de fichier *GCF_000001405.40_GRCh38.p14_genomic.fna.gz*.

Nous allons voir ensemble les informations sur ce génome de référence.

## Le nom du fichier

Nous allons décortiquer le nom de fichier du génome : 

- *GRCh38* (Genome Reference Consortium Human Build 38) correspond à la version du génome humain assemblé
- *p14* (patch release 14) correspond à la version 14 du patch, soit à la 14ème correction
- *GCF_000001405.40* est l'accession RefSeq (l'équivalent en accession GenBank est GCA_000001405.29)
- *_genomic.fna.gz* correspond au format FASTA comprennant les ADN présent dans le génome

## Accessions number

Quand l'on regarde les identifiants des séquences FASTA du génome de référence, nous pouvons constater que les numéros accessions débutent par *NT*, *NC* et *NW*. Les ids RefSeq débutant par *N* sont des *Known RefSeq*. Ainsi, ces séquences sont revues manuellement par des membres du NCBI ou des collaborateur. 
Sinon, si l'on se réfère à la Table 1 du [Chapter 18, The Reference Sequence (RefSeq) Database](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly), nous pouvont voir : 

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

## Quantification des séquences présentes

Nous pouvons voir que dans le fichier FASTA du génome de référence il y a 705 séquences dont :

- 25 *NC_* soit 24 chromosomes + mitochondrie
- 355 *NT_*, unlocalized genomic scaffold, assemblées mais leur emplacement précis n'est pas déterminé avec certitude
- 325 *NW_*, genomic scaffold, souvent utilisé comme révision ou affinage de région

Ce qui comprends environ : 

- 20000 à 21000 gènes codant pour des protéines
- 10000 à 20000 gènes non codants (ARN non codants, microARN, etc.) 