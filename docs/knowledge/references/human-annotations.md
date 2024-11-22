# Les annotations du génome humain

Le GTF humain (General Feature Format) est un fichier contenant des annotations génomiques spécifiques au génome humain. Ces fichiers fournissent des informations détaillées sur la localisation des gènes, des exons, des introns et d'autres éléments fonctionnels dans le génome humain. Le GTF est identique à la version 2 du GFF (General Feature Format). La différence entre le GTF et le GFF est le fait que le GTF contient les types d'intérêts pour de la transcriptomique.

Le fichier des annotations du génome assemblé [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) est disponible à l'adresse [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/) sous le nom de fichier *GCF_000001405.40_GRCh38.p14_genomic.gtf.gz*.

## La constitution du fichier

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

## Les types présents dans le GTF

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