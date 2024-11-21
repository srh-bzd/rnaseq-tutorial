# RefSeq

## A propos
[RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) est la base de données NCBI de séquences de référence ; un ensemble organisé et non redondant comprenant des contigs d'ADN génomique, des ARNm et des protéines pour des gènes connus, ainsi que des chromosomes entiers.

## Accession numbers et types de molécules
Les entrées RefSeq sont similaires au format d'entrée de GenBank. Attribuer une nouvelle entrée RefSeq inclus un unique préfixe d'accession suivi d'un underscore.
Le préfixe d'accession RefSeq a une signification implicite en termes de type de molécule qu'il représente. 
Voici tous les préfixes d'accession RefSeq et le type de molécule associé provenant de la Table 1 du  
[Chapitre 18 The Reference Sequence (RefSeq) Database](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly) retrouvé dans [The NCBI Handbook](https://www.ncbi.nlm.nih.gov/books/NBK21101/)

| Accession prefix | Molecule type                                                | Comment                                                      |
| :--------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
| AC_              | Genomic                                                      | Complete genomic molecule, usually alternate assembly        |
| NC_              | Genomic                                                      | Complete genomic molecule, usually reference assembly        |
| NG_              | Genomic                                                      | Incomplete genomic region                                    |
| NT_              | Genomic                                                      | Contig or scaffold, clone-based or WGSa                      |
| NW_              | Genomic                                                      | Contig or scaffold, primarily WGSa                           |
| NZ_             | Genomic                                                      | Complete genomes and unfinished WGS data                     |
| NM_              | [mRNA](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app114/) | Protein-coding transcripts (usually curated)                 |
| NR_              | [RNA](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app158/) | Non-protein-coding transcripts                               |
| XM_             | [mRNA](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app114/) | Predicted model protein-coding transcript                    |
| XR_             | [RNA](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app158/) | Predicted model non-protein-coding transcript                |
| AP_              | Protein                                                      | Annotated on AC_ alternate assembly                          |
| NP_              | Protein                                                      | Associated with an NM_ or NC_ accession                      |
| YP_             | Protein                                                      | Annotated on genomic molecules without an instantiated transcript record |
| XP_             | Protein                                                      | Predicted model, associated with an XM_ accession            |
| WP_              | Protein                                                      | Non-redundant across multiple strains and species            |

## Code status
Chaque enregistrement a un *COMMENT*, indiquant le niveau de conservation qu'il a reçu et l'attribution du groupe collaborateur. 
Ainsi, un enregistrement RefSeq peut être une copie essentiellement inchangée et validée de la soumission originale de l'INSDC, ou inclure des informations mises à jour ou supplémentaires fournies par des collaborateurs ou le personnel du NCBI.
Voici tous les codes status RefSeq provenant de la Table 2 du  
[Chapitre 18 The Reference Sequence (RefSeq) Database](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly) retrouvé dans [The NCBI Handbook](https://www.ncbi.nlm.nih.gov/books/NBK21101/)

| Code        | Description                                                  |
| :---------- | :----------------------------------------------------------- |
| MODEL       | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record is provided by the [NCBI](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app116/) Genome Annotation pipeline and is not subject to individual review or revision between annotation runs. |
| INFERRED    | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record has been predicted by genome sequence analysis, but it is not  yet supported by experimental evidence. The record may be partially  supported by homology data. |
| PREDICTED   | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record has not yet been subject to individual review, and some aspect of the [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record is predicted. |
| PROVISIONAL | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record has not yet been subject to individual review. The initial  sequence-to-gene association has been established by outside  collaborators or [NCBI](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app116/) staff. |
| REVIEWED    | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record has been reviewed by [NCBI](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app116/) staff or by a collaborator. The [NCBI](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app116/) review process includes assessing available sequence data and the literature. Some [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) records may incorporate expanded sequence and annotation information. |
| VALIDATED   | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record has undergone an initial review to provide the preferred  sequence standard. The record has not yet been subject to final review  at which time additional functional information may be provided. |
| WGS         | The [RefSeq](https://www.ncbi.nlm.nih.gov/books/n/handbook/A1237/def-item/app155/) record is provided to represent a collection of whole genome shotgun  sequences. These records are not subject to individual review or  revisions between genome updates.
