# L'inférence du strandedness en RNAseq

Les librairies RNAseq sont préparées de plusieurs façon ce qui fait que les lectures de séquençage peuvent présenter diverses caractéristiques. Les lectures présentent ainsi :

une orientation

- I = inward, vers l'intérieur (fr)
- O = outward, vers l'extérieur (rf) 
- M = matching (ff)

la spécificité de librairie

  - S = stranded
  - U = unstranded

la spécificité de brin (seulement dans le cas d'une librairie stranded)

  - F = forward, dans le cas où les lectures du fichier 1 dérivent du brin ADN sens et sont donc les copies du brin anti-sens (secondstrand)
  - R = reverse, dans le cas où les lectures du fichier 1 dérivent du brin ADN anti-sens et sont donc les copies du brin sens (firststrand)

![image](https://gist.github.com/user-attachments/assets/5620a856-5538-42fe-afc8-7016566cfe1f)
(voir la page de documentation de [Salmon](https://salmon.readthedocs.io/en/latest/library_type.html))

A titre d'exemple, une librairie ISR indique que :

- les lectures pairées vont vers l'intérieur l'une de l'autre, soit une polarité fr
- les lectures sont strand spécifiques
- la lecture 1 est reverse, elle dérive donc du brin ADN anti-sens et est donc la copie du brin sens. La spécifité de brin est firststrand
 
Les outils RNAseq d'assemblage et de comptage intègrent ces informations via un paramètre. Une utilisation incorrecte de ce paramètre peut avoir un impact sur le résultat des analyses RNA-Seq, notamment la perte de lectures lors d'alignement contre une référence et l'obtention de faux positifs et négatifs dans les résultats d'expression différentielle.

Selon l'outil utilisé, voici l'option à renseigner (https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) : 

![image](https://gist.github.com/user-attachments/assets/9acdcb99-db6a-44ed-82ff-7d8d22087d39)