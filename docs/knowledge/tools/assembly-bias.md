# Compréhension des biais d'assemblage

Lorsque l'on effectue un assemblage de novo genome guided avec un outil comme [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml), il arrive que les résultats ne correspondent pas entièrement aux attentes biologiques. 

Lors d'une analyse de transcrits issus d’un assemblage, nous avons pu observer que tous les transcrits obtenus correspondaient à une fusion de deux gènes (*shble* et *eGFP*) même si le niveau de couverture entre ceux-ci était très inégal. En effet, la couverture de *shble* était nettement plus élevée que celle de *eGFP*. 

![alt text](../../img/igv-cov-transcrits-assembly.png)

Qu'est-ce qui peut expliquer qu'aucun transcrit *shble* seul n'a été identifié, malgré sa forte couverture ? 

## Pourquoi StringTie fusionne-t-il ces transcrits ?
### Un modèle de transcription basé sur la continuité

StringTie en genome guided utilise les alignements des lectures pour reconstruire les transcrits en favorisant les connexions entre régions. *shble* et *eGFP* étant exprimés à partir d'un même locus, car sur un plasmide, StringTie peut supposer qu’ils forment un seul transcrit. Cela reste vrai même si la couverture entre ces régions est inégale. La présence de lectures chevauchant les deux gènes contribue à cette fusion.

### Faible couverture de *eGFP*

Bien que la couverture de *eGFP* soit faible, elle semble suffisante pour établir un lien avec *shble*. StringTie ne discrimine pas les régions à faible couverture si elles sont reliées par des lectures contiguës. Ainsi, tant qu’une continuité existe dans les alignements, StringTie préférera reconstruire un transcrit complet plutôt que de séparer les régions.

### Paramètres de StringTie

StringTie offre de nombreux paramètres configurables, mais certains d'entre eux peuvent accentuer le problème :

- `-f` : ce paramètre contrôle la proportion minimale d'expression requise pour qu'un isoforme soit reporté. Une valeur élevée peut favoriser les transcrits dominants (comme *shble-eGFP*) et exclure des formes minoritaires (comme *shble* seul)
- `--rf / --fr` : ces paramètres permettent de renseigner le strandedness de la librairie. Un mauvais choix peut favoriser des défauts d'assemblage


```
-f <0.0-1.0>	Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts. Default: 0.01

--rf	Assumes a stranded library fr-firststrand.

--fr	Assumes a stranded library fr-secondstrand. 
```

## Comment explorer et résoudre ce problème ?

### Visualisation des alignements

L’analyse visuelle des alignements est un premier pas crucial. Les outils comme [IGV](https://igv.org/) permettent d’observer si des lectures couvrent uniquement shble ou si elles sont toujours connectées à *eGFP*. Cela permet aussi d’évaluer si les régions à faible couverture de *eGFP* influencent l’assemblage.

### Ajustement des paramètres StringTie

Les paramètres de StringTie cités plus haut peuvent être modifiés pour explorer d'autres hypothèses :

- réduire `-f` pour inclure des isoformes minoritaires (même si déjà très faible)
- s'assurer que le strandedness est renseigné correctement

### Assembler différemment

En complément de StringTie, il peut être utile d'utiliser d’autres outils comme [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) pour assembler les transcrits. Ces outils ont des algorithmes différents et pourraient révéler des transcrits absents dans les résultats de StringTie.

## L'impact de la biologie

Il est important de se demander si l'absence de transcrits *shble* seul est une limite technique ou une réalité biologique. Dans certains systèmes, notamment avec des constructions artificielles comme des plasmides, tous les transcrits peuvent naturellement être co-transcrits. Si tel est le cas, il n’y aura biologiquement aucun *shble* isolé même si la couverture suggère qu’il pourrait exister.
C'est là où il est intéressant de compléter la bioinformatique par la paillasse.