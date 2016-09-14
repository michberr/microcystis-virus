# CRISPRs from Lake Erie *Microcystis*








# Microcystis spacer analyses 
These are analyses from the extracted spacers from the two Microcystis DR types








## Basic spacer statistics



<img src="Figs/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

<img src="Figs/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />
The maximum coverage of any spacer is 460. It is found in sample 
E0108-CNA
      
There are 455 spacers with more than 10x coverage


## Spacer alpha diversity 

#### Number of spacers in each sample

The number of spacers associated with the two Microcystis repeats vary a lot over time and fractions. The highest number of spacers are found during the two peak bloom dates in the 100um fraction. Is this a good repesentation of higher spacer richness, or just better sampling depth because Microcystis is more abundant? Relative abundances are not that different between 100um and 53um fractions on the peak dates, so i think the lowered spacer richness is a real result.

<img src="Figs/unnamed-chunk-6-1.png" style="display: block; margin: auto;" /><img src="Figs/unnamed-chunk-6-2.png" style="display: block; margin: auto;" />
Here I calculated the inverse simpson index of spacers for each DR type in each sample using coverage
<img src="Figs/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />


## Spacer overlap over time

How much spacer overlap is there over time? 

There is a lot of spacer overlap between the two peak bloom dates in the 100um fraction - about 400 spacers or half of the total spacers in those samples. Need to remember that E0048 (august 4) has two DR types but E0108 (Sep 29) has only one DR type
<img src="Figs/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />



# Virome spacer hits

These are analyses of the BLAST hits between viral contigs from our viromes and CRISPR spacers from our microbial metagenome





<img src="Figs/unnamed-chunk-11-1.png" style="display: block; margin: auto;" /><img src="Figs/unnamed-chunk-11-2.png" style="display: block; margin: auto;" />

<img src="Figs/unnamed-chunk-12-1.png" style="display: block; margin: auto;" /><img src="Figs/unnamed-chunk-12-2.png" style="display: block; margin: auto;" />

Hits to the same contig

```
## Source: local data frame [1 x 6]
## Groups: metaGDate, metaGFraction, virome, DR [1]
## 
##   metaGDate metaGFraction virome    DR   VirContig  hits
##       (chr)        (fctr)  (chr) (chr)       (chr) (int)
## 1     E0108           100  E0108     T contig_5127    11
```

# Literature

#### Rath et al 2015 Biochimie 

- length of DR's are 21-48     
- length of spacers are 26-72    
- number of spacers vary from a few to several hundred    
- Genomes can have single or multiple CRISPR loci

Methanocaldococcus has 18 CRISPRs and 191 spacers - makes up 1% of genome    

Not all CRISPR loci have adjacent cas genes - instead rely on trans encoded factors
leader - conserved sequence located upstream of CRISPr

**Cas proteins:**     
- Cas1 and Cas2 are virtually universal, others are not      
Type 1, II and III Crispr-cas systems, newly proprosed Type IV complex which contains Cascade genes but no crispr, cas1 or cas2. Guided by protein-DNA interaction, not by crRNA. Type IV could be ancestral innate immune system that became adaptive by associating with transposon like element containing cas1 and cas2

CRISPR-cas are probably mobile genetic elements that tranfer frequently horizontally

AMD system where no two cells had the same spacers, viruses used recombination to diversify rapidly, making any but the most recent spacers ineffectual



small number of spacers dominates but relative abundance oscillates rapidly [29]

**Primed spacer acquisition:**    
phages can still replicate in populations with one but not two spacers targeting them? see [30]

Adaptation has been shown to be coupled to the interference machinery through primed spacer acquisition, which occurs when there is a targeting spacer already present in the CRISPR array. The interference machinery and a pre-existing spacer accelerate the acquisition of subsequent spacers from the same target. Primed spacer acquisition was first described in the Type I-E system in E. coli [37], but has subsequently been reported for I-B in H. hispanica [34] and I-F in P. atrosepticum [39], but **so far not in any Type II or III system**. Priming seems to occur by slightly different processes in the described cases but the exact molecular mechanisms remain unknown. In Type I-F systems, Cas2 is fused to Cas3 [13], further indicating a direct connection between the adaptation and interference processes. Interestingly, spacers with several mismatches that are incapable of providing protection against the target still induce primed spacer acquisition [49]. It should be noted that although Cas9 is required for spacer acquisition in the Type II-A system, this is not an example of primed spacer acquisition as the requirement is not dependent on a pre-existing spacer against the target [20]. The advantages of primed spacer acquisition are obvious: multiple spacers provide increased resistance against invading DNA, and make it more difficult for target to evolve escape mutants as several sites would need to be changed simultaneously.     

In an unusual turn of events, it has been demonstrated that the CRISPR-Cas system can be used by viruses to promote infection. Vibrio cholerae ICP1 phages carry a Type I-F CRISPR–Cas system that targets a host locus, PLE, containing an anti-phage system [89]. After entry of the phage genome into the cell, the viral crRNAs and cas genes are expressed to enable infection of the V. cholerae host. If the host or the CRISPR is engineered so that the viral CRISPR-Cas system no longer matches the PLE, the ability of ICP1 to infect is largely lost. The few phages that manage to infect do so by having picked up a new spacer targeting the host locus, demonstrating that the viruses can use the full adaptive potential of the CRISPR-Cas system.



**Knowledge gaps: **      
- how are protospacers selected? They are next to a PAM, cas9 plays a role in identifying PAMs       
- not known whether protospacer is copied or cut out of the target      
- when and why are spacers deleted?       
- How does primed spacer acquisition occur ... link to adaptiveness     



