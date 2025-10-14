## Reviewer 1

In this paper, Chen and Xie et al present a novel single cell Hi-C method which they term HiChew. Single-cell Hi-C strategies generally fall into two categories: enrichment-based methods (e.g., standard scHi-C), which achieve high valid pair ratios but suffer from reduced library complexity due to DNA loss, and non-enrichment approaches (e.g., Dip-C), which preserve library complexity but yield very low valid pair ratios, necessitating deep sequencing. HiChew is pitched as a good middle ground between these strategies, where a dam methylase is used to mark ligation scars post-pcr before immunoprecipitation to isolate ligated read-pairs. They present evidence of increased efficiency (valid pair ratio) compared to Dip-C but more complexity than classic enrichment-based methods, allowing them to achieve up to 5kb resolution with very low cell numbers. They further validate HiChew by analysing B-compartment melting, CTCF depletion, and chromatin conformation in a primary cell type.

Overall, this work is impressive, combining a novel enrichment strategy with comprehensive comparison data to existing methods. There is also application to interesting biology at the end. I believe HiChew adds value to the growing list of single-cell HiC methods due to its consistent mid-range resolution with exceedingly low cell numbers. The majority of the concerns from previous revision rounds for which I was not involved in, seem to have been addressed. However, the paper is quite difficult to read in certain sections, and I believe there a few areas where clarification or additional analysis would strengthen the manuscript further.

### Major concerns:

·         The CTCF experiments and its relationship to melting are really interesting but I'm not sure about the validity of using the heterozygous knockout cell line.  Have the authors considered validating these findings with an acute depletion system (e.g., dTAG or siRNA)? If not feasible, further justification for the choice of model would be important. This leads onto my next points.
<br> <u>wetlab data to show siRNA is not efficient</u> <br>

·         Can we see the raw genome browser tracks of the CTCF ChIP-seq at selected peaks where it is claimed the signal is lower? The meta-plot is a good sanity check, but to be more convinced I would like to see some representative examples. Furthermore, given it’s quite a subtle altercation, were spike-in controls used?
<br> `Fig 4a label CTCF ChIP-seq WT/KD fold change and move to bottom` <br>

·         Can you quantify the degree of CTCF reduction at all the CTCF peaks and then use this to ask:  (i) Are genes which are transcriptionally dysregulated in CTCF -/+ likely to be closer to a reduced CTCF bound peak and  (ii) are those genes more/less likely to have a greater KS melting score. Rather than quantifying change in CTCF binding, you could also use number of CTCF peaks reduced within a specific window around the gene.
<br><br> `chip /research/xieyeming1/proj_2023/zhichao_ctcfkdChip_20240805/align/rep_map_old/peak/wt_kd_peak.pdf` 
<br> `ks /research/xieyeming1/proj_2023/hichew_paper_20230710/hichew_NM_resubmit_figs/code/ploting_scripts/long_gene_melting/final_ctcfkd` 
<br> `rna /research/xieyeming1/proj_2023/zhichao_RNAseqWTKD_2024/BGI_pipe` 

·         Although it looks like CTCF reduction impacts a large number of genes in terms of melting (563 out of 752) - only 62 of these (roughly 10%) are affected transcriptionally. Can the authors comment on this and hypothesise about the disconnect between melting and transcriptional output? Perhaps the analysis above may shed some light on this.
<br><br> `reduced melting + ctcf, 62 genes` <br>

### Minor points:

·         It would be nice to have a bit more data on the loop calling for HiChew compared to DipC/HiC. Are the loops called largely shared between the methods. What is the breakdown of loops in terms of Enh-prom, Enh-CTCF, prom-CTCF etc? Is this consistent with Dip-C?
<br>`combine fq, run hic_pro, /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulkHiChew_5k_500_50/raw_data` <br>
<br>`R intersect E-P, E-C, P-C, loop bin, /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulk_biotin_hic/HEK293T/hic_pro/bulkHiC_HEK293T/outdir/hic_results/matrix/bulkHiC_HEK293T/iced/5000 `<br>

·         The authors convincingly demonstrate HiChew’s improved performance in terms of valid pair ratio (compared to Dip-C) and complexity compared to more classic enrichment based Hi-C approaches. However, in all 3 application experiments (B compartment melting, CTCF depletion, and chromatin conformation in spermatocytes) I fail to see what biological insight HiChew uncovers compared to Dip-C or enrichment based HiC? Is the resolution better? Are the cell numbers going into these experiments lower? It would be good if this was clear, either in a figure or in the text.
<br> `emphasis: better resolution to uncover changes in TAD at 25k resolution and applicable in testis` <br>

·         Many recent Hi-C and Micro-C protocols employ double crosslinking (e.g., DSG + formaldehyde; see Oksuz et al., Nat Methods 2021), which reduces spurious trans interactions. Could the authors comment on whether such an approach has been tested in HiChew?
<br> `emphasis: double crosslinking in testis` <br>

·         The theoretical complexity arguments for HiChew compared to Dip-C are interesting – because in theory you would expect Dip-C to have higher complexity with unlimited sequencing (assuming losses of contact pairs during methylation and immunoprecipitation). Is there a reason why authors only performed the sequencing saturation analysis on the 50-cell dataset? What do these curves look like at the different cell numbers? Given HiChew’s more efficient performance compared to Dip-C (valid pair ratio per read depth), I don’t believe it needs to have the same level of theoretical complexity to be useful, I am just intrigued by that result.
<br> `sat curve bulk hichew vs pcr, 5k, 500, 50 /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulkHiChew_5k_500_50/` <br>

·         The correlation between digestion efficiency and valid pair ratio makes sense in the context of the dam methylase. Does the 16h tested represent full digestion of the chromatin or could the time be extended further yielding even better valid pair ratios? Many 3C methods now utilise a qPCR based estimation of digestion efficiency (see Downes et al, Nature protocols, 2022) – this may help with optimising new cell types in the future.
<br> <u>explain with wet lab data?</u> <br>

·         CNTNAP2 expression is said to be 1.5 fold increased on line 402, but 2 fold increased on line 417 – is this a different RNA-seq data set?  Please clarify.
<br> `cntnap2 1.5 fold is bulk RNA exp, 2 fold is single cell RNA exp` <br>

Grammer/style points:
·         Line 59 "Lin developed" should be Lin et. al. or Lin and colleagues.
<br> `rewrite` <br>

·         Line 62 MDA abbreviation is used before its definition on line 65.
<br> `rewrite` <br>

·         Line 180 - “Worth noting is that we should steer clear of the GATC motif when designing the barcode”. Rephrase to something like “It is worth noting that, we avoid the GATC motif when designing etc.”
<br> `rewrite` <br>

·         This manuscript is dense with detail in both the comparison sections (bulk and single cell). I wonder if a summary table or schematic of the key comparison metrics (valid pair ratio, theoretical complexity limits, cell numbers etc.) may help to get across the improvements of HiChew compared to Dip-C and enrichment Hi-C to the reader.
<br> `add key metrics to Extened Data Fig.1c` <br>

·         Furthermore, for readability, I think the section “evaluation of factors influencing snHiChew/HiChew efficiency” could be condensed or made into a supplementary note. It is interesting and useful data, but I think it’s placement in the middle of the manuscript doesn’t help the flow.
<br> `rewrite` <br>

## Reviewer 2

The authors present Hi-Chew, a Hi-C variant that employs DpnII digestion, methylation labeling at GATC sites, antibody pulldown, and a streamlined library preparation workflow. Hi-Chew captures more chromatin contacts per cell than Hi-C or Dip-C, thereby improving resolution and reducing sequencing costs. Comparisons between bulk and low-input datasets demonstrate strong concordance. At the single-cell level (snHi-Chew), the authors use split-pool barcoding. The manuscript highlights three applications: analysis of “TAD melting” associated with gene expression in B compartments, the impact of CTCF knockdown on melting, and spermiogenesis. Overall, the work provides convincing support for the efficacy of Hi-Chew and snHi-Chew. The applications are interesting and largely well supported, though several areas could be strengthened.

### Major comments

1. The section on melting in B compartments could be improved.
- The focus on the CNTNAP2 locus is not well motivated, and it is unclear whether this represents a typical region or a cherry-picked example. A more systematic approach would be preferable: analyze all long genes in B compartments, rank them, and then present CNTNAP2 as a top candidate alongside general results.
- It would also be important to check whether technical artifacts could account for the apparent clustering. For example, do QC metrics have similar distributions between the two clusters?
- At present the analysis considers one gene at a time, but embedding genome-wide contact data into a latent space and testing whether the identified “melting clusters” align with latent space clustering patterns could provide evidence of coordinated states rather than random effects.
- While the observation of TAD melting is not surprising in itself, the key biological question is how genes in repressive compartments become activated, and the single-cell data should be leveraged more directly to shed light on this.
<br>`rank gene; melt vs conc QC; melt vs conc 25k counts PCA; ctcf (reviewer 1)`<br>

2. The spermatogenesis analysis is interesting but one point is not well explained. The statement that pachytene spermatocytes demonstrate a lower X/Y chromosome percentage relative to haploid round spermatids is not fully convincing as presented, and would benefit from a clearer rationale or additional supporting evidence.
<br> `rewrite` <br>

3. Figure presentation could also be improved. Figure 5c is difficult to interpret with the current color choices, and colors that are more easily distinguished should be used. The observation of different interaction distance distributions is potentially interesting, but to establish its significance the authors should demonstrate consistency across independent biological replicates.
<br>`plot fig5c, two replicates adjust color /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_testis/benchmark`<br>

4. The testis cell type analysis would be strengthened by examining whether genes expressed in the respective cell types show compartment switching in the expected direction. This type of association is often observed and would provide an important cross-check.
<br>`the diploid to haploid transformation may not correlated with gene expression`<br>

## Reviewer 3

This study from Chen, Xie et al. present a new method for single-cell 3D genome analysis that the authors term HiChew. The basic premise of this is that they use a DNA methyltransferase to mark the sequence targeted by the restriction enzyme in an otherwise single cell Hi-C experiment. The value is that this allows for efficient sticky ended ligation while also being able to enrich for ligation products using an antibody against DNA methylation. This is a nice idea and it does hit on what I think is a clear barrier to the wider adoption for single-cell Hi-C methods, namely that you waste so much money sequencing DNA that isn’t informative for 3D genome organization. Perhaps the strongest statement I can make regarding publication of this method is that we will likely try to use this in our lab going forward for single-cell Hi-C experiments. In this sense, I think the methodological aspect of this study is interesting. I’m less convinced by the biological applications. They focus on a phenomenon they term “melting” of chromatin, but this is poorly defined. They also argue that this is linked to active transcription states but this is not well demonstrated either. I would frankly find this study more suitable for publication if many of the claims regarding novel biology in the second half of the manuscript were toned down or removed, as the method aspect I think it nice and will be of use. I do think that this study requires substantial revisions before it would be appropriate for publication. I have divided my comments below into major and minor comments:
### Major comments:
The authors uses HiC-Pro to characterize read pairs between Dip-C and HiChew (Ext. Fig. 1). To the best of my knowledge, HiC-Pro uses restriction enzyme cut site to characterize read pairs, and if reads align to two different restriction fragments these are characterized as valid pairs. The problem is that this cannot tell if this is a real ligation even or a read that crosses an uncut digestion site. I think that in addition to HiC-Pro, they should also simply quantify what percentage of reads align to the genome either to different chromosomes or greater than an arbitrary threshold (i.e. >1kb or >10kb) as an alternative measure of what fraction of reads are true “contacts”.
<br>`add hic pro QC of cis_long, cis_short, trans`<br>
<br>`/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulkHiChew_5k_500_50/hic_pro`<br>

The relationship between “melting” and expression is weak. I think to really claim this they would need to do some kind of DNA+nascent RNA FISH to show that the changes in spatial distance at the gene are associated with changes in active transcription. 
<br>`add fish data`<br>
The clearest experiment they have to this effect is the triptolide experiments in Figure 4g, but the effect on the insulation scores is very weak, so it isn’t clear that transcriptional is really affecting the “melting”. Further, these effects don’t really account for the single-cell nature of the experiments, as the RNA-seq analysis is using bulk RNA and not single cell, so knowing that the expression is linked to the specific subset of cells with changes in melting is not demonstrated. 
<br> `limca /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/LiMCA/GM12878_RNA` <br>

The word “melting” to me is vague and imprecise. Chromatin isn’t literally melting. I think it would be worthwhile to define this better. My impression is that it basically can be interpreted as “less insulation”. Along the same lines, there are relatively few focal “dips” in the insulation score apparently in Figure 3b,c at the CNTNAP2 gene, and the changes in insulation score (Fig. 3c) appear to be more broadly distributed, so it is not clear to me what this effect is. Is this a weakening of intra-TAD contacts?
<br> `rewrite` <br>

How are the single-cell clustering analysis performed? For Figures 3c, 4b, and 4g, is this done on the local insulation scores shown in the figure? This is fine, but the issue to me is that some of the differences appear to be driven by the insulation scores in the adjacent A compartment regions (148Mb-152Mb).
<br> `cluster by ins in ROI; show the melting in 145-148 alone in suppl` <br>

### Minor comments:
In the introduction, the authors write: “However, biotin enrichment resulted in significant DNA loss, most likely due to inefficient biotin incorporation or loss of biotin binding”. They cite Belton et al 2012 here to justify this statement. This seems odd. First, the Belton paper was written in 2012 before the first single-cell Hi-C experiments were published, so it is not really the right source for a statement regarding single-cell Hi-C. Second, I don’t think that failure of biotin incorporation or binding is the problem for single cell, it is much more that the blunt end ligations are very inefficient. At the very least, it isn’t well demonstrated here what the limitations of biotin incorporation really are here, so I don’t think this sentence is really appropriate.
<br> `rewrite` <br>

Related to the above comment, they also say that ChIA-PET was introduced as a means of circumventing this low biotin incorporation using the biotinylated oligos, but this isn’t true. ChIA-PET was really focused on targeted enrichment associated with specific TFs or chromatin bound proteins.
<br> `rewrite` <br>

One general comment regarding the introduction (lines 69-71), for methods such as snm3C-seq, the lack of enrichment of ligated fragments isn’t actually such a bad thing as these regions can then give genome wide information regarding methylation status.
The authors report (line 120) that conventional Hi-C libraries have valid pair rates of ~95%. That is really high. I don’t know if I have seen much data where it is that high. They need to cite a reference here supporting this claim.
They reference Figure 3i, but there is no panel for this. It seems this is referring to 3g or 3h?
<br> `rewrite, should be Fig. 3g` <br>
