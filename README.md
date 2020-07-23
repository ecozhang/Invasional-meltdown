# Invasional-meltdown

**Soil-microbes-mediated invasional meltdown in plants**\

**Zhijie Zhang, Yanjie Liu, Caroline Brunel and Mark van Kleunen**\
**Nature Ecology and Evolution, 2020.**\

**Codes, reports and data**




Please contact Zhijie Zhang (zhijie.zhang@uni-konstanz.de) if any questions.


We experimentally tested whether and how a third plant species affected the competitive outcomes between alien and native plants through its soil legacy. We first conditioned soil with one of ten species (six natives and four aliens; these are the soil-conditioning species) or without plants. Then, we grew on these 11 soils, five aliens and five natives (these are the test species) without competition, or with intra- or interspecific competition.



1.	The files started with ‘01’ are the **codes**. 
    - ‘01functions.R’ is the functions that are used in ‘01plant_NEE.Rmd’ and ‘01soil_NEE.Rmd’. 
    - ‘01plant_NEE.Rmd’ is the code for plant data, and it calls ‘03plant_dat.csv’ (see below for details of ‘03plant_dat.csv’). If you want to re-run the analyses, please run ‘01plant_NEE.Rmd’ before ‘01soil_NEE.Rmd’, as it will export data that will be used in ‘01soil_NEE.Rmd’
    - ‘01soil_NEE.Rmd’ is the code for the soil-microbe data, and it calls ‘04phyloseq_16s_its.RData’ and ‘05dfr_funguild_guilds.csv’ (see below for details).
2. The files started with ‘02’ are the reports of ‘01xxx.Rmd’ **I think, reading these two reports, after downloading, will solve most questions**. Or, you can simply visit https://ecozhang.github.io/Invasional-meltdown/02plant_NEE for plant results and https://ecozhang.github.io/Invasional-meltdown/02soil_NEE for soil-microbe results.
3. The file ‘03plant_dat.csv’ is the **biomass data of test plants**.
   - The variable ‘fix’ is identity of unit (pot), i.e. the plants with same ‘fix’ number were grown in the same pot.
   - The variables ‘sp_p1’, ‘family_p1’ and ‘origin_p1’ are the species, family and origin (native or alien, or empty if no soil-conditioning plant) of the soil-conditioning plant, respectively.
   - The variable ‘bio_p1’ is the aboveground biomass of soil-conditioning plant when harvest, and in ‘bio_p1_0’, we set the biomass as 0 for pots without plants.
   - The variables ‘target’, ‘family’ and ‘origin_p2’ are the species, family and origin of the focal test plant, respectively.
   - The variable ‘biomass’ is the aboveground biomass of focal test plant, and ‘root’ is the belowground biomass of focal test plant that was grown alone.
   - The variables ‘comp’ and ‘family_comp’ are the species and family of competitor test plant, respecitively.
   - ‘competition’ describe whether the focal test plant was grown alone, or with intra- or interspecific competition.
   - The variable ‘trans_date_c’ is the transplanting date [day as the unit]. The plants that are transplanted first are set as 0.
   - The variable ‘comb’ is the ID of the test species pair (e.g. A1 indicates that Dactylis glomerata was grown alone, H15 indicates that Plantago media and Lolium multiflorum were grown together).
4. The file ‘04phyloseq_16s_its.RData’ include two phyloseq objects of 16S (bacteria) and ITS (fungi). It contains the cleaned **sequencing data** for each soil sample, which was sampled when we harvested the soil-conditioning plants. We used three tables of each phyloseq object:
   - otu_table: a sample-by-ASV [Amplicon sequence variant] matrix, with number of reads as entries.
   - tax_table: a table of the taxonomy of each ASV.
   - sam_data: a table of sample information.
     - ‘fix’ can linked with the ‘03plant_dat.csv’.
     - ‘status’ is the origin of the soil-conditioning plant (n: native, a: alien, ck: no soil-conditioning plant).
     - ‘ster.not’ indicate whether soil inoculum was sterilized or not before the soil-conditioning.
5. The file ‘05dfr_funguild_guilds.csv’ is the FUNGuild (http://www.funguild.org/) information of all ITS ASVs, which is used to determine the functional groups of fungi (according to the ‘Guild’ variable.). The FUNGuild  was accessed on 2019-10-30.
6. The file ‘06tree.tre’ is the phylogenetic tree of all species in the project.
