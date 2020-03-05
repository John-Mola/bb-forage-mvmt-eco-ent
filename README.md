# Bumble bee foraging movements in montane meadow complex

![]("./graphical_abstract.png")

This repository is intended to document the code behind the project "Forests do not limit bumble bee foraging movements in a montane meadow complex". As of now this manuscript is in revision at Ecological Entomology. 

Anyone who downloads this repository should be able to reconstruct the results of the analysis following COLONY assignment. Should you wish you reconstruct from the beginning, you can download the raw genotype files at [DRYAD DATA REPO]. For most folks, I think this should cover what you're interested in reviewing/referencing for your own future needs. 

Overall, if you run all of the scripts in the "analyses" directory in order (02a, 02b, 02c, etc) it will create a series of files saved into "analyses_output". You can then run the R markdown file within "docs" to generate a summary of all the figures, model outputs, and odds and ends needed to review the results used in the manuscript. To conduct a clean run, you could delete all of the files in analyses_output, then run all of the analyses, then run the markdown file. 

Please note, this code is not written by a master of programming ability, but by a bumble bee expert who uses code as a means to an end. All of the code is functional, but it may not represent the best way to do something. Should you have any questions or encounter issues, please contact me at jmola@usgs.gov. 

## Description of analysis scripts

1. analyses/02a_capwire_colony_abundance.R - wrangling and analysis associated with the Capwire estimates of colony abundance (used in Table 1)

2. analyses/02b_recaptured_tagged_vosnesenskii.R - matching recaptured tagged bumble bees back to their original location and finding the separation distance

3. analyses/02c_sibling_separation_colony_distance.R - code associated with "Sibling separation distance and colony-specific foraging range" sections. (Generates Figure 3; and subsampling used in Supplemental)

4. analyses/02d_sibling_plant_usage_distance.R - code associated with "Plant use across distance" sections. (Generates Figure 4)

5. analyses/02e_habitat_randomization.R - code associated with the first half of the "Habitat connectivity" section. (Generates Figure 5)

6. analyses/02f_habitat_dissimilarity.R - code associated with the second half of the "Habitat connectivity" section. (Generates Figure 6)

7. analyses/02g_plant_captures_supplemental.R - creates the supplemental table of plant captures