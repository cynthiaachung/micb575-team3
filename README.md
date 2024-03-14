# MICB 575 PROJECT 2
**Team 3:** Cynthia Chung, Hannah Hauch, Ellie Kim, Negar Zaghi, Nicole Howes

# Upcoming Meeting Agenda ✨
### March 14, 2024
#### Agenda
* Discussing results of the correlation analyses (Cynthia)
#### Meeting Notes
* Correlation Analyses
  * Missing groups - used as a reference! Remember to specify reference for dependent variable in addition to independent variable
  * Interpretation: reference group is Soil A, for Soil B compared to Soil A, if it's a positive coefficient it means that there's a greater likeihood that the tree-type (for example) will occur in Soil B then Soil A
* Core Microbiome
  * (Ellie asked her questions)

# Archive 🗒️
### March 7, 2024
#### Agenda
* Going over team proposal feedback from Avril
* Troubleshooting analyses (i.e., code, approach, etc.)
#### Meeting Notes
* Avril was confused about soil classification and compactness - should define it better for him
* Move some of the previous paper details into the dataset overview section
* Title should be change from "within" to "and"
* Alpha-/beta-diversity analyses: emphasize the compactness more
  * Look at compaction as its own factor across soil classification
  * If there isn't interesting data, then could be an extra point to look at soil compactness within soil classification as well 
* Core microbiome analyses: describe it like it is in the module
  * Could do it between different compactness as well to strengthen approach
* Indicator taxa analysis: you would need a model to predict, but we're not making one
  * Could add the sentence that Avril suggested
  * Only looking at soil type but not compactness, wanted to use taxonomic data to predict unknown soil type
  * More so on the prediction front (if this is the bacteria population, likely this soil type) - less so of the other, Avril likely didn't pick up on that
  * Rephrase it as we know this unknown soil type, can we infer it using compactness, if we find multiple correlations between compactness and taxonomy, therefore move on to use bacteria taxonomy to predict it
  * Could do it for compactness as well to strengthen argument (?) 
* Differential abundance: Just focus on compactness
* Correlation analyses: (See Bessie's email)
### February 22, 2024
#### Agenda
* Rediscuss project topic and its relation to the dataset
* Confirmation of QIIME2 processing decisions, check trimming parameters...
* Discussion of correlation analysis (Wanted to confirm that it's correlation matrix)
* Discussion of logistical regression analysis (Code? How does it work?)
* Layout of project - is it okay if the aims are outlined by methods as opposed to specific research questions?
#### Meeting Notes
* Hannah asked about primers
  * Based on the paper, the primers that they used contained weird base pairs (i.e., W, M)
  * She ended up using the universal primers for V1-V3
  * Bessie said that there are universal/classic primers and improved primers, likely they used the improved ones
    * Bessie says that it was a good decision to use the universal ones - yay!
* Cynthia asked about correlation analysis
  * Confirmed correlation matrix is what we wanted to do
  * Evelyn hasn't sent the code for the logistical regression analysis, but Bessie will ask her to send it!
* It's okay to have the aim titles as the steps, but it's important to describe why you're doing it (explain the rationale)
  * Bessie: Make sure you understand what the code is for and what it's doing, include that in the proposal (Evelyn really likes to see that!)
* Hannah wanted to clarify our research objective and how that works with the dataset
  * The dataset is looking at LTSP and forest management, but we're looking at how soil classification impacts tree covers > would we have to take this into account?
  * Hannah found a separate paper that said that LSTP has minimal effects on the microbial community, but that there is is an impact
  * Bessie says that that's a fair concern, but something is still impacting the soil diversity so most of the research project is still valid
    * The key concern is how this may impact the correlation analysis
    * The first thought is to eliminate for this confounding factor, you could select for the reference for the correlation analysis
      * Briefly looking at the dataset, the sample size should be okay 
    * Hannah is proposing that for some of the other aims, would it be useful to perform a correlation analysis between soil compaction types 
    * The second thought is that since we know that it's likely because of soil treatment and not tree cover, then we could use that as a segway to other components of our report
* We've confirmed that our overarching research question is what variable(s) impacts soil classification
* Negar asked about if we should filter for alpha/beta diversity and trimming
  * Negar shared the screen with the demux.qzv file
  * Looking at the sequence length summary, it isn't perfect like in the modules (obviously) but it's pretty good at the beginning before dropping down
  * At 409, it goes from a quality of 29 to 26 and stays around here > this is where she chose
  * Bessie says that it looks good, the rest of us agree :)
* Negar asked about the checklist and how it asks us to use the stats.qzv file to justify her trimming choice
  * Bessie doesn't think that it's super necessary > could look at percentages and such, but since we decided it together it should be okay!
  * An alternative is to choose a "bad" one and show the differences
* Ellie asked about rarefaction
  * After filtering the chloroplast and mitochondria, the category she's looking at is soil classification but she was worried about the sampling depth and losing a group within soil classificaiton
  * Looking at the alpha rarefaction curve and histogram, Bessie thinks that 2526 is a good sampling depth so that we can include the ~8 Orthic Gray
* Bessie's key notes
  * Cite the bioinformatics tools you use
  * Detailed figure legends are important (title = takeaway, legend = all the information that the reader needs to know) 

### February 15, 2024
#### Agenda
* Discuss next steps for project
* Clarify objectives for proposal
#### Meeting Notes
* Data wrangling, processing metadata > we haven't done that yet but that's the immediate next steps
  * Metadata wrangling and correlation analysis should be separate aims
  * For both the correlation analysis and the PCA plot, Bessie says that we need numerical data
  * Logistical regression analysis is another route we can go because it works on categorical
    * Define which variables we want to correlate (so tree cover v. soil classification)
    * It's similar to a Chi square test, will tell us how significantly correlated it is
    * Can outline confounding variables and it will take that into account
* When separating samples, it is important to give them new names so that the QIIME2 pipeline won't get confused by it
  * i.e., Tree would become separated into Tree-1, Tree-2, Tree-3 instead of Tree, Tree, Tree
  * Bessie says that it should be easy code to separate
* **Updated Aim 1:** QIIME2 processing
  * Needs to be done by proposal time
  * No filtering and denoising > look at it as is, and then decide how the data should be subsetted for Aims 3-5
* **Updated Aim 2:** Basic alpha- and beta- diversity for both tree cover and soil 
  * Keep the data set as is (no wrangling)
* **Updated Aim 3:** Indicator taxa analysis on the soil
  * It's a form of predictive analysis which will be applicable for our research question
* **Updated Aim 4:** Core microbiome analysis
* **Updated Aim 5:** Differential abundance
* **Updated Aim 6:** Metadata wrangling for correlation analysis, logistical regression analysis
  * Evelyn will send us the code Simran made for the logistical regression analysis
* Other than the checklist, is there anything else that's required for the project proposal > the outline is under the assignment submission page
  * Checklist is only for the dataset 

### February 8, 2024
#### Agenda
* Discuss/confirm final project topic
* Next steps on getting started
#### Meeting Notes
* Soil dataset confirmed - yay!!!
  * Cellular degradation rates and the proteomics associated with it
* Dr. Sun pulled up the metadata to go through it and have a discussion with it
  * Bessie noted that some of the samples in the dataset don't have their soil classification listed down, so we can analyze the microbial taxa from another dataset to do it retrospectively
    * All of the NAs are from the same place, so we could confirm our conclusions in the end by looking through the literature
    * Could also see if certain soil types correlate with tree covers
  * We could filter the dataset to only look at unique soil types (i.e., remove ones that have two)
* This will be our plan of action:
  * **Aim 1:** Correlation analysis to see if there's any link between soil analysis and the trees (independent from QIIME)
    * Filter out all the samples that have mixed soils (will help decrease sample size)
    * For the tree covers, take the ones that have multiple and expand it
    * Can either do PCA or create counts to allow for a correlation matrix
  * **Aim 2:** QIIME2 pipeline process
    * Keep it as it is, don't have to separate it
  * **Aim 3:** General diversity metrics (i.e., alpha/beta) with both soil classication and tree covers
    * Could do in both R and QIIME
  * **Aim 4:** Indicator taxa
    * Can create a predictive model to see if one will predict the other
    * Once we have an established model, we can further segregate the data to see if it will improve (i.e., atmosphere vs. organic, herbicide, etc.)
  * **Aim 5:** Core microbiome
* As we go, we can figure out the balance of having a more robust dataset vs. sample size
* For the project proposal, have the QIIME2 pipeline completed (up to the denoising step)
* When updating the files, have a folder for QIIME and R > each will have the scripts, outputs, etc. 
