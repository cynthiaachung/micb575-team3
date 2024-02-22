# MICB 575 PROJECT 2
**Team 3:** Cynthia Chung, Hannah Hauch, Ellie Kim, Negar Zaghi, Nicole Howes

# Upcoming Meeting Agenda ✨
### February 22, 2024
#### Agenda
* Rediscuss project topic and its relation to the dataset
* Confirmation of QIIME2 processing decisions
* Discussion of correlation analysis (Spearman's Rank?)
* Discussion of logistical regression analysis (Code? How does it work?)
* Layout of project - is it okay if the aims are outlined by methods as opposed to specific research questions?

# Archive 🗒️
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
  * Evlyn will send us the code Simran made for the logistical regression analysis
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
