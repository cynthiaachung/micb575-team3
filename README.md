# MICB 575 PROJECT 2
**Team 3:** Cynthia Chung, Hannah Hauch, Ellie Kim, Negar Zaghi, Nicole Howes

# Upcoming Meeting Agenda ✨
### February 14, 2023
* 

# Archive 🗒️
### February 8, 2023
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
