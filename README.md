# HIV-phyloTSI
Estimate time since infection (TSI) from HIV deep-sequencing data. Amplicon-based and capture-based sequencing are supported. This code requires sequencing reads for each sample to be mapped using [shiver](https://github.com/ChrisHIV/shiver/ "shiver"), and the resulting BAM files batched and processed by [phyloscanner](https://github.com/BDI-pathogens/phyloscanner "phyloscanner"). We recommend batching 5-100 samples per phyloscanner run to avoid excessive runtimes for large datasets, and using the GTR+F+R6 (FreeRate) substitution model. Example input files are provided in ExampleInputs.

_Installation using conda:_
1. Clone this repository and cd to its location: `git clone git@github.com:BDI-pathogens/HIV-phyloTSI.git; cd HIV-phyloTSI`
2. Create the conda environment: `conda env create -f hivphylotsi.yml`
3. Activate it: `conda activate hivphylosti`
    * Alternatively, you can manually install the requirements: python3, pandas v0.25.3 and scikit-learn v0.22.
4. Test run using example input data:
`./HIVPhyloTSI.py -d Model -p ExampleInputs/testpatstats.csv -m ExampleInputs/testmaf.csv -o test.csv`

A csv file `test.csv` is created with TSI estimates for data in ExampleInputs/. The input and output data are described in more detail below.

_Input:_

Two CSV files:
1. **PatStats:** patStats.csv from phyloscanner, where host.id refers to individual samples (do not merge multiple timepoints); and
2. **MAF:** a csv file of minor allele frequencies (MAF), where rows match the host.id values in your **PatStats** file, and columns are labelled by position relative to HXB2. MAF is calculated from the shiver BaseFreqs_WithHXB2.csv, or from any pileup of reads, at each HXB2 position as follows:

            MAF = 1-max(A+C+G+T)/(A+C+G+T)

where A, C, G and T are counts of each base (ignoring gaps and Ns).

_Output:_

A CSV file is produced that contains the aggregated values for gag, pol, gp120 and gp41 regions for MAF at codon position 1/2 and codon position 3, largest.subgraph.rtt from phyloscanner (LRTT), number of tips and whether the sample looks like a dual infection. A binary indicator for amplicon-based sequencing is also added (is_mrc). The final columns are predictions from the model, as follows:

Square root TSI values (intermediate output):
*  RF_pred_sqrt - prediction of square-root transformed TSI
*  RF_std - standard error on the above, calculated as standard deviation of the full range of predictions from the 1000 decision trees in the forest
*  RF_cc025 and RF_cc975 - the 95% limits of the above distribution, ie the 2.5th and 97.5th percentile of predictions from the 1000 decision trees
*  RF_pred_MAE - the prediction error, calculated as the mean in square-root-space of the prediction from the error model

Linearised TSI values (final output):
*  RF_pred_linear - linearised TSI estimate ( model prediction in years ) - calculated as RF_pred_sqrt^2
*  RF_pred_min_linear and RF_pred_max_linear - limits of the prediction interval in linear TSI (ie the predicted seroconversion interval), calculated as 2 standard deviations of the error model above and below the point estimate.


