#!/usr/bin/env python
__doc__ ='''


~~~~~~~~~~~~~
HIV-PhyloTSI:
Estimate time since infection (TSI) for HIV samples using deep sequenced short-read data. Both amplicon-based and capture-based protocols are supported.
~~~~~~~~~~~~~

Input:
Two CSV files:
    (1) patStats.csv from phyloscanner, where host.id refers to individual samples (do not merge multiple timepoints); and
    (2) corresponding CSV where rows are the same host.id values as in (1), and columns are positions, containing values of cumulative minor allele frequency (MAF), calculated from the shiver BaseFreqs.csv at each position as:

            MAF = (1-(A+C+G+T))/(A+C+G+T)

Output:
A CSV file is produced that contains the aggregated values for gag, pol, gp120 and gp41 regions for MAF at codon position 1/2 and codon position 3, largest.subgraph.rtt from phyloscanner (LRTT), number of tips and whether the sample looks like a dual infection. A binary indicator of amplicon-based sequencing is also added (is_mrc). The final columns are the predictions from the model, as follows:

Square root TSI values:
  RF_pred_sqrt - prediction of square-root transformed TSI
  RF_std - standard error on the above, calculated as standard deviation of the full range of predictions from the 1000 decision trees in the forest
  RF_cc025 and RF_cc975 - the 95% limits of the above distribution, ie the 2.5th and 97.5th percentile of predictions from the 1000 decision trees
  RF_pred_MAE - the prediction error, calculated as the mean in square-root-space of the prediction from the error model
Linearised TSI values:
  RF_pred_linear - linearised TSI estimate ( model prediction in years ) - calculated simply as RF_pred_sqrt^2
  RF_pred_min_linear and RF_pred_max_linear - limits of the prediction interval in linear TSI (ie the predicted seroconversion interval), calculated as 2 standard deviations of the error model above and below the point estimate.

'''

import argparse, sys, os.path
import numpy as np
import pandas as pd
import pickle
import warnings
warnings.filterwarnings("ignore")

# ============================================================================ #
# GLOBAL VARIABLES                                                             #
# ============================================================================ #
_args = None
_progname=os.path.basename(sys.argv[0])

# ============================================================================ #
# LOGGING                                                                      #
# ============================================================================ #
def loginfo(s, verbose=True):
    if verbose:
        sys.stderr.write('  Info: {0}\n'.format(s))
def logerr(s):
    sys.stderr.write('  Warning: {0}\n'.format(s))
def stoperr(s, errcode=1):
    errword = 'Finished' if not errcode else 'Error'
    sys.stderr.write('  {0}: {1}\n'.format(errword, s))
    sys.exit(errcode)

# ============================================================================ #
# PROGRAM USAGE                                                                #
# ============================================================================ #
def Initialise():
    '''
    Parse command-line arguments.
    '''
    global _args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, epilog='Example: {progname} -p /path/to/patstats/file -m /path/to/maf/file -m /path/to/model/directory -o output.csv'.format(progname=_progname))
    parser.add_argument( '-p', '--patstatspath', required=True, help='PatStats.csv file from running phyloscanner on all samples of interest.' )
    parser.add_argument( '-m', '--mafpath', required=True, help='Minor allele frequency data with one row per sample and columns for each genomic position.' )
    parser.add_argument( '-o', '--outpath', required=True, help='Filename for saving predictions.' )
    parser.add_argument( '-d', '--modeldir', required=True, help='Path to directory containing models and reference data.')
    parser.add_argument( '--modelname', required=False, default='MEANS_feats_LRTT', help='Name of model (default: MEANS_feats_LRTT).' )
    parser.add_argument( '--amplicons', type=bool, required=False, default=False, help='Is the data amplicons? If in doubt, leave as False (default: False) ' )
    _args = parser.parse_args()
    return

def Clean_Commandline():
    '''
    Print errors, and exit if necessary, on bad input data.
    '''
    # Validate input file
    if not os.path.isfile(_args.patstatspath):
        stoperr('Input file {} does not exist or is not accessible'.format(_args.patstatspath))
    if not os.path.isfile(_args.mafpath):
        stoperr('Input file {} does not exist or is not accessible'.format(_args.mafpath))
    if not os.path.isdir(_args.modeldir):
        stoperr('Model directory {} does not exist or is not accessible'.format(_args.modeldir))
    return

# ============================================================================ #
# DATA PROCESSING FUNCTIONS                                                    #
# ============================================================================ #
def load_reference_data(modeldir):
    ''' Load HXB2 data.'''
    #  HXB2 positions
    hxb2 = pd.read_csv('{}/HXB2_refdata.csv'.format(modeldir))
    hxb2['position']=hxb2['HXB2 base position']
    hxb2.set_index('position' , inplace=True)
    # Third codon position sites
    rf3_3cp = hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].max()
    rf2_3cp = hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].max()
    rf1_3cp = hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].max()
    # Make into set
    t1 = set(rf1_3cp.reset_index()['HXB2 base position'])
    t2 = set(rf2_3cp.reset_index()['HXB2 base position'])
    t3 = set(rf3_3cp.reset_index()['HXB2 base position'])
    # Same for 1st/2nd position
    rf1_12 = set(hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    rf2_12 = set(hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    rf3_12 = set(hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    # Summarise
    first_second_codon_pos = set(rf1_12 | rf2_12 | rf3_12)
    third_codon_pos = set(t1 | t2 | t3) - first_second_codon_pos
    # Genes
    gp120 = hxb2[hxb2['RF3 protein']=='gp120'].index
    gp41 = hxb2[hxb2['RF3 protein']=='gp41'].index
    gag = hxb2[hxb2['RF1 protein']=='gag'].index
    pol = hxb2[hxb2['RF3 protein']=='pol'].index
    return first_second_codon_pos, third_codon_pos, gag, pol, gp120, gp41


def load_model(modelname, modeldir):
    ''' Load pre-fitted model and scaler. '''
    with open('{}/full_scaler.selfscale.pkl'.format(modeldir), 'rb') as inf:
        full_scaler = pickle.load(inf, fix_imports=True, encoding="latin1")
    with open('{}/full_model.selfscale.{}.pkl'.format(modeldir, modelname), 'rb') as inf:
        full_model_selfscale = pickle.load(inf, fix_imports=True, encoding="latin1")
    modelfeats = np.genfromtxt('{}/feature_list_{}.txt'.format(modeldir, modelname), dtype='str')
    with open('{}/err_model_mae.{}.pkl'.format(modeldir, modelname), 'rb') as inf:
        err_model = pickle.load(inf, fix_imports=True, encoding="latin1")
    return full_model_selfscale, full_scaler, modelfeats, err_model


def load_patstats(fpath):
    ''' Load phyloscanner output - PatStats.csv '''
    patstats = pd.read_csv(fpath)
    Xlrtt = patstats.groupby(['host.id', 'xcoord'])['normalised.largest.rtt'].mean().unstack()
    Xtips = patstats.groupby(['host.id', 'xcoord'])['tips'].mean().unstack()
    Xdual = patstats.groupby(['host.id', 'xcoord'])['solo.dual.count'].mean().unstack()
    try:
        assert (Xlrtt.index == Xtips.index).all()
        assert (Xlrtt.index == Xdual.index).all()
    except AssertionError:
        logerr('Index mismatch between phyloscanner outputs: {}, {}, {}').format(Xlrtt.shape, Xtips.shape, Xdual.shape)
    loginfo('Loaded phyloscanner data, shape={}'.format(Xlrtt.shape))
    return Xlrtt, Xtips, Xdual


def load_maf(fpath):
    ''' Load shiver output processed to give cumumative minor allele frequencies. '''
    Xmaf = pd.read_csv(fpath, index_col=0)
    Xmaf.columns = [int(float(c)) for c in Xmaf.columns]
    loginfo('Loaded MAF data, shape={}'.format(Xmaf.shape))
    return Xmaf


def generate_features(Xlrtt, Xtips, Xdual, Xmaf, modeldir, modelfeats):
    ''' Generate aggregated predictors. '''
    loginfo('Model Features: {}'.format(modelfeats))
    first_second_codon_pos, third_codon_pos, gag, pol, gp120, gp41 = load_reference_data(modeldir)
    Xgene_means = pd.DataFrame()
    Xgene_means['genome_lrtt'] = Xlrtt.mean(axis=1)
    Xgene_means['genome_maf12c'] = Xmaf.loc[:,list(first_second_codon_pos)].mean(axis=1)
    Xgene_means['genome_maf3c'] = Xmaf.loc[:,list(third_codon_pos)].mean(axis=1)
    Xgene_means['genome_tips'] = Xtips.mean(axis=1)
    Xgene_means['genome_dual'] = Xdual.mean(axis=1)
    for g in ['gag', 'pol', 'gp120', 'gp41']:
        Xgene_means['{}_lrtt'.format(g)] = Xlrtt.loc[:, eval(g)].mean(axis=1)
        Xgene_means['{}_tips'.format(g)] = Xtips.loc[:, eval(g)].mean(axis=1)
        Xgene_means['{}_dual'.format(g)] = Xdual.loc[:, eval(g)].mean(axis=1)
        Xgene_means['{}_maf12c'.format(g)] = Xmaf.loc[:, list( (set(eval(g)) & first_second_codon_pos) )].mean(axis=1)
        Xgene_means['{}_maf3c'.format(g)] = Xmaf.loc[:, list( (set(eval(g)) & third_codon_pos) )].mean(axis=1)
    return Xgene_means


def _impute_knn(X, k=3):
    ''' Use k nearest rows which have a feature to fill in each row's
    missing features. '''
    if np.isnan(X).sum() == 0:
        loginfo('No missing data to impute, continuing.')
        return X
    # Number of rows (samples) with missing data
    nm = (np.isnan(X).sum(axis=1)!=0).sum()
    loginfo('Using KNN (3 nearest neighbours) to impute {} missing values ({} sample{})'.format(np.isnan(X).sum(),
                                                    nm, ('','s')[nm>1]))
    from sklearn.impute import KNNImputer
    X_knn = KNNImputer(n_neighbors=k).fit_transform(X)
    return X_knn

def prepare_data(full_scaler, Xgene_means, modelfeats, is_amplicons=False):
    ''' Prepare features for model. '''
    seqtype = pd.DataFrame(index=Xgene_means.index)
    seqtype['is_mrc'] = int(is_amplicons)
    data_for_model = pd.concat((seqtype, Xgene_means), axis=1)[modelfeats]
    data_for_model_z = full_scaler.transform(data_for_model)
    data_for_model_imp = _impute_knn(data_for_model_z)
    return data_for_model_imp


def predict_from_model(Xgene_means, data_for_model_imp, full_model_selfscale, err_model):
    ''' Apply model to data and generate estimates and prediction intervals. '''
    loginfo('Collecting estimates...')
    preds = Xgene_means.copy()
    preds['RF_pred_sqrt'] = full_model_selfscale.predict(data_for_model_imp)
    cc_dist = np.asarray([est.predict(data_for_model_imp) for est in full_model_selfscale.estimators_])
    # Standard error (of mean estimate)
    preds['RF_std'] = cc_dist.std(axis=0)
    preds['RF_cc025'] = np.apply_along_axis(np.percentile, 0, cc_dist, q=2.5)
    preds['RF_cc975'] = np.apply_along_axis(np.percentile, 0, cc_dist, q=97.5)
    # Prediction error and 95% prediction interval in square-root space
    preds['RF_pred_MAE'] = err_model.predict(data_for_model_imp)
    # Linear predictions
    preds['RF_pred_linear'] = preds.RF_pred_sqrt**2
    preds['RF_pred_min_linear'] = (preds.RF_pred_sqrt - 1.96*preds.RF_pred_MAE).apply(lambda x: max(0, x)**2)
    preds['RF_pred_max_linear'] = (preds.RF_pred_sqrt + 1.96*preds.RF_pred_MAE).apply(lambda x: max(0, x)**2)
    loginfo('Done.')
    return preds


def save_predictions(preds, outf):
    ''' Write output. '''
    preds.to_csv(outf)
    loginfo('Output saved as {}.'.format(outf))

# ============================================================================ #
# MAIN                                                                         #
# ============================================================================ #
def main():
    ''' Run the pipeline. '''
    full_model_selfscale, full_scaler, \
                modelfeats, err_model = load_model(_args.modelname, _args.modeldir)
    Xlrtt, Xtips, Xdual = load_patstats(_args.patstatspath)
    Xmaf = load_maf(_args.mafpath)
    Xgene_means = generate_features(Xlrtt, Xtips, Xdual, Xmaf, _args.modeldir, modelfeats)
    data_for_model_imp = prepare_data(full_scaler, Xgene_means, modelfeats, _args.amplicons)
    preds = predict_from_model(Xgene_means, data_for_model_imp, full_model_selfscale, err_model)
    save_predictions(preds, _args.outpath)

if __name__ == '__main__':
    Initialise()
    Clean_Commandline()
    main()
# ============================================================================ #
