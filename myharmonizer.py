import json

from sklearn.preprocessing import MinMaxScaler
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics.pairwise import cosine_distances

import numpy as np
import pandas as pd

import tensorflow as tf
from tensorflow.keras import backend

from itertools import combinations
from scipy.spatial.distance import cityblock
from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

from pathlib import Path
import os
import re

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

class GlobalMinMaxScaler(BaseEstimator, TransformerMixin):
    """Min max scaling based on global min and max values"""
    def __init__(self, min_=None, max_=None):
        self.min_ = min_
        self.max_ = max_

    def fit(self, X, y=None):
        self.min_ = X.values.min()
        self.max_ = X.values.max()
        return self

    def transform(self, X, y=None):
        return (X - self.min_) / (self.max_ - self.min_)


def gene_length_scaling(datas, genelength):
    """Convenience function to scale features according to gene length"""

    # remove genes not in gene_length list
    allgenelength = datas.columns.isin(genelength.index)
    datas = datas.loc[:, allgenelength]
    genelength = genelength.loc[datas.columns]

    datas_scaled = datas.divide(genelength, axis='columns')
    return datas_scaled


class LSPreprocessing(BaseEstimator, TransformerMixin):
    """Sklearn-style class to perform library scale preprocessing a la Scanpy"""

    def fit(self, X, y=None):
        self.median_ = X.sum(axis=1).median()
        return (self)

    def transform(self, X, y=None):
        Xt = X.mul(self.median_ / X.sum(axis=1), axis="rows")
        return Xt.apply(np.log)


class TPMPreprocessing(BaseEstimator, TransformerMixin):
    """Class to calculate TPM from pseudocounts, gene length"""
    def __init__(self, genelength=None):
        self.genelength = genelength

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        rpk = gene_length_scaling(X, self.genelength)
        pkm = (rpk / 1E6).sum(axis=1)
        rpkm = rpk.div(pkm, axis=0)
        return rpkm.apply(np.log)


class QuantilePreprocessing(BaseEstimator, TransformerMixin):
    """Class to quantile preprocess dataset"""
    def __init__(self, means_=None):
        self.means_ = means_

    def fit(self, X, y=None):
        Xs = pd.DataFrame(np.sort(X.values, axis=1)).apply(np.log)
        means_ = Xs.mean(axis=0)
        means_.index = np.arange(1, len(means_) + 1)
        self.means_ = means_
        return self

    def transform(self, X, y=None):
        sorted = X.rank(method="min", axis=1).stack().astype(int).map(self.means_).unstack()
        return sorted


class RLEPreprocessing(BaseEstimator, TransformerMixin):
    """Class for relative log expression preprocessing, also called median of ratios method"""
    def __init__(self, pseudoref=None):
        self.pseudoref = pseudoref

    def fit(self, X, y=None):
        self.pseudoref = stats.gmean(X, axis=0)
        return self

    def transform(self, X, y=None):
        ls = X.div(self.pseudoref, axis=1).median(axis=1)
        return X.div(ls, axis=0).apply(np.log)


def rle_preprocessing(datas, geo_mean=None, floor=0.05, ceiling=0.95):
    """Modified version of rle scaling in DESeq2 - includes addition of 0.1 pseudocount since intended for large
    datasets - gene geometric means were predominated by 0s since at least one sample had a count of 0 for a given
    gene. Quantile capping of floor and ceiling also performed to mitigate effect of outliers."""

    datas, floors, ceilings = quantile_floorcap(datas + 0.1, floor=floor, ceiling=ceiling)

    if geo_mean is None:
        geo_mean = stats.gmean(datas, axis=0)
    # remove genes with geo_mean==0
    datas_z = datas.iloc[:, geo_mean > 0]
    geo_mean_z = geo_mean[geo_mean > 0]

    # calculate ratio of each sample to the reference
    datas_ratio = datas_z.div(geo_mean_z, axis="columns")

    # calculate median of sample-wise ratios
    datas_median = datas_ratio.median(axis=1)

    # scale samples
    datas_scaled = datas_z.div(datas_median, axis="index")

    # apply ln
    datas_ln = datas_scaled.apply(np.log)

    return datas_ln, geo_mean, floors, ceilings


def quantile_floorcap(datas, floor=0.05, ceiling=0.95):
    """Cap dataset at floor and ceiling quantiles to mitigate outlier effect"""
    if type(floor) is float and type(ceiling) is float:
        floors = datas.quantile(q=floor, axis="rows")
        ceilings = datas.quantile(q=ceiling, axis="rows")

    elif type(floor) is pd.Series and type(ceiling) is pd.Series:
        floors = floor
        ceilings = ceiling

    else:
        raise TypeError('Mismatch between floor and ceiling types or type is not a float/pandas Series.')

    for i, f in enumerate(floors):
        datas.iloc[datas.iloc[:, i] < f, i] = f

    for i, c in enumerate(ceilings):
        datas.iloc[datas.iloc[:, i] > c, i] = c

    return datas, floors, ceilings


def json_to_sklearn(data, modeltype):
    # with open(file, "r") as handle:
    #     data = json.load(handle)

    # Restore data types to init_params
    for name, p in data['init_params'].items():
        if data['init_params_types'][name] == 'ndarray':
            p = np.array(p)
        if data['init_params_types'][name] == 'Series':
            p = pd.Series(p, index=data['init_params_index'][name])

    match modeltype:
        case 'feature':
            model = MinMaxScaler(**data['init_params'])
        case 'global':
            model = GlobalMinMaxScaler(**data['init_params'])
        case 'LS':
            model = LSPreprocessing(**data['init_params'])
        case 'TPM':
            model = TPMPreprocessing(**data['init_params'])
        case 'QT':
            model = QuantilePreprocessing(**data['init_params'])
        case 'RLE':
            model = RLEPreprocessing(**data['init_params'])

    for name, p in data['model_params'].items():
        typep = data['model_params_types'][name]
        if typep == 'ndarray':
            setattr(model, name, np.array(p))
        elif typep == 'Series':
            # s = pd.Series()
            setattr(model, name,
                    pd.Series(p, index=data['model_params_index'][name]))
        else:
            setattr(model, name, p)
    return model


def chain(start, *funcs):
    res = start
    for func in funcs:
        res = func(res)
    return res

# Supporting functions and classes
class CenteredGaussianNoise(tf.keras.layers.Layer):
    """Apply additive non zero-centered Gaussian noise.

    As a regularization layer, it is only activate at training time."""

    def __init__(self, center, stddev, **kwargs):
        super().__init__(**kwargs)
        self.supports_masking = True
        self.stddev = stddev
        self.center = center

    def call(self, inputs, training=None):
        def noised():
            return inputs + backend.random_normal(
                shape=tf.shape(inputs),
                mean=self.center,
                stddev=self.stddev,
            )

        return backend.in_train_phase(noised, inputs, training=training)

    def get_config(self):
        config = {'center': self.center,
                  'stddev': self.stddev}
        base_config = super().get_config()
        return dict(list(base_config.items()) + list(config.items()))

def calculate_ccc(Xtest):
    df = Xtest.transpose()

    numerator = df.cov() * 2

    denominator = pd.DataFrame(np.zeros(numerator.shape), index=numerator.index, columns=numerator.columns)
    var = df.var()
    mu = df.mean()

    denominator = denominator.add(mu, axis='rows')
    denominator = (denominator.sub(mu, axis='columns')) ** 2

    denominator = denominator.add(var)
    denominator = denominator.add(var, axis='rows')

    ccc = numerator.div(denominator)
    return ccc

def calculate_euc(Xtest):
    # comb = pd.DataFrame(combinations(Xtest.index, 2))
    dis_m = np.zeros((Xtest.shape[0], Xtest.shape[0]))

    for i in range(Xtest.shape[0]):
        for j in range(Xtest.shape[0]):
            dis = np.linalg.norm(Xtest.iloc[i] - Xtest.iloc[j])
            dis_m[i, j] = dis
            dis_m[j, i] = dis

    return pd.DataFrame(dis_m, index=Xtest.index, columns=Xtest.index)

def calculate_manhattan(Xtest):
    comb = pd.DataFrame(combinations(Xtest.index, 2))
    dis_m = np.zeros((Xtest.shape[0], Xtest.shape[0]))

    for i in range(Xtest.shape[0]):
        for j in range(Xtest.shape[0]):
            dis = cityblock(Xtest.iloc[i], Xtest.iloc[j])
            dis_m[i, j] = dis
            dis_m[j, i] = dis

    return pd.DataFrame(dis_m, index=Xtest.index, columns=Xtest.index)




# Main class
class myHarmonizer:
    def __init__(self, myHarmonizerPath):

        with open(myHarmonizerPath, "r") as handle:
            harmonizer = json.load(handle)

        self.myHarmonizer_version = harmonizer['myHarmonizer_version']
        self.metadata = harmonizer['metadata']
        self.modelmeta = harmonizer['modelmeta']
        self.models = harmonizer['models']
        self.data = harmonizer['data']

    def niceify(self, data):
        """Clean a new dataset to bring it to a comparable representation as the raw data in the knowledge
        base corpus. Input data should be formatted as a pandas dataframe with samples as rows and features as columns.

        Features that cannot be mapped to the original representation will be excluded and features that are not present
        will be imputed with the median feature value from the base corpus before scaling and normalization.

        :param data: The input dataset
        :type data: :class:'pandas.DataFrame'
        """

        if not type(data) == pd.DataFrame:
            raise TypeError("Input data is not formatted as a pandas dataframe.")

        ## Remove genes not in feature list
        raw_medians = pd.Series(json.loads(self.models['raw_medians']))
        remove = data.columns.isin(raw_medians.index)
        data = data.loc[:, remove]
        print("{} gene features were excluded because they could not be mapped to the knowledge base.".format(
            sum(~remove)))

        ## Impute features not in feature list as well as NaN values
        include = ~raw_medians.index.isin(data.columns)
        dataz = data.reindex(columns=raw_medians.index)

        ## Add a pseudocount of 1 to all samples and cap at 99th percentile (sample-wise)
        dataz = dataz + 1
        scap = dataz.quantile(0.99, axis=1).apply(np.ceil)
        dataz = (dataz).clip(0, scap, axis=0).astype(int)



        for c in dataz.columns:
            dataz.loc[dataz[c].isnull(), c] = raw_medians[c]
        dataz = dataz.astype(int)

        print("{} gene features were imputed from the median value of the knowledge base.".format(sum(include)))

        return dataz

    def normalize(self, data):
        """Normalize a clean dataset to bring it to a comparable representation to the normalized data in the knowledge
        base corpus. Input data should be formatted as a pandas dataframe with samples as rows and features as columns
        and should already have passed through the niceify method.

        :param data: The input dataset
        :type data: :class:'pandas.DataFrame'
        """

        if not type(data) == pd.DataFrame:
            raise TypeError("Input data is not formatted as a pandas dataframe.")

        raw_medians = pd.Series(json.loads(self.models['raw_medians']))
        if not raw_medians.index.equals(data.columns):
            raise IndexError(
                "Expected gene features are not found in dataset. Try the niceify method before this method.")

        ## Load normalization model
        prep_method = self.modelmeta['normalize_scale']['preprocessing_method'][0]
        if prep_method in ['LS', 'QT', 'TPM', 'RLE']:
            normalize = json_to_sklearn(self.models['preprocessing'],
                                        prep_method)

            return normalize.transform(data)

        elif prep_method in ['VST', 'GeVST', 'TMM', 'GeTMM']:
            r = robjects.r

            # write temp .rds file
            poi = (Path.cwd() / 'temp.rds').as_posix()
            with open(poi, "w") as handle:
                handle.write(self.models['preprocessing'])

            # cast pandas df as r df
            with localconverter(robjects.default_converter + pandas2ri.converter):
                r_df = robjects.conversion.py2rpy(data)

                r_df.colnames = r.sub("-", "\\.", r.colnames(r_df))

            match prep_method:
                case 'VST':
                    # Read dispersion list
                    dispersionList = r.readRDS(poi)

                    # Get functions
                    r.source('VST_preprocessing_local.R')

                    # Run vst
                    r_prep = r.vst_preprocessing(r_df, geneCorr='none', dispersionList=dispersionList)

                    # Convert back to python object
                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        prep = robjects.conversion.rpy2py(r_prep)

                    # Remove temporary rds file
                    os.remove(poi)

                    return pd.DataFrame(prep, index=data.index, columns=r.colnames(r_prep))

                case 'GeVST':
                    dispersionList = r.readRDS(poi)

                    r.source('VST_preprocessing_local.R')
                    r_prep = r.vst_preprocessing(r_df, geneCorr='rpk', dispersionList=dispersionList)

                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        prep = robjects.conversion.rpy2py(r_prep)

                    # Remove temporary rds file
                    os.remove(poi)

                    return pd.DataFrame(prep, index=data.index, columns=r.colnames(r_prep))

                case 'TMM':
                    train_params = r.readRDS(poi)

                    r.source('TMM_preprocessing_local.R')
                    r_prep = r.tmm_preprocessing(r_df, referenceSample=train_params)

                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        prep = robjects.conversion.rpy2py(r_prep)

                    # Remove temporary rds file
                    os.remove(poi)

                    return pd.DataFrame(prep, index=data.index, columns=r.colnames(r_prep))

                case 'GeTMM':
                    train_params = r.readRDS(poi)

                    r.source('TMM_preprocessing_local.R')
                    r_prep = r.g_tmm_preprocessing(r_df, referenceSample=train_params)

                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        prep = robjects.conversion.rpy2py(r_prep)

                    # Remove temporary rds file
                    os.remove(poi)

                    return pd.DataFrame(prep, index=data.index, columns=r.colnames(r_prep))

        elif prep_method == 'None':
            return data

        else:
            raise ValueError(prep_method + ' preprocessing method from myHarmonizer object is unknown.')

    def scale(self, data):
        """Scale a clean dataset to bring it to a comparable representation to the scaled data in the knowledge
        base corpus. Input data should be formatted as a pandas dataframe with samples as rows and features as columns
        and should already have passed through the niceify method and the normalize method (opt).

        :param data: The input dataset
        :type data: :class:'pandas.DataFrame'
        """

        if not type(data) == pd.DataFrame:
            raise TypeError("Input data is not formatted as a pandas dataframe.")

        # raw_medians = pd.Series(json.loads(self.models['raw_medians']))
        # # When R imports data, it replaces - with . leading to changes in the feature names
        # datacolumns_r = pd.Index([re.sub(r"-", r".", d) for d in raw_medians.index])
        #
        # if (not raw_medians.index.equals(data.columns)) and (not datacolumns_r.equals(data.columns)):
        #     raise IndexError(
        #         "Expected gene features are not found in dataset. Try the niceify method before this method.")

        ## Load normalization model
        scaling = json_to_sklearn(self.models['scaling'],
                                  self.modelmeta['normalize_scale']['scaling_method'][0].lower())

        return pd.DataFrame(scaling.transform(data),
                            index=data.index,
                            columns=data.columns)
    def encode(self, data):
        """Encode a cleaned and scaled dataset to bring it to a comparable representation to the encoded data
         in the knowledge base corpus. Input data should be formatted as a pandas dataframe with samples as rows and
         features as columns and should have already passed through the niceify, normalize, and scale methods.

         :param data: The input dataset
         :type data: :class:'pandas.DataFrame'
         """
        if not type(data) == pd.DataFrame:
            raise TypeError("Input data is not formatted as a pandas dataframe.")

        # raw_medians = pd.Series(json.loads(self.models['raw_medians']))
        # # When R imports data, it replaces - with . leading to changes in the feature names
        # datacolumns_r = pd.Index([re.sub(r"-", r".", d) for d in raw_medians.index])
        #
        # if (not raw_medians.index.equals(data.columns)) and (not datacolumns_r.equals(data.columns)):
        #     raise IndexError(
        #         "Expected gene features are not found in dataset. Try the niceify method before this method.")

        ## Load encoder
        architecture = self.models['encoder_model_architecture']
        weights_unlisted = self.models['encoder_model_weights']
        weights = [np.array(w) for w in weights_unlisted]

        encoder = tf.keras.models.model_from_json(architecture,
                                                  custom_objects={'CenteredGaussianNoise': CenteredGaussianNoise})
        encoder.set_weights(weights)

        # Subset features in encoder input
        dataz = data[self.modelmeta['encoder_metrics']['features']]

        ## Get data representation
        rep = pd.DataFrame(encoder.predict(dataz),
                            index=data.index)
        rep.columns = rep.columns.astype('str')

        return rep

    def transform(self, data):
        """Run full harmonization pipeline to bring a new dataset to a comparable representation to the encoded data
                 in the knowledge base corpus. Input data should be formatted as a pandas dataframe with samples as rows and
                 features as columns and should have already passed through the niceify, normalize, and scale methods.

                 :param data: The input dataset
                 :type data: :class:'pandas.DataFrame'
                 """

        return chain(data, self.niceify, self.normalize, self.scale, self.encode)

### Define main similarity function
def similarity(dataset1, dataset2, metric='Pearson'):
    """Calculate similarity between dataset1 and dataset2 based on metric. Input datasets should be
    pandas DataFrames.

    Pearson and Spearman metrics are correlation coefficients and output values from [-1,1] where 1 indicates the highest
    positive inter-dataset correlation and -1 indicates the highest negative correlation between datasets. Standard
    definitions and methods have been used to calculate these distances (see :func:`pandas.DataFrame.corr()`)

    CCC is the intraclass correlation coefficient, which is a measure of agreement. Similar to traditional correlation
    coefficients, it has a range of [-1,1] with similar interpretation. The principle difference is that this metric
    also considers differences in scale between datasets. For this implementation, we used the following definition:
    .. math::
        \\rho_{ccc} = \\rho * \\frac{2 * \\sigma^{(1)}\\sigma^{(2)}}{\\sigma^{(1)2} + \\sigma^{(2)2} + (\\mu^{(1)} - \\mu^{(2)})^2}


    Euclidean distance, also called Pythagorean distance, is implemented using the standard mathematical form:
    .. math::
        d = \\sqrt{\\sum(x^{(1)} - x^{(2)})^2}

    as is the Manhattan distance (see :func:`scipy.spatial.distance.cityblock()`)
    .. math::
        d_{manhattan} = \\sum |x^{(1)} - x^{(2)}
    Both distances have a range of [0, inf).

    Cosine distance also has a standard definition with a range of [0,1] where 0 would indicate no distance between
    datasets (see :func:`sklearn.metrics.pairwise.cosine_distances()`)
    .. math::
        d_{cos} = 1 - \\frac{\\sum x^{(1)}x^{(2)}}{\\sqrt{\\sum x^{(1)2}} \\sqrt{\\sum x^{(2)2}}}


    :param dataset1: The first dataset to be compared.
    :type dataset1: :class:'pandas.DataFrame'

    :param dataset2: The second dataset to be compared.
    :type dataset2: :class:'pandas.DataFrame'

    :param metric: One of 'Pearson', 'Spearman', 'CCC', 'Euclidean', 'Manhattan', or 'Cosine'.
    """

    if not type(dataset1) == pd.DataFrame or not type(dataset2) == pd.DataFrame:
        raise TypeError("Input data is not formatted as a pandas dataframe.")

    if not dataset1.columns.equals(dataset2.columns):
        raise IndexError(
            "Datasets contain different feature lists (columns).")

    data = pd.concat([dataset1, dataset2])

    match metric:
        case 'Pearson':
            sim = data.transpose().corr()

        case 'Spearman':
            sim = data.transpose().corr(method="spearman")

        case 'CCC':
            sim = calculate_ccc(data)

        case 'Euclidean':
            sim = calculate_euc(data)

        case 'Manhattan':
            sim = calculate_manhattan(data)

        case 'Cosine':
            sim = pd.DataFrame(cosine_distances(data), index=data.index, columns=data.index)

        case other:
            raise KeyError(metric + " is not one of Pearson, Spearman, CCC, Euclidean, Manhattan, or Cosine.")



    # Subset similarity so that dataset1 is along rows and dataset2 along columns

    sim = sim.loc[dataset1.index, dataset2.index]

    return sim

### Build heatmap

def listmetadata(myHarmonizer):
    """Returns the metadata keys available in the knowledge base

    :param myHarmonizer: myHarmonizer object
    :type myHarmonizer: :class:'myHarmonizer'
    """

    if myHarmonizer.metadata:
        meta = pd.DataFrame(json.loads(myHarmonizer.metadata))

        return meta.columns.to_list()

    else:
        print('No sample metadata found for this myHarmonizer object.')

def heatmap(similaritydf, myHarmonizer, user_metadata=None, kb_metadata=None):
    """Convenience function to draw heatmap with optional metadata.

    :param similaritydf: Similarity dataframe.
    :type similaritydf: :class:'pandas.DataFrame'

    :param user_metadata: A Series of discrete metadata with sample names in the Index.
    :type user_metadata: :class:'pandas.Series'

    :param kb_metadata: A string representing one sample metadata column in the knowledge base. Options can be found with the :func:`listmetadata()` function.
    :type kb_metadata: :class:`str`

    """

    plt.clf()
    row_colors = None
    col_colors = None

    if kb_metadata is not None:
        kbmeta = pd.DataFrame(json.loads(myHarmonizer.metadata))[kb_metadata]

        # Assign colors to each unique metadata category for kb
        kbuniquemeta = kbmeta.unique()
        palette1 = dict(zip(kbuniquemeta, sns.color_palette("husl", len(kbuniquemeta))))
        kbhandles = [Patch(facecolor=palette1[name]) for name in palette1]
        col_colors = kbmeta.map(palette1)


    if user_metadata is not None:
        # Assign colors to each unique metadata category for user
        uniquemeta = user_metadata.unique()
        palette2 = dict(zip(uniquemeta, sns.color_palette("husl", len(uniquemeta))))
        userhandles = [Patch(facecolor=palette2[name]) for name in palette2]
        row_colors = user_metadata.map(palette2)

    cg = sns.clustermap(similaritydf.round(2),
                        cmap="mako",
                        row_colors=row_colors,
                        col_colors=col_colors)

    if kb_metadata is not None:
        ax = cg.ax_heatmap
        first_legend = ax.legend(kbhandles, palette1, bbox_to_anchor=(0.2, 1),
                                 bbox_transform=plt.gcf().transFigure, loc="upper right")
        ax.add_artist(first_legend)

    if user_metadata is not None:
        ax = cg.ax_heatmap
        ax.legend(userhandles, palette2, bbox_to_anchor=(0.35, 1),
                  bbox_transform=plt.gcf().transFigure, loc="upper right")

    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=30)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=30)

    return plt




