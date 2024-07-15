import numpy as np
from enum import Enum
from math import exp, pow
from osgeo import gdal
import os
import pandas as pd
import argparse

class IntegrationMethod(Enum):
    MINIMUM = 0
    MEAN = 1

def cov_simi(v1, v2, type, std, mean):
        simi = -1
        if type == 0: # 0 for continous type
            simi = exp(-pow((v1 - v2), 2) * 0.5 / pow(std, 4) * (pow(std, 2) + pow(mean - v1, 2)))
            if simi > 1:
                simi = 1
            if simi < 0:
                simi = 0
        elif type == 1: # 1 for categorical type
            if v1 == v2:
                simi = 1
            else:
                simi = 0
        return simi

def loc_simi(covs1, covs2, types, means,stds,integration = IntegrationMethod.MINIMUM):
    simi_tmp = []
    for i in range(len(covs1)):
        simi_tmp.append(cov_simi(covs1[i], covs2[i], types[i], stds[i], means[i]))
    if integration is IntegrationMethod.MEAN:
        return np.mean(simi_tmp)
    else:
        return np.min(simi_tmp)

def site_inference(covs,sample_covs,target_vals,means,stds,types,threshold = 0.):
    weight_sum = 0
    value_sum = 0
    uncertainty = 1
    for i in range(len(target_vals)):
        simi=loc_simi(sample_covs[i], covs,types,means,stds)
        uncertainty = 9999
        if simi > threshold:
            if simi>0 and 1-simi<uncertainty:
                uncertainty = 1-simi
            weight_sum = weight_sum + simi
            value_sum = value_sum + target_vals[i] * simi
    if weight_sum > 0:
        res = value_sum/weight_sum
    else:
        res = -1
    return res, uncertainty



def individualPredictiveMapping(covariateFilenames, sampleFile, target,outputTifName,uncerTifName = None,types=None):
    # sample file: csv, columns are: x,y,target
    #check data
    train = pd.read_csv(sampleFile)
    if 'X' in train.columns:
        col_x = 'X'
    elif 'x' in train.columns:
        col_x = 'x'
    else:
        return "sample file error"
    if 'Y' in train.columns:
        col_y = 'Y'
    elif 'y' in train.columns:
        col_y = 'y'
    else:
        return "sample file error"
    if target not in train.columns:
        return "sample file error"

    covs_fn = covariateFilenames.split('#')
    covs_num = len(covs_fn)
    covs_ar=[]
    for cov_fn in covs_fn:
        dataset = gdal.Open(cov_fn)
        band = dataset.GetRasterBand(1)
        cov_ar = band.ReadAsArray()
        na = band.GetNoDataValue()
        cov_ar[cov_ar==na]=np.nan
        covs_ar.append(cov_ar)
    covs_ar = np.array(covs_ar)
    if types is None:
        types = np.zeros(covs_num)
    if len(types) != covs_num:
        return "data types error"
    for i in range(covs_num):
        try:
            types[i] = int(types[i])
        except ValueError:
            return "data types error"
        if types[i] > 2 or types[i] < 0:
            return "data types error"
    mean = []
    std = []
    for i in range(covs_num):
        mean.append(np.nanmean(covs_ar[i]))
        std.append(np.nanstd(covs_ar[i]))

    # sample covs
    ysize, xsize = cov_ar.shape
    #get xy coordinates
    dataset = gdal.Open(covs_fn[0])
    geo = dataset.GetGeoTransform()
    xmin = geo[0]
    xgridsize = geo[1]
    ymin = geo[3]
    ygridsize = geo[5]
    sample_covs = []
    target_vals = []
    for i in range(len(train)):
        coordx = train[col_x][i]
        coordy = train[col_y][i]
        col_num = int((coordx-xmin)/xgridsize)
        row_num = int((coordy-ymin)/ygridsize)
        covs = covs_ar[:,row_num,col_num]
        if not np.isnan(covs).any():
            sample_covs.append(covs)
            target_vals.append(train[target][i])
        
    map_predict=np.zeros(covs_ar[0].shape)
    map_uncertainty=np.zeros(covs_ar[0].shape)
    NoData = -9999
    for i in range(covs_ar[0].shape[0]):
        for j in range(covs_ar[0].shape[1]):
            if np.isnan(covs_ar[:,i,j]).any():
                map_predict[i,j]=NoData
                map_uncertainty[i,j]=NoData
            else:
                p,u = site_inference(covs_ar[:,i,j],sample_covs,target_vals,mean,std,types)
                map_predict[i,j]=p
                map_uncertainty[i,j]=u

    map_predict_mask = map_predict
    map_predict_mask[map_predict_mask==-9999]=np.nan
    stat_min = np.nanmin(map_predict_mask)
    stat_max = np.nanmax(map_predict_mask)
    stat_mean = np.nanmean(map_predict_mask)
    stat_stdDev = np.nanstd(map_predict_mask)

    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.CreateCopy(outputTifName, dataset, strict=0)
    dst_ds.GetRasterBand(1).WriteArray(map_predict)
    dst_ds.GetRasterBand(1).SetStatistics(stat_min,stat_max,stat_mean,stat_stdDev)
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.GetRasterBand(1).ComputeBandStats(0)
    dst_ds = None
    
    if uncerTifName is not None:
        map_ucner_mask = map_uncertainty
        map_ucner_mask[map_ucner_mask==-9999]=np.nan
        stat_min = np.nanmin(map_ucner_mask)
        stat_max = np.nanmax(map_ucner_mask)
        stat_mean = np.nanmean(map_ucner_mask)
        stat_stdDev = np.nanstd(map_ucner_mask)

        driver = gdal.GetDriverByName('GTiff')
        dst_ds = driver.CreateCopy(uncerTifName, dataset, strict=0)
        dst_ds.GetRasterBand(1).WriteArray(map_uncertainty)
        dst_ds.GetRasterBand(1).SetStatistics(stat_min,stat_max,stat_mean,stat_stdDev)
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.GetRasterBand(1).ComputeBandStats(0)
        dst_ds = None
    

# ======function call example=========
#individualPredictiveMapping("C:\\work\\ipsm_neighbor_idw\\heshan\data\\env&mask_30m\\slope_30m.tif#C:\\work\\ipsm_neighbor_idw\\heshan\data\\env&mask_30m\\plan_30m.tif", "C:\\work\\ipsm_neighbor_idw\\heshan\data\\samples\\samples_purposive.csv", "C:\\work\\ipsm_neighbor_idw\\heshan\data\\python_test.tif","org")

parser = argparse.ArgumentParser(description='Individual predictive soil mapping')
parser.add_argument('-inlayers', dest = 'covariateFilenames', type=str, required=True,
                    help='filenames for covariates connected with #')
parser.add_argument('-sample', dest='sampleFile', type=str,required=True,
                    help='sample filename')
parser.add_argument('-predmap', dest='outputTifName', type=str,required=True,
                    help='predicted map filename')
parser.add_argument('-uncmap', dest='uncerTifName', type=str,required=False,
                    help='uncertainty map filename')
parser.add_argument('-target', dest='target', type=str,required=True,
                    help='field name for the target to be predicted in the sample file')
parser.add_argument('-datatypes', dest='datatypes', type=str,required=False,
                    help='datatypes for covariates, CONTINUOUS or CATEGORICAL, connected with #')
args = parser.parse_args()
#====parse para: datatype========
types = args.datatypes.split('#')
datatypes = []
for datatype in types:
    if datatype.upper()=="CONTINUOUS":
        datatypes.append(0)
    elif datatype.upper()=="CATEGORICAL":
        datatypes.append(1)
        
individualPredictiveMapping(args.covariateFilenames,args.sampleFile,args.target,args.outputTifName,args.uncerTifName,datatypes)

#file execution command example
#python "C:\\ipsm_demo\\individualPredictiveMapping.py" -inlayers "C:\\ipsm_demo\\slope_30m.tif#C:\\ipsm_demo\\plan_30m.tif" -sample C:\\ipsm_demo\\samples.csv -predmap C:\\ipsm_demo\\python_test.tif -target org -datatypes CONTINUOUS#CONTINUOUS