import numpy as np
import matplotlib.pyplot as plt
from NumericalMethods.NumericalMetods import NumericalMethods
from ReactorConstants import ReactorConstants
from PIL import Image
import pandas as pd

def datavsr(r):
    #Recycle Feed Rate + Composition Initial Conditions
    converged = False
    RecycleFeed = 0
    Recycle_zTHF = 0
    Recycle_zMEOH = 0
    
    #VLE Data & Other Constants
    AzeotropeLine = ReactorConstants.azeotrope
    oatm_VLE = ReactorConstants.VLE_1atm
    tatm_VLE = ReactorConstants.VLE_10atm
    oatm_xT = ReactorConstants.xT_1atm
    tatm_xT = ReactorConstants.xT_10atm
    oatm_xb = ReactorConstants.xb_col1
    tatm_xb = ReactorConstants.xb_col2
    Feed = ReactorConstants.feedRate
    Feed_zTHF = ReactorConstants.zf
    Feed_zMEOH = 1 - ReactorConstants.zf
    oatm_Reflux = ReactorConstants.refluxRatio_col1
    THF_VapHeat = ReactorConstants.THF_VapHeat
    MEOH_VapHeat = ReactorConstants.MEOH_VapHeat
    oatm_xd = 0.454 #temp, upper bound on first col
    tatm_xd = 0.295 #temp, lower bound on second col
    
    
    
    oatm_Reflux = r
    oatm_q = 1
    tatm_Reflux = r
    tatm_q = 1
    
    
    FeedArr = [0]
    while(True):
        
        oatm_zTHF = ((Feed * Feed_zTHF) + (RecycleFeed * Recycle_zTHF)) / (Feed + RecycleFeed) #composition of THF going into first col
        oatm_zMEOH = ((Feed * Feed_zMEOH) + (RecycleFeed * Recycle_zMEOH)) / (Feed + RecycleFeed) #composition of MEOH going into first col
        
        oatm_steps, oatm_data, oatm_lines = NumericalMethods.staircaseData1atm(oatm_q, oatm_zTHF, oatm_Reflux, oatm_VLE, oatm_xd, oatm_xb)
        oatm_qLine, oatm_oLine, oatm_sLine = oatm_lines
        oatm_xPinch, oatm_yPinch = NumericalMethods.getPolynomialIntercept(oatm_qLine, oatm_oLine)
        
        
        oatm_bComposition = [oatm_data[0][0], 1 - oatm_data[0][0]] #zTHF, zMEOH
        oatm_dComposition = [oatm_data[-1][0], 1 - oatm_data[-1][0]] #zTHF, zMEOH
        
        varM1 = [[1,1],[oatm_dComposition[0], oatm_bComposition[0]]] #variable Matrix
        ansM1 = [[Feed + RecycleFeed],[(Feed + RecycleFeed) * oatm_zTHF]] #answer Matrix
        
        oatm_dFeed, oatm_bFeed = np.linalg.solve(varM1, ansM1) #feed rate of top & bottom
        oatm_dFeed = oatm_dFeed[0]
        oatm_bFeed = oatm_bFeed[0]
        oatm_dTemperature = np.polyval(oatm_xT, oatm_dComposition[0])
        oatm_Boilup = -oatm_bComposition[0]/oatm_sLine[1]
        oatm_L = oatm_dFeed * oatm_Reflux
        oatm_VBar = oatm_bFeed * oatm_Boilup
        oatm_VapHeatCondenser = (oatm_dComposition[0]*THF_VapHeat)+(oatm_dComposition[1]*MEOH_VapHeat)
        oatm_VapHeatEvaporator = (oatm_bComposition[0]*THF_VapHeat)+(oatm_bComposition[1]*MEOH_VapHeat)
        oatm_Condenser = oatm_L * oatm_VapHeatCondenser * 0.0002777778 #convertion of kJ/kmol to kJ/second (aka kW)
        oatm_Evaporator = oatm_VBar * oatm_VapHeatEvaporator * 0.0002777778
        
        tatm_zTHF, tatm_zMEOH = oatm_dComposition
        
        tatm_steps, tatm_data, tatm_lines = NumericalMethods.staircaseData10atm(tatm_q, tatm_zTHF, tatm_Reflux, tatm_VLE, tatm_xd, tatm_xb)
        tatm_qLine, tatm_oLine, tatm_sLine = tatm_lines
        tatm_xPinch, tatm_yPinch = NumericalMethods.getPolynomialIntercept(tatm_qLine, tatm_oLine)
        
        tatm_bComposition = [tatm_data[0][0], 1 - tatm_data[0][0]] #zTHF, zMEOH
        tatm_dComposition = [tatm_data[-1][0], 1 - tatm_data[-1][0]] #zTHF, zMEOH
        
        varM10 = [[1,1],[tatm_dComposition[0], tatm_bComposition[0]]] #variable Matrix
        ansM10 = [[oatm_dFeed],[oatm_dFeed * tatm_zTHF]] #answer Matrix
        
        tatm_dFeed, tatm_bFeed = np.linalg.solve(varM10, ansM10) #feed rate of top & bottom, distilled = recycle
        tatm_dFeed = tatm_dFeed[0]
        tatm_bFeed = tatm_bFeed[0]
        tatm_dTemperature = np.polyval(tatm_xT, tatm_dComposition[0])
        tatm_Boilup = -tatm_bComposition[0]/tatm_sLine[1]
        tatm_L = tatm_dFeed * tatm_Reflux
        tatm_VBar = tatm_bFeed * tatm_Boilup
        tatm_VapHeatCondenser = (tatm_dComposition[0]*THF_VapHeat)+(tatm_dComposition[1]*MEOH_VapHeat)
        tatm_VapHeatEvaporator = (tatm_bComposition[0]*THF_VapHeat)+(tatm_bComposition[1]*MEOH_VapHeat)
        tatm_Condenser = tatm_L * tatm_VapHeatCondenser * 0.0002777778 #convertion of kJ/kmol to kJ/second (aka kW)
        tatm_Evaporator = tatm_VBar * tatm_VapHeatEvaporator * 0.0002777778
        
        if(tatm_dFeed - RecycleFeed < 0.5 and tatm_dFeed - RecycleFeed > -0.5):
            converged = True
            RecycleFeed = tatm_dFeed
            FeedArr.append(RecycleFeed)
            break
        else:
            RecycleFeed = tatm_dFeed
            FeedArr.append(RecycleFeed)
            Recycle_zTHF, Recycle_zMEOH = tatm_dComposition
        
    col1data = [r, oatm_Condenser, oatm_Evaporator, oatm_Boilup, oatm_steps]
    col10data = [r, tatm_Condenser, tatm_Evaporator, tatm_Boilup, tatm_steps]
    return col1data, col10data

col1DF = pd.DataFrame({
    "reflux" : [],
    "condenser" : [],
    "evaporator" : [],
    "boilup" : [],
    "steps" : []
})

col10DF = pd.DataFrame({
    "reflux" : [],
    "condenser" : [],
    "evaporator" : [],
    "boilup" : [],
    "steps" : []
})

refluxVals = range(1,6)
for reflux in refluxVals:
    col1data, col10data = datavsr(reflux)
    col1DF.loc[len(col1DF.index)] = col1data
    col10DF.loc[len(col10DF.index)] = col10data

col1DF.to_csv("Data/oatmData.csv", encoding='utf-8', index=False)
col10DF.to_csv("Data/tatmdata.csv", encoding='utf-8', index=False)
