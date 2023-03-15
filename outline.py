import numpy as np
import matplotlib.pyplot as plt
from NumericalMethods.NumericalMetods import NumericalMethods
from ReactorConstants import ReactorConstants
from PIL import Image

if(__name__ == '__main__'):
    ##########################################
    ##               Step Zero              ##
    ##########################################
    
    #Recycle Feed Rate + Composition Initial Conditions
    converged = False
    RecycleFeed = 0
    Recycle_zTHF = 0
    Recycle_zMEOH = 0
    
    #VLE Data & Other Constants
    AzeotropeLine = ReactorConstants.azeotrope
    oatm_VLE = ReactorConstants.VLE_1atm
    tatm_VLE = ReactorConstants.VLE_10atm
    oatm_xb = ReactorConstants.xb_col1
    tatm_xb = ReactorConstants.xb_col2
    Feed = ReactorConstants.feedRate
    Feed_zTHF = ReactorConstants.zf
    Feed_zMEOH = 1 - ReactorConstants.zf
    oatm_Reflux = ReactorConstants.refluxRatio_col1
    oatm_xd = 0.454 #temp, upper bound on first col
    tatm_xd = 0.295 #temp, lower bound on second col
    
    
    
    oatm_Reflux = 1.5
    oatm_q = 1.2
    tatm_Reflux = 1.5
    tatm_q = 1.2
    
    
    iteration = 1
    while(True):
        ##########################################
        ##               Step One               ##
        ##########################################
        
        oatm_zTHF = ((Feed * Feed_zTHF) + (RecycleFeed * Recycle_zTHF)) / (Feed + RecycleFeed) #composition of THF going into first col
        oatm_zMEOH = ((Feed * Feed_zMEOH) + (RecycleFeed * Recycle_zMEOH)) / (Feed + RecycleFeed) #composition of MEOH going into first col
        
        
        ##########################################
        ##               Step Two               ##
        ##########################################
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
        
        ##########################################
        ##               Step Three             ##
        ##########################################
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
        
        ##########################################
        ##               Step Four              ##
        ##########################################
        
        if(tatm_dFeed - RecycleFeed < 0.5 and tatm_dFeed - RecycleFeed > -0.5):
            converged = True
            break
        else:
            RecycleFeed = tatm_dFeed
            print(f"Iteration {iteration}, Feed is {RecycleFeed}")
            iteration += 1
            Recycle_zTHF, Recycle_zMEOH = tatm_dComposition
    
    ##########################################
    ##               Step Five              ##
    ##########################################
    
    fig, ax = plt.subplots(1,2,figsize=(14, 8))
    fig.subplots_adjust(bottom=0.3)
    
    xvals_1atm = np.linspace(0,1) #az_1atm
    ax[0].set_title('1atm Column')
    ax[0].plot(xvals_1atm, np.polyval(AzeotropeLine, xvals_1atm), label='Azeotrope Line') #y=x line
    ax[0].plot(xvals_1atm,np.polyval(oatm_VLE, xvals_1atm), label='VLE Line') #VLE line
    ax[0].plot([oatm_zTHF, oatm_xPinch],[oatm_zTHF, oatm_yPinch], label='q Line') #qline
    ax[0].plot([oatm_xd, oatm_xPinch],[oatm_xd, oatm_yPinch], label='Operating Line') #operatingLine
    ax[0].plot([oatm_xb, oatm_xPinch],[oatm_xb, oatm_yPinch], label='Stripping Line') #strippingLine
    ax[0].legend()
    i=0
    while(i<len(oatm_data)-1):
        ax[0].plot([oatm_data[i][0], oatm_data[i+1][0]],[oatm_data[i][1],oatm_data[i+1][1]],linestyle='dashed', c='black')
        i+=1
        
    
    xvals_10atm = np.linspace(0, 1) #az_10atm
    ax[1].set_title('10atm Column')
    ax[1].plot(xvals_10atm, np.polyval(AzeotropeLine, xvals_10atm), label='Azeotrope Line') #y=x line
    ax[1].plot(xvals_10atm,np.polyval(tatm_VLE, xvals_10atm), label='VLE Line') #VLE line
    ax[1].plot([tatm_zTHF, tatm_xPinch],[tatm_zTHF, tatm_yPinch], label='q Line') #qline
    ax[1].plot([tatm_xd, tatm_xPinch],[tatm_xd, tatm_yPinch], label='Operating Line') #operatingLine
    ax[1].plot([tatm_xb, tatm_xPinch],[tatm_xb, tatm_yPinch], label='Stripping Line') #strippingLine
    ax[1].legend()
    i=0
    while(i<len(tatm_data)-1):
        ax[1].plot([tatm_data[i][0], tatm_data[i+1][0]],[tatm_data[i][1],tatm_data[i+1][1]],linestyle='dashed', c='black')
        i+=1
    
    print()
    print(f"Col1Distilate : xTHF = {oatm_dComposition[0]}, xMEOH = {oatm_dComposition[1]}, Flow = {oatm_dFeed}")
    print(f"Col1Bottoms : xTHF = {oatm_bComposition[0]}, xMEOH = {oatm_bComposition[1]}, Flow = {oatm_bFeed}")
    print(f"Col2Distilate : xTHF = {tatm_dComposition[0]}, xMEOH = {tatm_dComposition[1]}, Flow = {tatm_dFeed}")
    print(f"Col2Bottoms : xTHF = {tatm_bComposition[0]}, xMEOH = {tatm_bComposition[1]}, Flow = {tatm_bFeed}")
    
    plt.figtext(0.13, 0.22, "Distillate", fontsize='x-large', fontweight='bold')
    plt.figtext(0.13, 0.19, f"Flow = " + str(oatm_dFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.13, 0.16, f"xTHF = " + str(oatm_dComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.13, 0.13, f"xTHF = " + str(1-oatm_dComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.13, 0.10, "Reflux Ratio = ")
    plt.figtext(0.13, 0.07, "Condenser Duties = [kW]")
    
    plt.figtext(0.3, 0.22, "Bottoms", fontsize='x-large', fontweight='bold')
    plt.figtext(0.3, 0.19, f"Flow = " + str(oatm_bFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.3, 0.16, f"xTHF = " + str(oatm_bComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.3, 0.13, f"xTHF = " + str(1-oatm_bComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.3, 0.10, "Boilup Ratio = ")
    plt.figtext(0.3, 0.07, "Evaporator Duties = [kW]")
    
    plt.figtext(0.55, 0.22, "Recycle", fontsize='x-large', fontweight='bold')
    plt.figtext(0.55, 0.19, f"Flow = " + str(tatm_dFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.55, 0.16, f"xTHF = " + str(tatm_dComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.55, 0.13, f"xTHF = " + str(1-tatm_dComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.55, 0.10, "Reflux Ratio = ")
    plt.figtext(0.55, 0.07, "Condenser Duties = [kW]")
    
    plt.figtext(0.72, 0.22, "Bottoms", fontsize='x-large', fontweight='bold')
    plt.figtext(0.72, 0.19, f"Flow = " + str(tatm_bFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.72, 0.16, f"xTHF = " + str(tatm_bComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.72, 0.13, f"xTHF = " + str(1-tatm_bComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.72, 0.10, "Boilup Ratio = ")
    plt.figtext(0.72, 0.07, "Evaporator Duties = [kW]")
    
    
    img = np.array(Image.open('img/MPLDiagram.png'))
    fig2, ax2 = plt.subplots(1,1,figsize=(14, 8))
    ax2.axis('off')
    
    imgplot = plt.imshow(img)
    plt.show()