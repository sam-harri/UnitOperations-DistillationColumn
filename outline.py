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
    Recycle_zTHF = 0.3
    Recycle_zMEOH = 0.7
    
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
    Feed_Temperature = ReactorConstants.Feed_Temperature
    THF_CPLiq = ReactorConstants.THF_CPLiq
    MEOH_CPLiq = ReactorConstants.MEOH_CPLiq
    
    oatm_Reflux = 1.35
    oatm_q = 1
    tatm_Reflux = 1.32
    tatm_q = 1
    
    
    FeedArr = [0]
    while(True):
        ##########################################
        ##               Step One               ##
        ##########################################
        oatm_zTHF = ((Feed * Feed_zTHF) + (RecycleFeed * Recycle_zTHF)) / (Feed + RecycleFeed) #composition of THF going into first col
        oatm_zMEOH = ((Feed * Feed_zMEOH) + (RecycleFeed * Recycle_zMEOH)) / (Feed + RecycleFeed) #composition of MEOH going into first col
        
        
        RecycleCP = (Recycle_zTHF*THF_CPLiq)+(Recycle_zMEOH*MEOH_CPLiq)
        FeedCP = (Feed_zTHF*THF_CPLiq)+(Feed_zMEOH*MEOH_CPLiq)
        oatm_CP = (oatm_zTHF*THF_CPLiq)+(oatm_zMEOH*MEOH_CPLiq)
        Recycle_Temperature = np.polyval(tatm_xT, Recycle_zTHF)
        oatm_PreTemperature = ((Feed * FeedCP * Feed_Temperature)+(RecycleFeed * RecycleCP * Recycle_Temperature))/((Feed * FeedCP)+(RecycleFeed * RecycleCP))
        oatm_WantedTemperature = np.polyval(oatm_xT, oatm_zTHF)
        oatm_TempertureDifference = oatm_PreTemperature - oatm_WantedTemperature
        oatm_EntryDuty = oatm_TempertureDifference * oatm_CP * (Feed + RecycleFeed) * 0.0002777778
        
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
        oatm_dTemperature = np.polyval(oatm_xT, oatm_dComposition[0])
        oatm_Boilup = -oatm_bComposition[0]/oatm_sLine[1]
        oatm_L = oatm_dFeed * oatm_Reflux
        oatm_VBar = oatm_bFeed * oatm_Boilup
        oatm_VapHeatCondenser = (oatm_dComposition[0]*THF_VapHeat)+(oatm_dComposition[1]*MEOH_VapHeat)
        oatm_VapHeatEvaporator = (oatm_bComposition[0]*THF_VapHeat)+(oatm_bComposition[1]*MEOH_VapHeat)
        oatm_Condenser = oatm_L * oatm_VapHeatCondenser * 0.0002777778 #convertion of kJ/kmol to kJ/second (aka kW)
        oatm_Evaporator = oatm_VBar * oatm_VapHeatEvaporator * 0.0002777778
        
        ##########################################
        ##               Step Three             ##
        ##########################################
        tatm_zTHF, tatm_zMEOH = oatm_dComposition
        tatm_Temperature = np.polyval(oatm_xT, tatm_zTHF)
        tatm_q1Temperature = np.polyval(tatm_xT, tatm_zTHF)
        tatm_VapHeatEntry = oatm_VapHeatCondenser
        tatm_CpEntry = (tatm_zTHF*THF_CPLiq)+(tatm_zMEOH*MEOH_CPLiq)
        tatm_newq = (tatm_VapHeatEntry+(tatm_CpEntry*(tatm_q1Temperature-tatm_Temperature)))/(tatm_VapHeatEntry)
        tatm_q = tatm_newq
        
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
        tatm_Boilup = -tatm_bComposition[0]/tatm_sLine[1]
        tatm_L = tatm_dFeed * tatm_Reflux
        tatm_VBar = tatm_bFeed * tatm_Boilup
        tatm_VapHeatCondenser = (tatm_dComposition[0]*THF_VapHeat)+(tatm_dComposition[1]*MEOH_VapHeat)
        tatm_VapHeatEvaporator = (tatm_bComposition[0]*THF_VapHeat)+(tatm_bComposition[1]*MEOH_VapHeat)
        tatm_Condenser = tatm_L * tatm_VapHeatCondenser * 0.0002777778 #convertion of kJ/kmol to kJ/second (aka kW)
        tatm_Evaporator = tatm_VBar * tatm_VapHeatEvaporator * 0.0002777778
        
        
        ##########################################
        ##               Step Four              ##
        ##########################################
        
        if(tatm_dFeed - RecycleFeed < 0.5 and tatm_dFeed - RecycleFeed > -0.5):
            converged = True
            RecycleFeed = tatm_dFeed
            FeedArr.append(RecycleFeed)
            break
        else:
            RecycleFeed = tatm_dFeed
            FeedArr.append(RecycleFeed)
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
    
    plt.figtext(0.13, 0.22, "Distillate", fontsize='x-large', fontweight='bold')
    plt.figtext(0.13, 0.19, "Flow = " + str(oatm_dFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.13, 0.16, "zTHF = " + str(oatm_dComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.13, 0.13, "zMEOH = " + str(1-oatm_dComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.13, 0.10, "Reflux Ratio = " + str(oatm_Reflux))
    plt.figtext(0.13, 0.07, "Condenser Duties =" + str(oatm_Condenser)[:str(oatm_Condenser).index('.')]+"[kW]")
    
    plt.figtext(0.3, 0.22, "Bottoms", fontsize='x-large', fontweight='bold')
    plt.figtext(0.3, 0.19, "Flow = " + str(oatm_bFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.3, 0.16, "zTHF = " + str(oatm_bComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.3, 0.13, "zMEOH = " + str(1-oatm_bComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.3, 0.10, "Boilup Ratio = " + str(oatm_Boilup)[0:4])
    plt.figtext(0.3, 0.07, "Evaporator Duties ="+ str(oatm_Evaporator)[:str(oatm_Evaporator).index('.')] +"[kW]")
    
    plt.figtext(0.55, 0.22, "Recycle", fontsize='x-large', fontweight='bold')
    plt.figtext(0.55, 0.19, "Flow = " + str(tatm_dFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.55, 0.16, "zTHF = " + str(tatm_dComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.55, 0.13, "zMEOH = " + str(1-tatm_dComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.55, 0.10, "Reflux Ratio = " + str(tatm_Reflux))
    plt.figtext(0.55, 0.07, "Condenser Duties = "+ str(tatm_Condenser)[:str(tatm_Condenser).index('.')]+"[kW]")
    
    plt.figtext(0.72, 0.22, "Bottoms", fontsize='x-large', fontweight='bold')
    plt.figtext(0.72, 0.19, f"Flow = " + str(tatm_bFeed)[0:5] +"[kmol/h]")
    plt.figtext(0.72, 0.16, f"zTHF = " + str(tatm_bComposition[0])[0:6] + "[molTHF/mol]")
    plt.figtext(0.72, 0.13, f"zMEOH = " + str(1-tatm_bComposition[0])[0:6] + "[molMEOH/mol]")
    plt.figtext(0.72, 0.10, "Boilup Ratio = " + str(tatm_Boilup)[0:4])
    plt.figtext(0.72, 0.07, "Evaporator Duties = " + str(tatm_Evaporator)[:str(tatm_Evaporator).index('.')] +"[kW]")
    
    
    img = np.array(Image.open('img/TestDiagram.png'))
    fig2, ax2 = plt.subplots(1,1,figsize=(14, 8))
    fig2.subplots_adjust(right=0.8)
    ax2.axis('off')
    #inlet
    plt.figtext(0.07, 0.58, str(Feed) + " [kmol/h]", fontsize = 'small')
    plt.figtext(0.07, 0.55, str(Feed_zTHF) + " [molTHF/mol]", fontsize = 'small')
    plt.figtext(0.07, 0.52, str(Feed_zMEOH) + " [molMEOH/mol]", fontsize = 'small')
    
    #col1 bottoms
    plt.figtext(0.46, 0.1, str(oatm_bFeed)[0:5] +"[kmol/h]", fontsize = 'small')
    plt.figtext(0.46, 0.07, str(oatm_bComposition[0])[0:6] + "[molTHF/mol]", fontsize = 'small')
    plt.figtext(0.46, 0.04, str(1-oatm_bComposition[0])[0:6] + "[molMEOH/mol]", fontsize = 'small')
    
    #col1 distills
    plt.figtext(0.52, 0.73, str(oatm_dFeed)[0:5] +"[kmol/h]", fontsize = 'small')
    plt.figtext(0.52, 0.7, str(oatm_dComposition[0])[0:6] + "[molTHF/mol]", fontsize = 'small')
    plt.figtext(0.52, 0.67, str(1-oatm_dComposition[0])[0:6] + "[molMEOH/mol]", fontsize = 'small')
    
    #col2 bottoms
    plt.figtext(0.75, 0.1, str(oatm_dFeed)[0:5] +"[kmol/h]", fontsize = 'small')
    plt.figtext(0.75, 0.07, str(oatm_dComposition[0])[0:6] + "[molTHF/mol]", fontsize = 'small')
    plt.figtext(0.75, 0.04, str(1-oatm_dComposition[0])[0:6] + "[molMEOH/mol]", fontsize = 'small')
    
    #col2 distils
    plt.figtext(0.78, 0.73, str(tatm_dFeed)[0:5] +" [kmol/h]", fontsize = 'small')
    plt.figtext(0.78, 0.7, str(tatm_dComposition[0])[0:6] + " [molTHF/mol]", fontsize = 'small')
    plt.figtext(0.78, 0.67, str(1-tatm_dComposition[0])[0:6] + " [molMEOH/mol]", fontsize = 'small')
    
    #mixing point
    plt.figtext(0.25, 0.47,str(Feed + RecycleFeed)[0:5] +" [kmol/h]", fontsize = 'small')
    plt.figtext(0.25, 0.44,str(oatm_zTHF)[0:6] + " [molTHF/mol]", fontsize = 'small')
    plt.figtext(0.25, 0.41,str(1-oatm_zTHF)[0:6] + " [molMEOH/mol]", fontsize = 'small')
    
    #condenser and evaporots
    plt.figtext(0.423, 0.265, str(oatm_Evaporator)[:str(oatm_Evaporator).index('.')] +"[kW]", fontsize = 'small') #col1eva
    plt.figtext(0.423, 0.695, str(oatm_Condenser)[:str(oatm_Condenser).index('.')] +"[kW]", fontsize = 'small') #col1cond
    plt.figtext(0.72, 0.265, str(tatm_Evaporator)[:str(tatm_Evaporator).index('.')] +"[kW]", fontsize = 'small') #col2eva
    plt.figtext(0.72, 0.695, str(tatm_Condenser)[:str(tatm_Condenser).index('.')] +"[kW]", fontsize = 'small') #col2con
    plt.figtext(0.29, 0.52, str(oatm_EntryDuty)[:str(oatm_EntryDuty).index('.')] +"[kW]", fontsize = 'small') #entry
    
    imgplot = plt.imshow(img)
    
    fig3, ax3 = plt.subplots(1,1, figsize=(14, 8))
    FeedArr = [i/max(FeedArr) for i in FeedArr]
    ax3.plot(range(1,len(FeedArr)+1), FeedArr)
    ax3.set_title("Final Recycle Flow Rate Ratio with Respect to Iteration Number")
    ax3.set_xlabel("Iteration [#]")
    ax3.set_ylabel("Ratio of Final Flow []")
    
    plt.show()