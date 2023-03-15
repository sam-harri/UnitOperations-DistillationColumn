#https://www.youtube.com/watch?v=FwN3LHN1vpc
#https://www.youtube.com/watch?v=MeBoAGZ6W6s

from typing import List
import numpy as np
import matplotlib.pyplot as plt
from NumericalMethods.NumericalMetods import NumericalMethods
from ReactorConstants import ReactorConstants

if(__name__=='__main__'):
    convergingRun = True
    #change
    ######################
    q_1atm = 1 #1.2
    xd_1atm = 0.454 #0.454
    refluxRatio_col1 : float = ReactorConstants.refluxRatio_col1
    
    q_10atm = 0.5 #0.5
    xd_10atm = 0.295 #0.295
    refluxRatio_col2 : float = ReactorConstants.refluxRatio_col2
    
    #1atm 
    ######################
    xb_1atm = ReactorConstants.xb_col1
    VLE_1atm : List[float] = ReactorConstants.VLE_1atm
    AzeotropeLine : List[int] = ReactorConstants.azeotrope
    zf_1atm : float = ReactorConstants.zf
    refluxRatio_col1 : float = ReactorConstants.refluxRatio_col1
    
    az_1atm = NumericalMethods.getPolynomialIntercept(VLE_1atm, AzeotropeLine)[0]
    qLine_1atm : List[float] = NumericalMethods.qLinePolynomial(q_1atm, zf_1atm)
    oLine_1atm : List[float] = NumericalMethods.operatingLinePolynomial(refluxRatio_col1, xd_1atm)
    sLine_1atm : List[float] = NumericalMethods.strippingLinePolynomial(qLine_1atm, oLine_1atm, xb_1atm)
    x_qo_1atm, y_qo_1atm = NumericalMethods.getPolynomialIntercept(qLine_1atm, oLine_1atm)
    
    # if(y_qo_1atm > np.polyval(VLE_1atm, x_qo_1atm)):
    #     print("Limiting case [LowP Col]: Intersection of oLine, qLine, and sLine is above the VLE Line")
    #     convergingRun = False
    # else:
    #     steps_1atm, data_1atm = NumericalMethods.staircaseData1atm(oLine_1atm, sLine_1atm, VLE_1atm, xd_1atm, xb_1atm)
    steps_1atm, data_1atm, lines = NumericalMethods.staircaseData1atm(q_1atm, zf_1atm, refluxRatio_col1, VLE_1atm, xd_1atm, xb_1atm)
    
    
    #10atm
    ######################
    xb_10atm = ReactorConstants.xb_col2
    print(xb_10atm)
    print(xd_10atm)
    VLE_10atm : List[float] = ReactorConstants.VLE_10atm
    zf_10atm : float = 0.5*(data_1atm[-1][0]+data_1atm[-1][1])
    az_10atm = NumericalMethods.getPolynomialIntercept(VLE_10atm, AzeotropeLine)[0]
    
    qLine_10atm : List[float] = NumericalMethods.qLinePolynomial(q_10atm, zf_10atm)
    oLine_10atm : List[float] = NumericalMethods.operatingLinePolynomial(refluxRatio_col2, xd_10atm)
    sLine_10atm : List[float] = NumericalMethods.strippingLinePolynomial(qLine_10atm, oLine_10atm, xb_10atm)
    x_qo_10atm, y_qo_10atm = NumericalMethods.getPolynomialIntercept(qLine_10atm, oLine_10atm)
    
    if(y_qo_10atm < np.polyval(VLE_10atm, x_qo_10atm)):
        print("Limiting case [HighP Col]: Intersection of oLine, qLine, and sLine is above the VLE Line")
        convergingRun = False
    else:
        steps_10atm, data_10atm = NumericalMethods.staircaseData10atmp(oLine_10atm, sLine_10atm, VLE_10atm, xd_10atm, xb_10atm)
    
    if(convergingRun):
        fig, ax = plt.subplots(1,2)
        
        xvals_1atm = np.linspace(0,1) #az_1atm
        ax[0].set_title('1atm Column')
        ax[0].plot(xvals_1atm, np.polyval(AzeotropeLine, xvals_1atm), label='Azeotrope Line') #y=x line
        ax[0].plot(xvals_1atm,np.polyval(VLE_1atm, xvals_1atm), label='VLE Line') #VLE line
        ax[0].plot([zf_1atm, x_qo_1atm],[zf_1atm, y_qo_1atm], label='q Line') #qline
        ax[0].plot([xd_1atm, x_qo_1atm],[xd_1atm, y_qo_1atm], label='Operating Line') #operatingLine
        ax[0].plot([xb_1atm, x_qo_1atm],[xb_1atm, y_qo_1atm], label='Stripping Line') #strippingLine
        ax[0].legend()
        
        i=0
        while(i<len(data_1atm)-1):
            ax[0].plot([data_1atm[i][0], data_1atm[i+1][0]],[data_1atm[i][1],data_1atm[i+1][1]],linestyle='dashed', c='black')
            i+=1
        
        xvals_10atm = np.linspace(0, 1) #az_10atm
        ax[1].set_title('10atm Column')
        ax[1].plot(xvals_10atm, np.polyval(AzeotropeLine, xvals_10atm), label='Azeotrope Line') #y=x line
        ax[1].plot(xvals_10atm,np.polyval(VLE_10atm, xvals_10atm), label='VLE Line') #VLE line
        
        # ax[1].plot(np.linspace(zf_10atm, x_qo_10atm, 2), np.polyval(qLine_10atm, np.linspace(zf_10atm, x_qo_10atm, 2)),  label='q Line') #qline
        # ax[1].plot(np.linspace(xb_10atm, x_qo_10atm,2),np.polyval(sLine_10atm, np.linspace(xb_10atm, x_qo_10atm,2)), label='Operating Line') #strippingLine
        # ax[1].plot(np.linspace(xd_10atm, x_qo_10atm, 2),np.polyval(oLine_10atm, np.linspace(xd_10atm, x_qo_10atm, 2)), label='Stripping Line') #operatingLine
        
        ax[1].plot([zf_10atm, x_qo_10atm],[zf_10atm, y_qo_10atm], label='q Line') #qline
        ax[1].plot([xd_10atm, x_qo_10atm],[xd_10atm, y_qo_10atm], label='Operating Line') #operatingLine
        ax[1].plot([xb_10atm, x_qo_10atm],[xb_10atm, y_qo_10atm], label='Stripping Line') #strippingLine
        
        ax[1].legend()
        
        i=0
        while(i<len(data_10atm)-1):
            ax[1].plot([data_10atm[i][0], data_10atm[i+1][0]],[data_10atm[i][1],data_10atm[i+1][1]],linestyle='dashed', c='black')
            i+=1
            
    

plt.show()