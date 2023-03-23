from typing import List

class ReactorConstants:
    VLE_1atm : List[float] = [7.44200335, -18.71029508, 18.98056828, -9.66489131, 2.95151852, 0]
    VLE_10atm : List[float] = [3.31669035, -5.09330499, 0.43491977, 4.12563412, -3.3635295, 1.57874402, 0]
    refluxRatio_col1 : float = 1.35
    refluxRatio_col2 : float = 1.35
    azeotrope : List[int] = [1,0]
    zf : float = 0.25
    xb_col1 : float = 0.0100
    xb_col2 : float = 0.9925
    feedRate : float = 110
    xT_1atm = [210.63330852, -622.23229707, 776.82882779, -522.34139505, 213.5536741, -55.04248631, 64.44918253]
    xT_10atm = [109.91581372, -170.78369919, 111.56103898, -28.29272664, 136.53283177]
    THF_VapHeat = 29510 #kJ/kmol
    MEOH_VapHeat = 37340 #kJ/kmol
    Feed_Temperature = 25
    THF_CPLiq = 124.1 #kJ/kmol K
    MEOH_CPLiq = 81.465 #kJ/kmol K
    
    
    ###density
    MEOH_Density_liq = 792 #kg/m^3
    MEOH_Density_vap = 1.44 #kg/m^3
    MEOH_SurfaceTension = 22.07 #mN/m or dyn/cm
    THF_Density_liq = 888 #kg/m^3
    THF_Density_vap = 3.23 #kg/m^3
    THF_SurfaceTension = 26.4 #mN/m or dyn/cm
    Water_HeatCapacity = 75.35 #kJ/kmol K
    Water_Temperature = 12 #C