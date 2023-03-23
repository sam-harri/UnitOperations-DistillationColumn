import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

oatmData = pd.read_csv("Data/oatmData.csv")
tatmData = pd.read_csv("Data/tatmData.csv")

fig, ax = plt.subplots(1,2,figsize=(14,8))

boilupTwin = ax[0].twinx()
ax[0].set_title("1atm Column")
ax[0].set_xlabel("Reflux Ratio []")
p01, = ax[0].plot(oatmData['reflux'], oatmData['condenser'], label = 'Condenser Duty', c='dodgerblue')
p02, = ax[0].plot(oatmData['reflux'], oatmData['evaporator'], label = 'Evaporator Duty', c='royalblue')
p03, = boilupTwin.plot(oatmData['reflux'], oatmData['boilup'], label = 'Boilup Ratio', c='lightgrey')
p04, = boilupTwin.plot(oatmData['reflux'], oatmData['steps'], label = 'Steps', c='dimgrey')

ax[0].set_ylabel('Duties [kW]')
ax[0].legend(handles =[p01,p02,p03,p04])

boilupTwin2 = ax[1].twinx()
ax[1].set_title("10atm Column")
ax[1].set_xlabel("Reflux Ratio []")
p11, = ax[1].plot(tatmData['reflux'], tatmData['condenser'], label = 'Condenser Duty', c='crimson')
p12, = ax[1].plot(tatmData['reflux'], tatmData['evaporator'], label = 'Evaporator Duty', c='lightcoral')
p13, = boilupTwin2.plot(tatmData['reflux'], tatmData['boilup'], label = 'Boilup Ratio', c='lightgrey')
p14, = boilupTwin2.plot(tatmData['reflux'], tatmData['steps'], label = 'Steps', c='dimgrey')

ax[1].set_ylabel('Duties [kW]')
ax[1].legend(handles =[p11,p12,p13,p14])

plt.show()