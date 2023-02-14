import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_validate
from sklearn.model_selection import train_test_split

VLE1atm = pd.read_csv('VLE/1atmVLEData.csv')
VLE10atm = pd.read_csv('VLE/10atmVLEData.csv')

x_1atm = np.array(VLE1atm['x1'])
y_1atm = np.array(VLE1atm['y1'])
x_1atm = x_1atm.reshape(-1, 1)
x_train_1atm, x_test_1atm, y_train_1atm, y_test_1atm = train_test_split(x_1atm, y_1atm, train_size=0.8)

x_10atm = np.array(VLE10atm['x1'])
y_10atm = np.array(VLE10atm['y1'])
x_10atm = x_10atm.reshape(-1, 1)
x_train_10atm, x_test_10atm, y_train_10atm, y_test_10atm = train_test_split(x_10atm, y_10atm, train_size=0.8)

maxdegree = 7

training_error_1atm = []
cross_validation_error_1atm = []
for d in range(1, maxdegree):
    x_poly_train_1atm = PolynomialFeatures(degree=d).fit_transform(x_train_1atm)
    x_poly_test_1atm = PolynomialFeatures(degree=d).fit_transform(x_test_1atm)
    lr_1atm = LinearRegression(fit_intercept=False)
    model_1atm = lr_1atm.fit(x_poly_train_1atm, y_train_1atm)
    y_train_pred_1atm = model_1atm.predict(x_poly_train_1atm)
    mse_train_1atm = mean_squared_error(y_train_1atm, y_train_pred_1atm)
    cve_1atm = cross_validate(lr_1atm, x_poly_train_1atm, y_train_1atm, scoring='neg_mean_squared_error', cv=5,return_train_score=True)
    training_error_1atm.append(mse_train_1atm)
    cross_validation_error_1atm.append(np.mean(np.absolute(cve_1atm['test_score'])))

training_error_10atm = []
cross_validation_error_10atm = []
for d in range(1, maxdegree):
    x_poly_train_10atm = PolynomialFeatures(degree=d).fit_transform(x_train_10atm)
    x_poly_test_10atm = PolynomialFeatures(degree=d).fit_transform(x_test_10atm)
    lr_10atm = LinearRegression(fit_intercept=False)
    model_10atm = lr_10atm.fit(x_poly_train_10atm, y_train_10atm)
    y_train_pred_10atm = model_10atm.predict(x_poly_train_10atm)
    mse_train_10atm = mean_squared_error(y_train_10atm, y_train_pred_10atm)
    cve_10atm = cross_validate(lr_10atm, x_poly_train_10atm, y_train_10atm, scoring='neg_mean_squared_error', cv=5, return_train_score=True)
    training_error_10atm.append(mse_train_10atm)
    cross_validation_error_10atm.append(np.mean(np.absolute(cve_10atm['test_score'])))

poly_features_1atm = PolynomialFeatures(degree=5, include_bias=False).fit_transform(x_1atm)
poly_reg_model_1atm = LinearRegression(fit_intercept=False)
poly_reg_model_1atm.fit(poly_features_1atm, y_1atm)
y_predicted_1atm = poly_reg_model_1atm.predict(poly_features_1atm)

poly_features_10atm = PolynomialFeatures(degree=6, include_bias=False).fit_transform(x_10atm)
poly_reg_model_10atm = LinearRegression(fit_intercept=False)
poly_reg_model_10atm.fit(poly_features_10atm, y_10atm)
y_predicted_10atm = poly_reg_model_10atm.predict(poly_features_10atm)

print(poly_reg_model_10atm.coef_)
print(poly_reg_model_10atm.intercept_)


fig, ax = plt.subplots(2, 2, figsize=(14, 8))
fig.tight_layout()
plt.subplots_adjust(top=0.95, bottom=0.07)

ax[0][0].plot(VLE1atm['x1'], VLE1atm['y1'], label='VLE Data')
ax[0][0].plot(VLE1atm['x1'], y_predicted_1atm, label='VLE Model')
ax[0][0].plot([0, 1], [0, 1], label='Azeotropic Line', linestyle='dashed')
ax[0][0].set_title('1atm VLE Data and Model')
# ax[0][0].set_xlabel('x1')
# ax[0][0].set_ylabel('y1')
ax[0][0].legend()
ax[0][1].plot(range(1, maxdegree), cross_validation_error_1atm)
ax[0][1].set_title('5-Fold 1atm Model CV')
# ax[0][1].set_xlabel('Polynomial Degree')
# ax[0][1].set_ylabel('MSE')

ax[1][0].plot(VLE10atm['x1'], VLE10atm['y1'], label='VLE Data')
ax[1][0].plot(VLE10atm['x1'], y_predicted_10atm, label='VLE Model')
ax[1][0].plot([0, 1], [0, 1], label='Azeotropic Line', linestyle='dashed')
ax[1][0].set_title('10atm VLE Data and Model')
# ax[1][0].set_xlabel('x1')
# ax[1][0].set_ylabel('y1')
ax[1][0].legend()
ax[1][1].plot(range(1, maxdegree), cross_validation_error_10atm)
ax[1][1].set_title('5-Fold 10atm Model CV')
# ax[1][1].set_xlabel('Polynomial Degree')
# ax[1][1].set_ylabel('MSE')

plt.show()
