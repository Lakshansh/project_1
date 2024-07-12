% Solving for H2O
dataH2O = readmatrix("C:\Users\laksh\OneDrive\Desktop\MATLAB\data_pvap.xlsx", 'sheet', 'H2O');
dataPartH2O = cvpartition(size(dataH2O, 1), 'HoldOut', 0.3);
idxH2O = dataPartH2O.test;
trainDataH2O = dataH2O(~idxH2O, :);
testDataH2O = dataH2O(idxH2O, :);
P_H2O = trainDataH2O(:, 2);
T_H2O = trainDataH2O(:, 1);
logP_H2O = log(P_H2O);
invT_H2O = 1 ./ T_H2O;
lm_H2O = fitlm(invT_H2O, logP_H2O, 'poly1');
coef_H2O = lm_H2O.Coefficients.Estimate;
A_H2O = coef_H2O(1);
B_H2O = -coef_H2O(2);
logPCalcH2O = (-B_H2O .* invT_H2O + A_H2O);
diffH2O = logP_H2O - logPCalcH2O;
rmseTrainH2O = sqrt(sum(diffH2O .* diffH2O) / numel(diffH2O));
invTTestH2O = 1 ./ testDataH2O(:, 1);
logPTestH2O = log(testDataH2O(:, 2));
logPCalcTestH2O = (-B_H2O .* invTTestH2O + A_H2O);
diffTestH2O = logPTestH2O - logPCalcTestH2O;
rmseTestH2O = sqrt(sum(diffTestH2O .* diffTestH2O) / numel(diffTestH2O));
logPCalcAllH2O = -B_H2O .* (1 ./ dataH2O(:, 1)) + A_H2O;
figure(1);
plot(dataH2O(:, 1), log(dataH2O(:, 2)));
hold on;
plot(dataH2O(:, 1), logPCalcAllH2O);
xlabel('T');
ylabel('log(P)');
legend('Data', 'Calculated');
title('P vs T for H2O');
hold off;

% Solving for CH3OH
dataCH3OH = readmatrix("C:\Users\laksh\OneDrive\Desktop\MATLAB\data_pvap.xlsx", 'sheet', 'CH3OH');
dataPartCH3OH = cvpartition(size(dataCH3OH, 1), 'HoldOut', 0.3);
idxCH3OH = dataPartCH3OH.test;
trainDataCH3OH = dataCH3OH(~idxCH3OH, :);
testDataCH3OH = dataCH3OH(idxCH3OH, :);
P_CH3OH = trainDataCH3OH(:, 2);
T_CH3OH = trainDataCH3OH(:, 1);
logP_CH3OH = log(P_CH3OH);
invT_CH3OH = 1 ./ T_CH3OH;
lm_CH3OH = fitlm(invT_CH3OH, logP_CH3OH, 'poly1');
coef_CH3OH = lm_CH3OH.Coefficients.Estimate;
A_CH3OH = coef_CH3OH(1);
B_CH3OH = -coef_CH3OH(2);
logPCalcCH3OH = (-B_CH3OH .* invT_CH3OH + A_CH3OH);
diffCH3OH = logP_CH3OH - logPCalcCH3OH;
rmseTrainCH3OH = sqrt(sum(diffCH3OH .* diffCH3OH) / numel(diffCH3OH));
invTTestCH3OH = 1 ./ testDataCH3OH(:, 1);
logPTestCH3OH = log(testDataCH3OH(:, 2));
logPCalcTestCH3OH = (-B_CH3OH .* invTTestCH3OH + A_CH3OH);
diffTestCH3OH = logPTestCH3OH - logPCalcTestCH3OH;
rmseTestCH3OH = sqrt(sum(diffTestCH3OH .* diffTestCH3OH) / numel(diffTestCH3OH));
logPCalcAllCH3OH = -B_CH3OH .* (1 ./ dataCH3OH(:, 1)) + A_CH3OH;
figure(2);
plot(dataCH3OH(:, 1), log(dataCH3OH(:, 2)));
hold on;
plot(dataCH3OH(:, 1), logPCalcAllCH3OH);
xlabel('T');
ylabel('log(P)');
legend('Data', 'Calculated');
title('P vs T for CH3OH');
hold off;

% Solving for the second part (RK and PR)
T = 210 + 273;
P = 78 * 1.013;
R = 0.08206 / 1000;
Pc_H2 = 12.93;
Tc_H2 = 33.18;
w_H2 = -0.22;
Tc_CO2 = 304.15;
w_CO2 = 0.225;
Pc_CO2 = 7.38;
RK_sigma = 1;
RK_epsilon = 0;
RK_w = 0.08664;
RK_psi = 0.42748;
PR_sigma = 1+sqrt(2);
PR_epsilon = 1-sqrt(2);
PR_w = 0.0778;
PR_psi = 0.45724;
Tr_H2 = T/Tc_H2;
Tr_CO2 = T/Tc_CO2;
Pr_H2 = P/Pc_H2;
Pr_CO2 = P/Pc_CO2;
bt_CO2_RK = RK_w*Pr_CO2/Tr_CO2;
bt_H2_RK = RK_w*Pr_H2/Tr_H2;
bt_CO2_PR = PR_w*Pr_CO2/Tr_CO2;
bt_H2_PR = PR_w*Pr_H2/Tr_H2;
al_CO2_RK = RK(Tr_CO2);
al_H2_RK = RK(Tr_H2);
al_CO2_PR = PR(Tr_CO2,w_CO2);
al_H2_PR = PR(Tr_H2,w_H2);
q_CO2_RK = RK_psi*RK(Tr_CO2)/(RK_w*Tr_CO2);
q_H2_RK = RK_psi*RK(Tr_H2)/(RK_w*Tr_H2);
q_CO2_PR = PR_psi*PR(Tr_CO2,w_CO2)/(PR_w*Tr_CO2);
q_H2_PR = PR_psi*PR(Tr_H2,w_H2)/(PR_w*Tr_H2);
b_CO2_RK = R*T*bt_CO2_RK/P;
b_H2_RK = R*T*bt_H2_RK/P;
b_CO2_PR = R*T*bt_CO2_PR/P;
b_H2_PR = R*T*bt_H2_PR/P;
a_CO2_RK = q_CO2_RK*b_CO2_RK*R*T;
a_H2_RK = q_H2_RK*b_H2_RK*R*T;
a_CO2_PR = q_CO2_PR*b_CO2_PR*R*T;
a_H2_PR = q_H2_PR*b_H2_PR*R*T;
options = optimset('Display', 'off');
V_ideal = R*T/P;
V_CO2_RK = fsolve(@(V) mol_vol(R,T,P,b_CO2_RK,a_CO2_RK,V,RK_epsilon,RK_sigma),1,options);
V_H2_RK = fsolve(@(V) mol_vol(R,T,P,b_H2_RK,a_H2_RK,V,RK_epsilon,RK_sigma),1,options);
V_CO2_PR = fsolve(@(V) mol_vol(R,T,P,b_CO2_PR,a_CO2_PR,V,PR_epsilon,PR_sigma),1,options);
V_H2_PR = fsolve(@(V) mol_vol(R,T,P,b_H2_PR,a_H2_PR,V,PR_epsilon,PR_sigma),1,options);
molVolIdeal = num2str(1/V_ideal);
molVolCO2_RK = num2str(1/V_CO2_RK);
molVolCO2_PR = num2str(1/V_CO2_PR);
molVolH2_RK = num2str(1/V_H2_RK);
molVolH2_PR = num2str(1/V_H2_PR);

function al = RK(Tr)
    al = Tr^(-0.5);
end

function al = PR(Tr,w)
    al = (1+(0.37464 + 1.54226*w - 0.26992*w^2)*(1 - Tr^0.5))^2;
end

function eqn = mol_vol(R,T,P,b,a,V,epsilon,sigma)
    eqn = V - R*T/P - b + (a/P)*((V-b)/((V+epsilon*b)*(V+sigma*b)));
end
