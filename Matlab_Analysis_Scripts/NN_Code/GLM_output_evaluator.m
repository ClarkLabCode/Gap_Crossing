% minDevLoc = FitInfo.IndexMinDeviance;
minDevLoc = 1;
GLM_Test_Pred = ((ParamWeights(:,minDevLoc)'*GLM_Input_Mat_Rand_Test')+FitInfo.Intercept(minDevLoc))>=0;
PredAndTestDiff = sum(abs(GLM_Test_Pred'-GLM_Resp_Vec_Rand_Test));
PredSuccRate = 1 - PredAndTestDiff/length(GLM_Test_Pred);