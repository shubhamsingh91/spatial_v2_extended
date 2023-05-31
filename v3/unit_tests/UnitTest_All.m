printHeader('CMM'); 
UnitTest_CMM

printHeader('Coriolis'); 
UnitTest_Coriolis

printHeader('Derivatives'); 
UnitTest_Derivatives

printHeader('Derivatives (Floating)'); 
UnitTest_Derivatives_FB

printHeader('Main Dynamics'); 
UnitTest_Dynamics

printHeader('Orientation'); 
UnitTest_Orientation

printHeader('Orientation Rates'); 
UnitTest_OrientationRates

printHeader('ID SO expressions'); 
UnitTest_IDSO_expressions

printHeader('SVA Tensor Identites'); 
UnitTest_SVATensor_iden

printHeader('SVA Tensor Properties'); 
UnitTest_SVATensor_prop

printHeader('Multi-Dof Identites'); 
UnitTest_multiDoF_iden

function printHeader(st)
    fprintf('************************************\n');
    fprintf('%s\n',st);
    fprintf('************************************\n');
end