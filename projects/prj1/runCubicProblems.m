function [] = runCubicProblems(CubicFun)

tol = 1e-6;

load cubicProblems.mat;
n = length(problems);
% disp(n);

disp(sprintf('-------- FIRST TEST SET --------'));
for i=1:n
    coeff = problems{i}.coef;
%     disp(coeff);
    exact = problems{i}.rts;
%     disp(exact);
    solution = CubicFun(coeff); if iscolumn(solution), solution=solution'; end
    
    disp(sprintf('Test problem %d',i));
    [errs,rrts,rexact,flg] = check_roots(solution,exact);
    if (flg == 1)
        disp('   Incorrect number of solutions.');
        disp(' ');
        disp('   CubicFun solutions:');
        disp(solution);
        disp('   GSI provided solutions:');
        disp(exact);
%         disp('   Inputs:');
%         disp(coeff);
    elseif (flg == 2)
        disp('   Some roots are NaN.');
        disp(' ');
        disp('   CubicFun solutions:');
        disp(solution);
        disp('   GSI provided solutions:');
%         disp(exact);
%         disp('   Inputs:');
        disp(coeff);
    elseif (max(abs(errs)) > tol)
        disp(sprintf('   Solutions not within acceptable rel err tol. Max rel err = %f',max(abs(errs))));
        disp(' ');
        disp('   CubicFun solutions:');
        disp(rexact);
        disp('   GSI provided solutions:');
        disp(rrts);
        disp('   Errors:');
        disp(errs);
%         disp('   Inputs:');
%         disp(coeff);
    else
        disp('  Test passed.');
    end
    disp(' ');
end

load cubicProblems2.mat
n = length(problems);

disp(sprintf('-------- SECOND TEST SET --------'));
for i=1:n
    coeff = problems{i}.coef;
%     disp(coeff);
    exact = problems{i}.rts;
%     disp(exact);
    solution = CubicFun(coeff); if iscolumn(solution), solution=solution'; end
    
    disp(sprintf('Test problem %d',i));
    [errs,rrts,rexact,flg] = check_roots(solution,exact);
    if (flg == 1)
        disp('   Incorrect number of solutions.');
        disp(' ');
        disp('   CubicFun solutions:');
        disp(solution);
        disp('   GSI provided solutions:');
        disp(exact);
    elseif (flg == 2)
        disp('   Some roots are NaN.');
        disp(' ');
        disp('   CubicFun solutions:');
        disp(solution);
        disp('   GSI provided solutions:');
        disp(exact);
    elseif (max(abs(errs)) > tol)
        disp(sprintf('   Solutions not within acceptable rel err tol. Max rel err = %f',max(abs(errs))));
        disp(' ');
        disp('   CubicFun solutions:');
        disp(rexact);
        disp('   GSI provided solutions:');
        disp(rrts);
        disp('   Errors:');
        disp(errs);
    else
        disp('  Test passed.');
    end
    disp(' ');
end


end