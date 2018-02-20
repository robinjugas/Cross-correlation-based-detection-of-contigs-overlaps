function [corr] = pearson(A,B)
% https://stackoverflow.com/questions/5644981/pearsons-coefficient-and-covariance-calculation-in-matlab
% Elapsed time is 0.013214 seconds. (2x50 vzorku)
% C=cov(A,B);
% p=C(2)/(std(A)*std(B));

% https://www.mathworks.com/help/matlab/ref/corrcoef.html
% Elapsed time is 0.009052 seconds.(2x50 vzorku)

corr=1/(length(A)-1)*sum(((A-mean(A))/std(A)).*((B-mean(B))/std(B)));

% tic
% xa = A - sum(A)./length(A);
% xb = B - sum(B)./length(B);
% c = sum(xa.*xb)./ (length(xa)-1);
% corr = c/(std(A)*std(B));
% toc

end

