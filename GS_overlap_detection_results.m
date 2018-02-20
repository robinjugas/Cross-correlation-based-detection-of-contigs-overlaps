function [accuracy,sensitivity,positive_predictivity,joint,disjoint,missing] = GS_overlap_detection_results( reference, result )
%
N=length(reference);
TP=0;    %correctly detected overlaps
FP=0;   %incorretly detected overlaps
TN=0;    %correctly undetected overlaps
FN=abs(length(reference)-length(result)); % missing undetected overlaps
ovrlps_num=factorial(N) / factorial(N-2); %total overlaps ~ "population"

%searching TP
for i=2:length(result)
    if result(i)==result(i-1)+1 || result(i)==result(i-1)-1
        TP=TP+1;    
    else
        FP=FP+1;
    end
end

%TN=n-others
TN=ovrlps_num-TP-FP-FN;
%ACC=TP+TN/total number of overlaps
accuracy=(TP+TN)/ovrlps_num;
%sensitivity
sensitivity=TP/(TP+FN);
%positive predictive  value
positive_predictivity=TP/(TP+FP);

%precentage of correct and incorrect joins
joint=TP;
disjoint=FP;
missing=FN;









