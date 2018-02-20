function [ signal ] = rozbalenafaze( fasta )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
a=1+1i;c=-1-1i;g=-1+1i;t=1-1i; %numerical representation formula
zapis=zeros(1,length(fasta),'double');
zapis(fasta=='A')=angle(a);
zapis(fasta=='C')=angle(c);
zapis(fasta=='G')=angle(g);
zapis(fasta=='T')=angle(t);
signal=unwrap(zapis); %variable of numerical reads coming into overlap detection
%     SequencesNumS{i,1}=SequencesNumS{i,1}-mean(SequencesNumS{i,1});


end

