function [signaly] = genomic_signal(fasta,reprezentace); %resample1,resample2)
%resample1=0 nic se neprovádí
fasta=upper(fasta);
N=length(fasta);

%% Converting into signals
signaly=cell(N,1);
switch reprezentace
    case 'kumulovana'
        for i=1:N
            signaly{i,1}=kumulovanafaze(fasta{i});
        end
    case 'rozbalena'
        for i=1:N
            signaly{i,1}=rozbalenafaze(fasta{i});
        end
    case 'faze'
        for i=1:N
            signaly{i,1}=faze(fasta{i});
        end
end

% %% resample
% for i=1:N
%     signaly{i,1}=resample( signaly{i,1},resample1,resample2);
% end

end

