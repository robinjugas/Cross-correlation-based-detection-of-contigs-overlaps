function [order] = overlap_detection_NN(seq)
%% INICIALIZACE

N=length(seq);
%% Overlap detection

correlation_matrix=zeros(N,N);
lags_matrix=zeros(N,N);
for i=1:N
    for j=1:N
        if j > i
            
            %který ze 2 readù je delší?
            if length(seq{i})>=length(seq{j})
                sig_delsi=seq{i};
                sig_kratsi=seq{j};
                LD=length(sig_delsi);LK=length(sig_kratsi);
            else
                sig_delsi=seq{j};
                sig_kratsi=seq{i};
                LD=length(sig_delsi);LK=length(sig_kratsi);
            end
            
            kor=0; %vektor hodnot korelací
            
            % for cyklus pro levý pøekryv
            for a=1:LK
                x=sig_kratsi(end-a+1:end);
                y=sig_delsi(1:a);
                kor=[kor pearson(x,y)];
            end
            kor(1:10)=0; %nulovaní prvních 10 vzorkù
            
            %for cyklus pro kratší v delším, jen pokud se liší v délce,
            %jinak beze smyslu
            if LD==LK
            else
                for b=1:abs(LD-LK)
                    x=sig_kratsi; %beru celý
                    y=sig_delsi(b:b+LK-1); %posouvám jen delší
                    kor=[kor pearson(x,y)];
                end
            end
            
            % for cyklus pro pravý pøekryv
            for c=1:LK
                x=sig_kratsi(1:end-c+1);
                y=sig_delsi(end-LK+c:end);
                kor=[kor pearson(x,y)];
                
            end
            kor(end-10:end)=0; %nulovaní posledních 10 vzorkù
            kor(1)=[];
            
            % uložím do matice
            [correlation_matrix(i,j),lags_matrix(i, j)]=max(kor);
            % uložit si pozice 
            
        else
            correlation_matrix(i, j) = correlation_matrix(j, i);
            lags_matrix(i, j) = lags_matrix(j, i);
        end
    end
end

%% metoda NN nejbližší soused
% symetrická matice
reseni=zeros(N,N); % zde budu ukladat jednotliva reseni (na kazdem radku jedno reseni)

for j=1:N
    % v kazdem kroku tohoto cyklu bude provedena metoda nejblizsiho souseda
    % pro pocatecni mesto = ind1
    matice=correlation_matrix;
    poc_pozice=j;
    akt_pozice=j;
    cesta=[poc_pozice];
    Skore=0;
    for i=1:N-1
            %najdu maxima %ošetøit, pokud jich je víc (možný exploit)
            [korelace_radek,pozice_radek]=max(abs(matice(akt_pozice,:))); %najdu maximum
            akt_pozice=pozice_radek;             
            cesta=[cesta akt_pozice]; %pøidám do cesty            
            Skore(end+1)=korelace_radek; %uložím skore=korelaci            
            
            %VYMAZAT ODPOVÍDAJÍCÍ ØÁDEK A SLOUPEC
            matice(cesta(end-1),:)=NaN;
            matice(:,cesta(end-1))=NaN;            
            
    end
    reseni(j,:)=cesta;
    SkoreReseni(j)=sum(Skore);
    
end

% find best solution
[hodnota, pozice]=max(SkoreReseni);
order=reseni(pozice,:);


end %konec funkce