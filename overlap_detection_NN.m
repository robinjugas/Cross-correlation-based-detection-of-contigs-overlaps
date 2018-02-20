function [order] = overlap_detection_NN(seq)
%% INICIALIZACE

N=length(seq);
%% Overlap detection

correlation_matrix=zeros(N,N);
lags_matrix=zeros(N,N);
for i=1:N
    for j=1:N
        if j > i
            
            %kter� ze 2 read� je del��?
            if length(seq{i})>=length(seq{j})
                sig_delsi=seq{i};
                sig_kratsi=seq{j};
                LD=length(sig_delsi);LK=length(sig_kratsi);
            else
                sig_delsi=seq{j};
                sig_kratsi=seq{i};
                LD=length(sig_delsi);LK=length(sig_kratsi);
            end
            
            kor=0; %vektor hodnot korelac�
            
            % for cyklus pro lev� p�ekryv
            for a=1:LK
                x=sig_kratsi(end-a+1:end);
                y=sig_delsi(1:a);
                kor=[kor pearson(x,y)];
            end
            kor(1:10)=0; %nulovan� prvn�ch 10 vzork�
            
            %for cyklus pro krat�� v del��m, jen pokud se li�� v d�lce,
            %jinak beze smyslu
            if LD==LK
            else
                for b=1:abs(LD-LK)
                    x=sig_kratsi; %beru cel�
                    y=sig_delsi(b:b+LK-1); %posouv�m jen del��
                    kor=[kor pearson(x,y)];
                end
            end
            
            % for cyklus pro prav� p�ekryv
            for c=1:LK
                x=sig_kratsi(1:end-c+1);
                y=sig_delsi(end-LK+c:end);
                kor=[kor pearson(x,y)];
                
            end
            kor(end-10:end)=0; %nulovan� posledn�ch 10 vzork�
            kor(1)=[];
            
            % ulo��m do matice
            [correlation_matrix(i,j),lags_matrix(i, j)]=max(kor);
            % ulo�it si pozice 
            
        else
            correlation_matrix(i, j) = correlation_matrix(j, i);
            lags_matrix(i, j) = lags_matrix(j, i);
        end
    end
end

%% metoda NN nejbli��� soused
% symetrick� matice
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
            %najdu maxima %o�et�it, pokud jich je v�c (mo�n� exploit)
            [korelace_radek,pozice_radek]=max(abs(matice(akt_pozice,:))); %najdu maximum
            akt_pozice=pozice_radek;             
            cesta=[cesta akt_pozice]; %p�id�m do cesty            
            Skore(end+1)=korelace_radek; %ulo��m skore=korelaci            
            
            %VYMAZAT ODPOV�DAJ�C� ��DEK A SLOUPEC
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