B = [-1,0;-1/2,1;1/2,1;1,0];
vred = bezierPts(B,1/3)

figure;
hold on; axis equal;
plot(B(:,1),B(:,2),".-",'MarkerSize',15);
vecPts = bezierPts(B,linspace(0,1));
plot(vecPts(:,1),vecPts(:,2));
hold off;

function b = bezierPts(B,t)
	% Opis:
	%  bezier vrne tocke na Bezierjevi krivulji pri danih parametrih
	%
	% Definicija:
	%  b = bezier(B,t)
	%
	% Vhodna podatka:
	%  B    matrika velikosti n+1 x d, ki predstavlja kontrolne tocke
	%       Bezierjeve krivulje stopnje n v d-dimenzionalnem prostoru,
	%  t    seznam parametrov dolzine k, pri katerih racunamo vrednost
	%       Bezierjeve krivulje
	%
	% Izhodni podatek:
	%  b    matrika velikosti k x d, kjer i-ta vrstica predstavlja tocko na
	%       Bezierjevi krivulji pri parametru iz t na i-tem mestu
	dim = size(B);
	k = size(t,2);
	b = zeros(k,dim(2));
	for i=1:k
		D = decasteljau(B,t(i));
		b(i,:) = D{1,dim};
	end
end

function D = decasteljau(b,t)
	% Opis:
	%  decasteljau vrne shemo de Casteljaujevega postopka za dan seznam
	%  koordinat b pri danem parametru t
	%
	% Definicija:
	%  D = decasteljau(b,t)
	%
	% Vhodna podatka:
	%  b    seznam koordinat kontrolnih tock Bezierjeve krivulje stopnje n,
	%  t    parameter, pri katerem racunamo koordinato (Bezierjeve krivulje)
	%
	% Izhodni podatek:
	%  D    tabela velikosti n+1 x n+1, ki predstavlja de Casteljaujevo shemo
	%       za koordinate b pri parametru t (element na mestu (1,n+1) je
	%       koordinata Bezierjeve krivulje pri parametru t)

	stTock = size(b,1);
	D = cell(stTock);
	for i=1:stTock
		D{i,1} = b(i,:);
	end

	for r=1:(stTock-1)
		for i=0:(stTock-r-1)
			D{i + 1,r + 1} = (1 - t)*D{i + 1,r} + t*D{i + 2,r};
		end
	end
end
	