format long;

f = @(x) 1./(x + 1);
fOdvod = @(x) -1./(x + 1).^2;
x0 = 1/2;
hs = [.01,.1,.15,.2];
stH = size(hs,2);
tocnaVred = fOdvod(x0);

METODA = struct("prema",0,"obratna",1,"simetricna",2);

pravilo = zeros(1,stH);
prema = zeros(1,stH);
obratna = zeros(1,stH);
simetricna = zeros(1,stH);
for i = 1:stH
	h = hs(i);
	pravilo(i) = praviloZaOdv(f,h,x0);
	prema(i) = diferencaFunkcije(f,1,h,x0,METODA.prema);
	obratna(i) = diferencaFunkcije(f,1,h,x0,METODA.obratna);
	simetricna(i) = diferencaFunkcije(f,1,h,x0,METODA.simetricna);
end

priblizkiSPravilom = pravilo

napakaPravilo = abs(tocnaVred - pravilo)

razlikaPraviloPrema = abs(prema - pravilo)
razlikaPraviloObratna = abs(obratna - pravilo)
razlikaPraviloSimetricna = abs(simetricna - pravilo)

function odvod = praviloZaOdv(f,h,x0)
	% Opis:
	%  diferenca izračuna odvod funkcije f v točki x0 s pomocjo nekega
	%  pravila
	%
	% Definicija:
	%  odvod = praviloZaOdv(f,h,x0)
	%
	% Vhodni podatki:
	%  f        funkcija,
	%  h        korak,
	%  x0       točka, v kateri nas zanima odvod funkcije f. 
	%
	% Izhodni podatek:
	%  odvod    priblizek odvoda f v tocki x0

	tocke = [x0 - 2*h,x0 - h,x0 + h,x0 + 2*h];
	modelVec = [1/12,-2/3,2/3,-1/12];
	dotProd = dot(f(tocke),modelVec);

	odvod = dotProd./h;
end 

function odvod = diferencaFunkcije(f,k,h,x0,metoda)
	% Opis:
	%  diferenca izračuna odvod funkcije f v točki x0 na podlagi diference, 
	%  določene z metodo
	%
	% Definicija:
	%  odvod = diferenca(f, k, metoda, h, x0)
	%
	% Vhodni podatki:
	%  f        funkcija,
	%  k        stopnja odvoda (k=1, razen pri simetrični diferenci za 2 odvod, kjer je k=2),
	%  metoda   metoda za izračun odvoda (prema, obratna ali simetrična diferenca),
	%  h        korak,
	%  x0       točka, v kateri nas zanima odvod funkcije f. 
	%
	% Izhodni podatek:
	%  odvod    odvod funkcije f v točki x0 izračunan z ustrezno metodo.
	METODA = struct("prema",0,"obratna",1,"simetricna",2);
	if(metoda == METODA.prema)
		if(k == 1)
			odvod = (f(x0 + h) - f(x0))/h;
		else
			printf("Neveljana stopna odvoda za motodo - %s",metoda);
		end
	elseif(metoda == METODA.obratna)
		if(k == 1)
			odvod = (f(x0) - f(x0 - h))/h;
		else
			printf("Neveljana stopna odvoda za motodo - %s",metoda);
		end
	elseif(metoda == METODA.simetricna)
		if(k == 1)
			odvod = (-f(x0-h) + f(x0 + h))/(2*h);
		elseif(k == 2)
			odvod = (f(x0 - h) - 2*f(x0) + f(x0 + h))/(h^2);
		else
			printf("Neveljana stopna odvoda za motodo - %s",metoda);
		end
	else
		printf("Neveljavna metoda");
		odvod = inf;
	end
end
