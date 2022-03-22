
dx = diferencaFunkcije(@(x) sin(x),1,2,0,"simetricna")

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
	if(metoda == "prema")
		if(k == 1)
			odvod = (f(x0 + h) - f(x0))/h;
		else
			printf("Neveljana stopna odvoda za motodo - %s",metoda);
		end
	elseif(metoda == "obratna")
		if(k == 1)
			odvod = (f(x0) - f(x0 - h))/h;
		else
			printf("Neveljana stopna odvoda za motodo - %s",metoda);
		end
	elseif(metoda == "simetricna")
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

