format long;
f = @(x) sin(x);

pribPoStirih = RombergovaShema(f,4,0,pi)

pribPoDveh = RombergovaShema(f,2,0,pi);
tocnaVred = integral(f,0,pi);
napakaPoDveh = abs(pribPoDveh - tocnaVred)

function Rf = RombergovaShema(fun,k,a,b)
	shema = zeros(1,k + 1);

	% poracunamo traprezne priblizke
	stInt = 1;
	for i = 1:k + 1
		h = (b - a)/stInt;
		tocke = fun(linspace(a,b,stInt+1));
		shema(i) = trapez(h,tocke);
		stInt = stInt * 2;
	end

	nova = zeros(1,k + 1);
	for i = 1:k
		koef = 2^(2*i);
		for j = i+1:k + 1
			nova(j) = (koef*shema(j) - shema(j - 1))/(koef - 1);
		end
		shema = nova;
	end

	Rf = shema(end);
end

function Tf = trapez(h,f)
	% Opis:
	%  Funkcija trapez izracuna priblizek za integral funkcije s sestavljenim
	%  trapeznim pravilom.
	%
	% Definicija:
	%  Tf = trapez(h,f)
	%
	% Vhod:
	%  h    dolzina koraka,
	%  f    vrstica vrednosti funkcije v delilnih tockah.
	%
	% Izhod:
	%  Tf   priblizek za integral funkcije, izracunan s sestavljenim trapeznim
	%       pravilom.
	dim = size(f,2);
	modelVec = [1,linspace(2,2,dim-2),1];
	dotProd = modelVec.*f;
	prib = sum(dotProd);
	Tf = h/2*prib;
end