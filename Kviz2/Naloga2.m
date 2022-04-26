format long;
a = 243/37;
f = @(x) 3*exp(a*x);

p3 = RombergovaShema(f,3,0,1);
p4 = RombergovaShema(f,4,0,1);

absRazlikaPrib = abs(p3 - p4)


tocnaVred = integral(f,0,1);
r = 0;
while(true)
	h = 1/2^r;
	Tf = trapez(h,f(0:h:1));
	if abs(Tf - tocnaVred) < abs(p3 - tocnaVred)
		pribPrikriteriju = Tf
		break
	end
	r = r + 1;
end

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