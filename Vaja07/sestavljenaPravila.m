format long;
m = 4;
fun = @(x) sin(x);
a = 0;
b = pi;

tocke = fun(linspace(a,b,m+1));
h = (b - a)/m; % (b - a)/m
pribTrapez = trapez(h,tocke)
pribSimpson = simpson(h,tocke)
pribShema = RombergovaShema(fun,m,a,b)

exact = integral(fun,a,b);
minMTrap = floor(fminsearch(@(m) mTrap(m,exact,1e-2,fun,a,b),100,optimset('TolX',.1)) + 1)
minMSimpson = floor(fminsearch(@(m) mSimpson(m,exact,1e-2,fun,a,b),100,optimset('TolX',.1)) + 1)

function d = mTrap(m,goal,tol,f,a,b)
	m = floor(m);
	stTock = m + 1;
	
	tocke = f(linspace(a,b,stTock));
	h = (b - a)/m;


	d = trapez(h,tocke);

	err = abs(goal - d); 
	if err < tol
		d = d + goal;
	end 
end

function d = mSimpson(m,goal,f,numPoints,a,b)
	m = floor(m);
	% m mora biti sod
	if rem(m,2) == 1
		m = m + 1
	end
	stTock = m + 1;
	
	tocke = f(linspace(a,b,stTock));
	h = (b - a)/m;

	d = simpson(h,tocke);
	
	err = abs(goal - d); 
	if err < tol
		d = d + goal;
	end 
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

function Sf = simpson(h,f)
	% Opis:
	%  Funkcija simpson izracuna priblizek za integral funkcije s
	%  sestavljenim Simpsonovim pravilom.
	%
	% Definicija:
	%  Sf = simpson(h,f)
	%
	% Vhod:
	%  h    dolzina koraka,
	%  f    vrstica vrednosti funkcije v delilnih tockah.
	%
	% Izhod:
	%  Sf   priblizek za integral funkcije, izracunan s sestavljenim
	%       Simpsonovim pravilom.
	dim = size(f,2);
	% zip 2 vectors to create alternating one
	twos = linspace(2,2,(dim-2)/2 + 1);
	fours = linspace(4,4,(dim-2)/2 + 1);
	model = [fours;twos];
	model = model(:)';
	modelVec = [1,model(1:end-1),1];
	dotProd = modelVec.*f;
	prib = sum(dotProd);
	Sf = h/3*prib;
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