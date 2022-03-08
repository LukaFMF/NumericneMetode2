f = @(x) sin(2*pi*x);
minM = floor(fminsearch(@(m) diffFromNorm(m,1e-6,f,100001),50,optimset('TolX',.1)) + 1)

function d = diffFromNorm(m,goal,f,numPoints)
	m = floor(m)
	intM = linspace(0,1,m + 1);

	
	x = linspace(0,1,numPoints);
	vredSin = f(x);
	err = max(abs(vredSin - linzlepek(intM,f(intM),x)));
	d = err - goal;
	if err < goal
		d = d + goal;
	end 
end

function zx = linzlepek(X,Y,x)
	% Opis:
	%  Funkcija linzlepek vrne vrednosti zveznega linearnega zlepka, ki
	%  interpolira dane vrednosti v delilnih tockah definicijskega intervala.
	%
	% Definicija:
	%  zx = linzlepek(X,Y,x)
	%
	% Vhod:
	%  X    vrstica, ki predstavlja delilne tocke definicijskega intervala
	%       zlepka,
	%  Y    interpolacijske vrednosti,
	%  x    vrstica tock na definicijskem intervalu.
	%
	% Izhod:
	%  zx   vrstica vrednosti linearnega zlepka v tockah iz x.
	p_s = [0 0]; 
	for i = 2:size(X,2)
		clen = (Y(i) - Y(i - 1))/(X(i) - X(i - 1));
		p_s = [p_s; clen, Y(i - 1) - clen*X(i - 1)];
	end
	p_s = p_s(2:end,:);

	zx = zeros(size(x));
	for k = 1:size(x,2)
		vred = x(k);
		for i = 1:(size(X,2) - 1)
			meja = X(i + 1);
			if vred <= meja
				zx(k) = polyval(p_s(i,:),vred);
				break
			end
		end
	end
end