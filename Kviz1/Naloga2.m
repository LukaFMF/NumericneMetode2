format long;
f = @(x) cos(pi*x);
df = @(x) -sin(pi*x)*pi;
stT = 50 + 1; % v navodilu: x = ...
interval = [0 1];

x0 = 39/61;
intM = linspace(0,1,stT);

linVred = linZlepek(intM,f(intM),x0);
absRazlikaLin = abs(linVred - f(x0))

kubVred = kubZlepek(intM,f(intM),df(intM),x0)

function zx = linZlepek(X,Y,x)
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

function zx = kubZlepek(X,Y,dY,x)
	p_s = [0 0 0 0]; 
	for i = 2:size(X,2)
		koef = trikotnaZaKubZlepek(X(i-1),X(i),Y(i-1),Y(i),dY(i-1),dY(i));

		% izpeljano iz teorije, pomagamo si s trikotnisko shemo
		vec = [
			koef(1), % *x^3
			koef(2) - 2*koef(1)*X(i-1) - koef(1)*X(i), % *x^2
			koef(3) - 2*koef(2)*X(i - 1) + koef(1)*X(i-1)^2 + 2*koef(1)*X(i-1)*X(i), % *x
			koef(4) - koef(3)*X(i-1) + koef(2)*X(i-1)^2 - koef(1)*X(i-1)^2*X(i) % *1
		]';
		p_s = [p_s;vec];
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

function koef = trikotnaZaKubZlepek(x1,x2,fx1,fx2,dfx1,dfx2)
	% x1 	fx1 	dfx1		...		...
	% x1 	fx1 	[x1,x2]f	...
	% x2 	fx2 	dfx2		
	% x2 	fx2
	koef = [0 0 dfx1 fx1];

	secondDiff = (fx2 - fx1)/(x2 - x1); % [x1,x2]f
	
	thirdDiff1 = (secondDiff - dfx1)/(x2 - x1);
	thirdDiff2 = (dfx2 - secondDiff)/(x2 - x1);
	koef(2) = thirdDiff1;

	fourthDiff = (thirdDiff2 - thirdDiff1)/(x2 - x1);
	koef(1) = fourthDiff;
end