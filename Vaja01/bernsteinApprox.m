format long;
n = 5;
x = linspace(0,1,1001);
figure(1);
hold on;
for j = 0:n
	plot(x,bernstein(n,j,x));
end
hold off;
vrednost32 = bernstein(10,5,sqrt(3)/2)

f = @(x) 1./(3*x + 1);
y = f(x);
figure(2);
hold on;
plot(x,y,"color","red");
for m = 3:10
	Bnf = bernpoly(f,m,x);
	plot(x,Bnf);
	fprintf("Napaka(norma inf) za n = %d: %f\n",m,max(abs(Bnf - y)));
end
hold off;

 
function b = bernstein(n,i,x)
	% bernstein vrne vrednosti Bernsteinovega baznega polinoma
	%
	% b = bernstein(n,i,x)
	%
	% Vhodni podatki:
	%  n    stopnja Bernsteinovega baznega polinoma,
	%  i    indeks Bernsteinovega baznega polinoma (med 0 in n),
	%  x    seznam parametrov, v katerih racunamo vrednost polinoma.
	%
	% Izhodni podatek:
	%  b    seznam vrednosti i-tega Bernsteinovega baznega polinoma stopnje n
	%       pri parametrih iz seznama x.
	b = nchoosek(n,i)*x.^i.*(1 - x).^(n - i);
end

function Bf = bernpoly(f,n,x)
	% bernpoly vrne vrednosti Bernsteinovega polinoma stopnje n za funkcijo f v
	% tockah x.
	%
	% Bf = bernpoly(f,n,x)
	%
	% Vhod:
	%  f    funkcija: @(x) f(x),
	%  n    stopnja Bernsteinovega polinoma,
	%  x    seznam abscis.
	%
	% Izhod:
	%  Bf   seznam vrednosti Bernsteinovega polinoma stopnje n za funkcijo f v
	%       tockah iz seznama x.
	Bf = 0;
	for j = 0:n
		Bf = Bf + f(j/n).*bernstein(n,j,x);
	end
end

