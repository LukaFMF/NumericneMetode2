f = @(x,y) x + y - 1;
a = 0;
b = 1;
h = 1/4;
m = round((b - a)/h);
y0 = 1;

xs = linspace(a,b,m + 1);
ys = explicitEuler(xs,f,y0)

hold on;
axis equal;
exactY = @(x) exp(x) - x;

xx = linspace(0,1,100);
plot(xx,exactY(xx));
plot(xs,ys,".-","MarkerSize",12);
hold off;


function y = explicitEuler(x,f,y0)
	% Opis:
	%  Funkcija expeuler vrne numericno resitev navadne diferencialne enacbe
	%  y' = f(x,y) pri pogoju y(x(1)) = y0, ki je izracunana z eksplicitno
	%  Eulerjevo metodo.
	%
	% Definicija:
	%  y = expeuler(x,f,y0)
	%
	% Vhodni podatki:
	%  x    vrstica delilnih tock,
	%  f    funkcija f v obliki @(x,y) f(x,y),
	%  y0   zacetna vrednost resitve.
	%
	% Izhodni podatek:
	%  y    vrstica numericnih priblizkov za vrednosti tocne resitve v delilnih
	%       tockah.
	dim = size(x,2);
	y = zeros(1,dim);
	y(1) = y0;

	for i = 2:dim
		y(i) = y(i-1) + (x(i) - x(i-1))*f(x(i-1),y(i-1));
	end
end