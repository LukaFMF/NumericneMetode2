format long;
f = @(x) sin(2*pi*x);
df = @(x) cos(2*pi*x)*2*pi;
stT = 10^5 + 1; % v navodilu: x = ...
interval = [0 1];

kubNapaka10 = napakaKubZlepek(f,df,10,stT,interval)
kubNapaka20 = napakaKubZlepek(f,df,20,stT,interval)
kubNapaka50 = napakaKubZlepek(f,df,50,stT,interval)

function err = napakaKubZlepek(f,df,m,stTock,interval)
	intM = linspace(interval(1),interval(end),m + 1);
	
	x = linspace(interval(1),interval(end),stTock);
	vredF = f(x);
	err = max(abs(vredF - kubZlepek(intM,f(intM),df(intM),x)));
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