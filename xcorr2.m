function xc = xcorr2(x,l)

xc = zeros(2*length(x),1);

y = zeros(2*length(x),1);
offset = floor(length(x)/2);
y(offset:offset+length(x)-1) = x;

lag = 1;
for i=-(l-1):l-1
    sum = 0;
    for k=0:l-1
        sum = sum + ( y(k+offset) * conj(y(k+i+offset)) );
    end
    xc(lag) = sum;
    lag = lag + 1;
end

end