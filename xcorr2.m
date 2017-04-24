function xc = xcorr(x,l)

xc = zeros(2*length(x),1);

y = zeros(2*length(x),1);
offset = floor(length(x)/2);
y(offset:offset+length(x)-1) = x;

for i=0:l-1
    sum = 0;
    for k=0:l-1
        sum = sum + ( y(k+offset) * conj(y(k+i+offset)) );
    end
    xc(i+1) = sum;
end

end