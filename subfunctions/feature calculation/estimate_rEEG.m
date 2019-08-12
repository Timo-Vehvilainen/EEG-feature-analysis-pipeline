function [t, reeg,reeg_nl,aref] = estimate_rEEG(nd1, fs2);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

epl = fs2*2; N = length(nd1);
block_no = floor(N/epl);
for ii = 1:block_no;
    q1 = (ii-1)*epl+1; q2 = q1+epl-1;
    reeg(ii) = max(nd1(q1:q2))-min(nd1(q1:q2));
end
t = 1:2:block_no*2;
reeg_nl = reeg;
reeg_nl(1,reeg(1,:)>10) = 10*log10(reeg(1,reeg(1,:)>10));

ef1 = [0 5 10 20 50 100 200 500];
aref(1,:) = ef1; aref(2,:) = ef1;
aref(2, ef1>10) = 10*log10(ef1(ef1>10));

end

