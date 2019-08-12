function z2 = nlin_map(tp1, tp2, s1, s2, s3, z1);

z1 = abs(z1);
ref1 = find(z1<tp1);
ref2 = find(z1>=tp1 & z1<tp2);
ref3 = find(z1>=tp2);
ref4 = find(z1>=100);

z2(ref1) = s1*z1(ref1);
z2(ref2) = s2*(z1(ref2)-tp1)+(s1*tp1);
z2(ref3) = s3*(z1(ref3)-tp2)+s2*(tp2-tp1)+(s1*tp1);
z2(ref4) = s3*(100-tp2)+s2*(tp2-tp1)+(s1*tp1);
