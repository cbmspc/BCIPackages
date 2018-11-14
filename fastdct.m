% Fast discrete cosine transform
% based on https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms

function [out, N] = fastdct (in)
if size(in,1) == 1
    in = in.';
end
for ch = size(in,2):-1:1
    [out(:,ch), N] = fastdct_2(in(:,ch));
end

function [out, N] = fastdct_2 (in)
in = in(:);
N = length(in);
N8 = ceil(N/8)*8;
if N < N8
    %in = [in; flipud(in(end - (N8-N) + 1 : end))];
    in = [in; zeros(N8-N,1)];
end
out = in*0;
for i = 0:8:length(in)-8
    out(i+(1:8)) = fastdct_1(in(i+(1:8)));
end


function v = fastdct_1 (v)
S = [
        0.353553390593273762200422
        0.254897789552079584470970
        0.270598050073098492199862
        0.300672443467522640271861
        0.353553390593273762200422
        0.449988111568207852319255
        0.653281482438188263928322
        1.281457723870753089398043
];

A = [
        NaN
        0.707106781186547524400844
        0.541196100146196984399723
        0.707106781186547524400844
        1.306562964876376527856643
        0.382683432365089771728460
];


v0 = v(1) + v(8);
v1 = v(2) + v(7);
v2 = v(3) + v(6);
v3 = v(4) + v(5);
v4 = v(4) - v(5);
v5 = v(3) - v(6);
v6 = v(2) - v(7);
v7 = v(1) - v(8);

v8 = v0 + v3;
v9 = v1 + v2;
v10 = v1 - v2;
v11 = v0 - v3;
v12 = -v4 - v5;
v13 = (v5 + v6) * A(4);
v14 = v6 + v7;
v15 = v8 + v9;
v16 = v8 - v9;
v17 = (v10 + v11) * A(2);
v18 = (v12 + v14) * A(6);
v19 = -v12 * A(3) - v18;
v20 = v14 * A(5) - v18;
v21 = v17 + v11;
v22 = v11 - v17;
v23 = v13 + v7;
v24 = v7 - v13;
v25 = v19 + v24;
v26 = v23 + v20;
v27 = v23 - v20;
v28 = v24 - v19;
v(1) = S(1) * v15;
v(2) = S(2) * v26;
v(3) = S(3) * v21;
v(4) = S(4) * v28;
v(5) = S(5) * v16;
v(6) = S(6) * v25;
v(7) = S(7) * v22;
v(8) = S(8) * v27;
