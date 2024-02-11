%% CG method for solving linear system
% code copy form internet; 
% code adress:https://blog.csdn.net/lusongno1/article/details/68942821
function [x,steps] = CG_method(A,b,x0,eps)
r0 = b - A*x0;
p0 = r0;
if nargin == 3
    eps = 1.0e-6;
end
steps = 0;
while 1
    if abs(p0) < eps
        break;
    end
    steps = steps + 1;
    a0 = r0'*r0/(p0'*A*p0);%多次用到可以存一步。
    x1 = x0 + a0*p0;

    r1 = r0 -a0*A*p0;

    b0 = r1'*r1/(r0'*r0);
    %这里的r'* r虽然后面可能还会用到，但是由于计算量不大，没有必要再设个新变量将
    %其存下了，内存上的开销划不来。
    p1 = r1 + b0*p0;

    %只是用到前后两层的向量，所以节省内存开销，计算完后面一层，可以往回覆盖掉没用
    %的变量。

    x0 = x1;
    r0 = r1;
    p0 = p1;

end
x = x0;
end
% ————————————————
% 版权声明：本文为CSDN博主「lsec小陆」的原创文章，遵循 CC 4.0 BY-SA 版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/lusongno1/article/details/68942821