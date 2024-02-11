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
    a0 = r0'*r0/(p0'*A*p0);%����õ����Դ�һ����
    x1 = x0 + a0*p0;

    r1 = r0 -a0*A*p0;

    b0 = r1'*r1/(r0'*r0);
    %�����r'* r��Ȼ������ܻ����õ����������ڼ���������û�б�Ҫ������±�����
    %������ˣ��ڴ��ϵĿ�����������
    p1 = r1 + b0*p0;

    %ֻ���õ�ǰ����������������Խ�ʡ�ڴ濪�������������һ�㣬�������ظ��ǵ�û��
    %�ı�����

    x0 = x1;
    r0 = r1;
    p0 = p1;

end
x = x0;
end
% ��������������������������������
% ��Ȩ����������ΪCSDN������lsecС½����ԭ�����£���ѭ CC 4.0 BY-SA ��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
% ԭ�����ӣ�https://blog.csdn.net/lusongno1/article/details/68942821