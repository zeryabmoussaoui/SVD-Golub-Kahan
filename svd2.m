%svd2 - returns the singular values of matrix A in descending order using Kahan-Golub method
% Reference : Golub-Kahan-1965 Calculating the singular values and pseudo-inverse of a matrix
%
% Syntax:  [sv] = svd2(A)
%
% Inputs:
%    A - n-by-m - input matrix
%
% Outputs:
%    sv - 1-by-min(n,m)-  singulars values in 
%
% Example: 
%    A=[3 2;1 4;0.5 0.5];
%    s=svd2(A);
%    s=
%
% Other m-files required: none
% Subfunctions: sturm
% MAT-files required: none
%
% See also: svd, gsvd
% Author: Zeryab Moussaoui
% email: zeryab.moussaoui@gmail.com
% Website: www.github.com/zeryabmoussaoui
% May 2016; Last revision: 12-June-2019

%------------- BEGIN CODE --------------
function [sv] = svd2(A)
 
%%% 1 - Transform into bi-diagonal matrix 
 
[m , n ] = size (A);

%% 1.1 - Transpose to get an along rectangular matrix 
if (m < n)
    A=A';
    tmp=m;
    m=n;
    n=tmp;
end
 
B=A;

%% 1.2 - Householder Transformation
% TODO : Case null element in diagional   
for k=1:n
 d=norm(B(k:m,k));
 if d~=0
    vkk=B(k,k)+d*sign(B(k,k));      
    h=d*(d+abs(B(k,k)));
    v=[zeros(k-1,1);vkk;B(k+1:m,k)];
    S= eye(m) - v*v'/ h;
    B=S*B;
 end
 if (k<n)
     r=norm(B(k,k+1:n));
     if r ~= 0
       vkk=B(k,k+1)+r*sign(B(k,k+1)); % v(k,k+1)
       h=r*(r+abs(B(k,k+1)));
       v=[zeros(k,1);vkk;B(k,k+2:n)'];
       P=eye(n) - v*v'/h;
       B=B*P;
     end
 end
end
 
 
%%% 2 - Transform into square matrix
 
[m , n ] = size (B);
if (m >= n)
    C=B(1:n,1:n);
    l=n;
else
    C=[ B(1:m,1:m+1) ; zeros(1,m+1)];
    l=m+1;
end
 
 
%%% 3 - Transform into real tridiagonal matrix
s=abs(diag(C));
t=abs(diag(C(1:l,2:l)));
 
K=zeros(l,l);
% i=1
K(1,1)= s(1)^2;
K(2,1)= s(1)*t(1);
for i=2:l-1
    K(i,i)   = s(i)^2 + t(i-1)^2 ;
    K(i-1,i) = s(i-1)*t(i-1);
    K(i+1,i) = s(i)*t(i);
    
end
% i=l
K(l,l) =  s(l)^2 + t(l-1)^2;
K(l-1,l) = s(l-1)*t(l-1);
 
 
%%% 4 - Compute eigenvalues of K using Strum method
 
%% 4.0 - Define Sturm function to compute sign change

function nb = strum(M,x)
 
% M : square matrix
% x : evaluation point
% return : nb of sign change
 
%TODO : Use the fraction to avoid l'overflow
 
n = length(M) ; % TODO : non square case ; test if tridiagonal
nb=0;
p=ones(n+1,1);
s=ones(n+1,1);
p(1)=1;
s(1)=1;
for i=2:n+1
   if (i==2)
       p(i) = M(1,1) - x;
    else
        d=M(i-1,i-1);
        e=M(i-1,i-2);
        p(i) = (d - x)*p(i-1) - (e^2)*p(i-2);
   end
   if (p(i)==0)
       s(i)=s(i-1);
   else
       s(i)=sign(p(i));
   end
   
    if ( s(i)*s(i-1) == -1)
        nb=nb+1;
    end
end
 
%% 4.1 - Look for an [a,b] interval with all eingenvalues
 
% Define diagonals 
al=(diag(K));
bet=abs(diag(K(1:l,2:l))); % take absolute value for Disk formula
 
% Gerschgorin's disks
 
I=zeros(l,2); % array with intervals
 
I(1,:)= [ -bet(1)+al(1)  bet(1)+al(1)];
for i=2:l-1
    I(i,:)=[ -bet(i)-bet(i-1)-al(i) bet(i)+bet(i-1)+al(i) ];
end
I(l,:)= [-bet(l-1)+al(l)  bet(l-1)+al(l) ];
 
% Gerschgorin's interval bounds
a=min(I(:,1));
b=max(I(:,2));
L=[a b];

%% 4.2 Bisectrix method 
 
R=zeros(l,2);  % Interval with only one eigenvalue
N=0;
while ( N ~= l )
    N=0;
    R=zeros(l,2);
    for i=2:length(L)
        neig=strum(K,L(i)) - strum(K,L(i-1));
        if(neig > 1)
            med=(L(i)+L(i-1))/2;
            L=[L med];
        elseif( neig ==1)
            N=N+1;
            R(N,:)=[L(i-1) L(i)];
        else
           %do nothing
        end 
    end
    L=sort(L);
end

end 
% Refine by bissection
 
eps=10^(-2);
 
   for i=1:l
       delta=eps*2; % To get around from do-while
        while(delta > eps)
            lmax=R(i,2);
            lmin=R(i,1);
            lmid=(lmax+lmin)/2;
            neig=strum(K,lmax) - strum(K,lmid);
            if ( neig == 1)
               R(i,1)=lmid; 
            else % no eigenvalue
               R(i,2)=lmid; 
            end 
         delta=norm( (R(i,2)-R(i,1))./( R(i,2)+R(i,1) ) );
        end
   end
   eigK=(R(:,1)+R(:,2))/2; % Eigenvalues of K
   sv=sort(sqrt(eigK));
   
end

%------------- END OF CODE --------------
%Please send suggestions for improvement of the above code 
%to Zeryab Moussaoui at this email address: zeryab.moussaoui@gmail.com
%Your contribution towards improving this file will be acknowledged in
%the "Changes" section of the TEMPLATE_HEADER web page on the Matlab
%Central File Exchange
