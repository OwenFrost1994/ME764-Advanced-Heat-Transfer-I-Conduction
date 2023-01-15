function [x,y]=geometry(bs,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
W=1.5e-2; %width of rectangle
H=1.5e-2; %height of rectangle
a=0.5e-2; %minor axis half-width
b=1e-2; %major axis half-width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbs=10; %number of boundary segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0,
 x=nbs; % number of boundary segments
 return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameterized boundaries
d=[
 0 a/H 0 1-b/W 0 a/H 0 1-b/W 0 0 % start parameter value
 a/H 1 1-b/W 1 a/H 1 1-b/W 1 1 1 % end parameter value
 0 0 0 0 0 0 0 0 2 2% left hand region
 2 1 1 2 2 1 1 2 1 1% right hand region
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs1=bs(:)';
if find(bs1<1 | bs1>nbs),
 error('PDE:rectangle:InvalidBs', 'Non existent boundary segment number.')
end
if nargin==1,
 x=d(:,bs1);
 return
end
x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
 bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
 error('PDE:rectangle:SizeBs', 'bs must be scalar or of same size as s.');
end
if ~isempty(s),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%boundary segment coordinates
 % boundary segment 1
 ii=find(bs==1);
 if length(ii)
 x(ii)=0;
 y(ii)=s(ii)*H;
 end

 ii=find(bs==2);
 if length(ii)
 x(ii)=0;
 y(ii)=s(ii)*H;
 end
 % boundary segment 2
 
 ii=find(bs==3);
 if length(ii)
 x(ii)=s(ii)*W;
 y(ii)=H;
 end

 ii=find(bs==4);
 if length(ii)
 x(ii)=s(ii)*W;
 y(ii)=H;
 end
 % boundary segment 3
 ii=find(bs==5);
 if length(ii)
 x(ii)=W;
 y(ii)=(1-s(ii))*H;
 end

 ii=find(bs==6);
 if length(ii)
 x(ii)=W;
 y(ii)=(1-s(ii))*H;
 end
 % boundary segment 4
 ii=find(bs==7);
 if length(ii)
 x(ii)=(1-s(ii))*W;
 y(ii)=0;
 end

 ii=find(bs==8);
 if length(ii)
 x(ii)=(1-s(ii))*W;
 y(ii)=0;
 end
 % boundary segment 5
 ii=find(bs==9);
 if length(ii)
 x(ii)=b*cos(s(ii)*pi/2);
 y(ii)=a*sin(s(ii)*pi/2);
 end

 % boundary segment 6
 ii=find(bs==10);
 if length(ii)
 x(ii)=W+b*cos(pi+s(ii)*pi/2);
 y(ii)=H+a*sin(pi+s(ii)*pi/2);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end