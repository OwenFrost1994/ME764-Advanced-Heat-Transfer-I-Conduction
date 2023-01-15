function [x,y]=HW14geometry6(bs,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
Rh=0.02;
W=0.01;
Rt=0.005;
th=0.006;
L=0.025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbs=10;  %number of boundary segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameterized boundaries
d=[
  0 0 0 0 0 0 0 0 0 0% start parameter value
  1 1 1 1 1 1 1 1 1 1% end parameter value
  0 0 0 0 0 0 1 1 1 1% left hand region
  1 1 1 1 1 1 2 2 2 2% right hand region
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
    %% boundary segments for outer circle
    ii=find(bs==1);
    if length(ii)
        x(ii)=Rh*cos(s(ii)*pi);
        y(ii)=Rh*sin(s(ii)*pi);
    end

    % outer circle segment 2
    ii=find(bs==2);
    if length(ii)
        x(ii)=Rh*cos(pi+s(ii)*pi);
        y(ii)=Rh*sin(pi+s(ii)*pi);
    end
    
    %% boundary segments for left circle
    ii=find(bs==3);
    if length(ii)
        x(ii)=-W+Rt*cos(s(ii)*pi);
        y(ii)=Rt*sin(s(ii)*pi);
    end

    % left circle boundary segment 2
    ii=find(bs==4);
    if length(ii)
        x(ii)=-W+Rt*cos(pi+s(ii)*pi);
        y(ii)=Rt*sin(pi+s(ii)*pi);
    end
    
    %% boundary segments for right circle
    ii=find(bs==5);
    if length(ii)
        x(ii)=W+Rt*cos(s(ii)*pi);
        y(ii)=Rt*sin(s(ii)*pi);
    end
    
    % right circle boundary segment 2
    ii=find(bs==6);
    if length(ii)
        x(ii)=W+Rt*cos(pi+s(ii)*pi);
        y(ii)=Rt*sin(pi+s(ii)*pi);
    end
    
    %% boundary segments for rectangle
    ii=find(bs==7);
    if length(ii)
        x(ii)=-th/2;
        y(ii)=-L/2+s(ii)*L;
    end
    
    ii=find(bs==8);
    if length(ii)
        x(ii)=-th/2+s(ii)*th;
        y(ii)=L/2;
    end
    
    ii=find(bs==9);
    if length(ii)
        x(ii)=th/2;
        y(ii)=L/2-s(ii)*L;
    end
    
    ii=find(bs==10);
    if length(ii)
        x(ii)=th/2-s(ii)*th;
        y(ii)=-L/2;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end