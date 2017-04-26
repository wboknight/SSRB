function img=fanFBP(p,theta,DetWidth,StoT,StoD,N,PixSize)
% p: projection data
% theta：projection angles, e.g. 1:360, unit: degree
% DetWidth: Width of dectector
% StoT: distance from source to target ( object)
% StoD: distance from source to detector
% N: reconstructed image size
% PixSize: pixel size of reconstructed image

ND=size(p,1); %number of detector
NA=size(p,2); %projection angles
if NA~=length(theta) 
    error('inconsistency of angles!!!');
end;
DetSize=DetWidth/ND; %size of single detector size
AngInt=abs(theta(1)-theta(2))*pi/180;%angle increment

% correction
s=((0:ND-1)-0.5*ND+0.5)*DetSize;
cor=StoD./sqrt(StoD^2+s.^2);
p=p.*repmat(cor',1,NA);

% filteration
h=zeros(2*ND-1,1);
h(ND)=1/8/DetSize^2;
h(1:2:2*ND-1)=-1/2/pi^2/DetSize^2./(-ND+1:2:ND).^2;
h=[h;zeros(ND-1,1)];
H=fft(h);
p(length(H),1)=0; 
p = fft(p); 
for i = 1:NA
   p(:,i) = p(:,i).*H; 
end
p = ifft(p,'symmetric');     
p=p(ND:2*ND-1,:);   




img = zeros(N); 
sinbeta=sin(theta*pi/180+pi/2);
cosbeta=cos(theta*pi/180+pi/2);

% back projection
for j=1:N
    y=(j-0.5*N-0.5)*PixSize;
    for i=1:N
        x=(i-0.5*N-0.5)*PixSize;
        for t=1:NA
            U=1+(sinbeta(t)*x-cosbeta(t)*y)/StoT;
            c=AngInt/U^2;
            s=StoD*(cosbeta(t)*x+sinbeta(t)*y)/(StoT+sinbeta(t)*x-cosbeta(t)*y)/DetSize;            
            if (abs(s)>0.5*ND-0.5)
                continue;
            end
            ind=s+0.5*ND+0.5;
            index=floor(ind);
            lam2=ind-index;lam1=1-lam2; 
            img(N+1-j,i)=img(N+1-j,i)+(lam1*p(index,t)+lam2*p(index+1,t))*c;
        end
    end
end

img=img*StoD/StoT*DetSize;
