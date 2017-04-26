readrawdata3; %Read in data

%%
%Rebinning
sz_or = -1.2;
dz = 0.3/360;
D = 8;
R = 4;

prj = zeros(256,360,256);

for theta = 1:size(prj,2)
    for z = 1:size(prj,3)
        oz = -1+2/255*(z-1); %z of object
        n = round((oz-sz_or-theta*dz)/0.3);
        sz = sz_or + n*0.3+theta*dz; %z of source
        Dz = oz - sz;
        for x = 1:size(prj,1)
            u = (x-128)/256.0*4.0;
            v = (u*u+D*D)/D/R*Dz;
            prj(x,theta,z) = sqrt((u*u+D*D)/(u*u+v*v+D*D))*...
                p(x,ceil(v/0.01+64),n*360+theta);
        end
    end
end

%%
%Reconstruction
img = zeros(256,256,256);
for z = 1:size(img,3)
    z
    im=fanFBP(prj(:,:,z),1:360,4,4,8,256,4/512);
    img(:,:,z)=im;
end
%To make reconstructed image same as phantom
im = img;
for i = 1:size(img,3)
    im(size(img,3)-i+1,end:-1:1,:) = img(:,:,i);
end
