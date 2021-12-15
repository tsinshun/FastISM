%% Efficient ISM image reconstruction
%  Author: shun qin
%  Contact: shun.qin@outlook.com
%
% imstack: raw image of real sample, F: index for each frame (from 1 to N),
% (X,Y) is the coordinate of each scanning position
% scale: image scale factor, not necessary but still available to use
% Npixel: the window size (Npixel*Npixel) used to cut each light-spot 
% Bg: background 

% ref: [1] Shun Qin,Fast image scanning microscopy with efficient image
% reconstruction,Optics and Laser in Engineering,2022.
%% 
function [Iism,Ipat] = getISM(imstack, F, X, Y,scale, Npixel,Bg)
% Npixel: Width of PSF

[M,N,Lz] = size(imstack);
% imstack(imstack-Bg<0)=0;


Nframe = max(F(:));
Iism = 0;
% Isup = zeros(M*scale*2,N*scale*2,Lz);
Ipat = 0;
for i = 1:Nframe  
    
    if i==37
       1; 
    end

   idx = find(F==i);
%    F1 = F(idx);
   X1 = X(idx);
   Y1 = Y(idx);
   
   frame = imstack(:,:,i);
   frame = frame - Bg; frame(frame<0) = 0; % remove background
   if scale>1
      I = imresize(frame,scale,'Method','lanczos2'); % lanczos2 bicubic
   else
       I = frame;
   end
%    I = frame-frame + 100;
%    I = I/max(I(:));
   [Iism0,Iism1] = copyPSF2(I,Npixel*scale,[X1 Y1]*scale,scale);
%    Isup(:,:,i) = Iism0;
%    Ipat(:,:,i) = Iism1;
   Iism = Iism + Iism0;
   Ipat = Ipat + Iism1;
   fprintf('Frame number = %d \n',i);
   imagesc(Iism0);title(['nFrame: ' num2str(i)]); colormap hot
   axis image;drawnow
   
end
if scale >1
    Iism = imresize(Iism,2*[M N]);
end
% Iism = sum(Isup,3);
% Ipat = sum(Ipat,3);

% imshow(mat2gray(Iism))
% 
% Iavg = mean(imstack1,3);
% Iavg = imresize(Iavg,2/scale,'Method','lanczos2');





function [Iism, Iism1]= copyPSF2(I,Npixel,Position,Scale)

[M, N] = size(I);
Np = floor(Npixel/2);
% Inew = zeros(M+Np,N+Np);
Iism = zeros(2*M,2*N);
Iism1 = zeros(2*M,2*N);

Inew = padarray(I, [Np Np]);
Iism = padarray(Iism, [Np Np]);  %Iism = Iism + nan;
Iism1 = padarray(Iism1, [Np Np]);

X = Position(:,1);
Y = Position(:,2);
 
L = length(X);

% PSF = gaussian2d(2*Np,3); PSF = PSF./max(max(PSF));%imagesc((PSF));

for i=1:L
    n = round(X(i)+0.5); m = round(Y(i)+0.5);
    x0 = n + Np;  % plus 0.5 to get real pixel index,i.e. (m,n)
    y0 = m + Np; 
    x1 = 2*n + Np;  
    y1 = 2*m + Np;
    V1 = 0.5+[X(i)-n  Y(i)-m];
    
    if x0>0 && y0>0 && x1>0 && y1>0
        if (x1-Np<1 || x1+Np>2*N+2*Np|| y1-Np<1|| y1+Np>2*M+2*Np) || (x0-Np<1 || x0+Np>2*N+2*Np|| y0-Np<1|| y0+Np>2*M+2*Np)
%             disp('window for PSF copy out of idx...');
            continue;
        end
        
%         if y0-Np>0 && y0-Np<=M && x0-Np>0 && x0-Np<=N && G(y0-Np,x0-Np)<=1            
%         end    
         %Iism(y1-Np:y1+Np,x1-Np:x1+Np) =  PSF;
%         W = fspecial('gaussian',2*Np+1,1.09);
%         W = ones(2*Np+1,2*Np+1);

%         [x,y] = meshgrid(x0-Np:x0+Np,y0-Np:y0+Np);
%         sigma = 1.09;
%         r3 = (x-(X(i)+0.5 + Np)).^2 + (y-(Y(i)+0.5 + Np)).^2;
%         W = exp(-r3/(2*sigma^2));

        W = Inew(y0-Np:y0+Np,x0-Np:x0+Np);
        [x,y] = meshgrid(x0-Np:x0+Np,y0-Np:y0+Np);
        r2 = (x-x0).^2 + (y-y0).^2;
        Id = r2>Np^2;
        W(Id)=0;               
                
        W1 = W; 
        W1= biShiftPSF(W,V1); % subpixel shifting
        Id = Id | (W1<0);     
        W1(Id)=0;     
        Iism(y1-Np:y1+Np,x1-Np:x1+Np) = W1;
        
%         V2 = [x1 y1] - (2*[X(i)+0.5 Y(i)+0.5]+Np);
%         Iism1(y1-Np:y1+Np,x1-Np:x1+Np) =  ShiftPSF(getPSF1(2*Np+1, 1.09),-V2);
        
        [x,y] = meshgrid(x1-Np:x1+Np,y1-Np:y1+Np);
        sigma = 1.09*Scale;
        r3 = (x-(2*(X(i)+0.5) + Np)).^2 + (y-(2*(Y(i)+0.5) + Np)).^2;
        W = exp(-r3/(2*sigma^2));       
        Iism1(y1-Np:y1+Np,x1-Np:x1+Np) = W;
    end

end
Iism = Iism(Np+1:end-Np,Np+1:end-Np);
Iism1 = Iism1(Np+1:end-Np,Np+1:end-Np);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1,3) (2,3)
% (1,4) (2,4)
% f(i+u,j+v),u=0.2,v=0.4, i=1, j=3
% x direction: f(R1)=u(f(Q21)-f(Q11))+f(Q11)
% y direction:
% or direct way: 
% f(i+u,j+v) = (1-u)(1-v)f(i,j) + (1-u)vf(i,j+1) + u(1-v)f(i+1,j) + uvf(i+1,j+1)

%%%%%%%%%%%%%%%%% shift light-spot by bilinear interpolation %%%%%%%%%%%%%%%%%%%%%%%
function Ish= biShiftPSF(W,V)

% W1 = padarray(W0,[10 10]);
% W = circshift(W1,[5 10]); 
[M,N] = size(W);
[X,Y] = meshgrid(1:N,1:M);
X=X-V(1);
Y=Y-V(2);

X1 = floor(X); % j
Y1 = floor(Y); % i
Dx=X-X1; % v
Dy=Y-Y1; % u


% (X_ISM,Y_ISM)=(X,Y)-V
% i.e. (Dx,Dy) = (X_ISM,Y_ISM)-(floor(X_ISM),floor(Y_ISM))
if V(1)>0
    Dx = 1-V(1);
else
    Dx = -V(1);   % v is dx
end

if V(2)>0
    Dy = 1-V(2);
else
    Dy =-V(2);  % u is dy
end

Cy = 1-Dy; % 1-dy
Cx = 1-Dx; % 1-dx

C1 = Cy.*Cx; % (1-dy)(1-dx)
C2 = Cy.*Dx; % (1-dy)dx
C3 = Cx.*Dy; % (1-dx)dy
C4 = Dx.*Dy; %  dxdy

% F1([Y1 X1]) = W([Y1 X1]);
% W1 = padArray(W,1);
Ish = zeros(M,N);
for m=1:M
    for n=1:N
        i = Y1(m,n);
        j = X1(m,n);
        if i>0 && j>0 && i+1<=M && j+1<=N
            Ish(m,n) = C1*W(i,j) + C2*W(i,j+1) + C3*W(i+1,j) + C4*W(i+1,j+1);
        end
        
    end
end
1;
% imagesc(Ish)
