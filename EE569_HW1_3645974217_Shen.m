%EE569 HOMEWORK ASSGINMENT #1
%Date: 2020/01/28
%Name: Shen Zhihong
%ID: 3645974217
%email: shenzhih@usc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem1: Image demosaicing and histogram manipulation
%Implementation1: BilinearDemosaicing
%Outputimage:demosaicing_a.raw
%Implementation2:MHC Demosaicing
%OutputImage:demosaicing_b.raw
%(c)histogram
%Implementation: transfer-function-based histogram manipulation
%OutputImage:histogram_c.raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem2:Image Denoising
%Implementation:Basic Denoising methods(uniform)
%OutputImage:basic_uniform.raw
%Implementation:Basic Denoising methods(Gaussian)
%OutputImage:basic_gaussain.raw
%Implementation:Bilateral Filtering
%OutputImage: bilateral.raw
%Implementation:NLM Filtering
%OutputImage:nlm_Filtering.raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Problem1(a)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row = 600;
col = 532;
G  = readraw('Dog.raw',600,532,1);
G_o = readraw('Dog_ori.raw',600,532,3);
maxvalue = max(max(G_o));
G_padded = zeros(col+2,row+2);
G_padded(2:end-1,2:end-1) = G;
G_padded(2:end-1,1) = G_padded(2:end-1,3);
G_padded(2:end-1,end) = G_padded(2:end-1,end-2);
G_padded(1,:) = G_padded(3,:);
G_padded(end,:) = G_padded(end-2,:);
imshow(G/255);
figure(1);
imshow(G_o/255);
G_r = zeros(col,row);
G_g = zeros(col,row);
G_b = zeros(col,row);
% 
for i = 1:2:col
    for j = 1:2:row
        G_r(i,j) = 0.5*(G_padded(i+1,j)+G_padded(i+1,j+2));
        G_r(i,j+1) = G_padded(i+1,j+2);
        G_r(i+1,j+1) = 0.5*(G_padded(i+1,j+2)+G_padded(i+3,j+2));
        G_r(i+1,j) = 0.25*(G_padded(i+1,j)+G_padded(i+1,j+2)+G_padded(i+3,j+2)+G_padded(i+3,j));
    end
end
% G_r(1:2:end,1:2:end) = 0.5*(G_padded(2:2:end-2,1:2:end-3)+G_padded(2:2:end-2,3:2:end-1));
% G_r(1:2:end,2:2:end) = G_padded(2:2:end-2,3:2:end-1);
% G_r(2:2:end,2:2:end) = 0.5*(G_padded(2:2:end-2,3:2:end-1)+G_padded(4:2:end,3:2:end-1));
% G_r(2:2:end,1:2:end) = 0.25*(G_padded(2:2:end-2,1:2:end-3)+G_padded(2:2:end-2,3:2:end-1)+G_padded(4:2:end,3:2:end-1)+G_padded(4:2:end,1:2:end-3));
%imshow(G_r/255);


for i =1:2:col
    for j = 1:2:row
        G_b(i,j) = 0.5*(G_padded(i,j+1)+G_padded(i+2,j+1));
        G_b(i+1,j) = G_padded(i+2,j+1);
        G_b(i,j+1) = 0.25*(G_padded(i,j+1)+G_padded(i+2,j+1)+G_padded(i,j+3)+G_padded(i+2,j+3));
        G_b(i+1,j+1) = 0.5*(G_padded(i+2,j+1)+G_padded(i+2,j+3));
    end
end

% G_b(1:2:end,1:2:end) = 0.5*(G_padded(1:2:end-3,2:2:end-2)+G_padded(3:2:end-1,2:2:end-2));
% G_b(2:2:end,1:2:end) = G_padded(3:2:end-1,2:2:end-2);
% G_b(1:2:end,2:2:end) = 0.25*(G_padded(1:2:end-3,2:2:end-1)+G_padded(3:2:end-1,2:2:end-2)+G_padded(1:2:end-3,4:2:end)+G_padded(3:2:end-1,4:2:end));
% G_b(2:2:end,2:2:end) = 0.5*(G_padded(3:2:end-1,2:2:end-2)+G_padded(3:2:end-1,4:2:end));
% 
for i =1:2:col
    for j =1:2:row
        G_g(i,j) = G_padded(i+1,j+1);
        G_g(i+1,j+1) = G_padded(i+2,j+2);
        G_g(i,j+1) = 0.25*(G_padded(i+2,j+2)+G_padded(i+1,j+1)+G_padded(i+1,j+3)+G_padded(i,j+1));
        G_g(i+1,j) = 0.25*(G_padded(i+2,j+2)+G_padded(i+1,j+1)+G_padded(i+3,j+1)+G_padded(i+2,j));
    end
end

% % G_g(1:2:end,1:2:end) = G_padded(2:2:end-2,2:2:end-2);
% % G_g(2:2:end,2:2:end) = G_padded(3:2:end-1,3:2:end-1);
% % G_g(1:2:end,2:2:end) = 0.25*(G_padded(3:2:end-1,3:2:end-1)+G_padded(2:2:end-2,2:2:end-2)+G_padded(2:2:end-2,4:2:end)+G_padded(1:2:end-3,2:2:end-2));
% % G_g(2:2:end,1:2:end) = 0.25*(G_padded(3:2:end-1,3:2:end-1)+G_padded(2:2:end-2,2:2:end-2)+G_padded(4:2:end,2:2:end-2)+G_padded(3:2:end-1,1:2:end-3));
% % 
G_bm = zeros(col,row,3);
G_bm(:,:,1) = G_r;
G_bm(:,:,2) = G_g;
G_bm(:,:,3) = G_b;
figure(2);
imshow(G_bm/255);
writeraw(G_bm,"demosaicing_a.raw",3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Problem1(b)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_padded2 = zeros(col+4,row+4);
G_padded2(3:end-2,3:end-2) = G;
G_padded2(1:2,3:end-2) = G_padded2(5:-1:4,3:end-2);
G_padded2(end-1:end,3:end-2) = G_padded2(end-3:-1:end-4,3:end-2);
G_padded2(:,1:2) = G_padded2(:,5:-1:4);
G_padded2(:,end-1:end) = G_padded2(:,end-3:-1:end-4);
% 
arfa = 0.5;
beta = 5/8;
gama = 3/4;
% 
G_r2 = zeros(col,row);
G_g2 = zeros(col,row);
G_b2 = zeros(col,row);
% 
for i =1:2:col
    for j = 1:2:row
        G_g2(i,j) = G_padded2(i+2,j+2);
        G_g2(i+1,j+1) = G_padded2(i+3,j+3);
        G_g2(i,j+1) = 1/8*(2*G_padded2(i+2,j+2)+2*G_padded2(i+2,j+4)+...
            2*G_padded2(i+1,j+3)+2*G_padded2(i+3,j+3)+4*G_padded2(i+2,j+3)-...
            G_padded2(i+2,j+1)-G_padded2(i+2,j+5)-G_padded2(i,j+3)-G_padded2(i+4,j+3));
        G_g2(i+1,j) = 1/8*(2*G_padded2(i+2,j+2)+2*G_padded2(i+3,j+1)+...
            2*G_padded2(i+4,j+2)+2*G_padded2(i+3,j+3)+4*G_padded2(i+3,j+2)-...
            G_padded2(i+1,j+2)-G_padded2(i+3,j)-G_padded2(i+2,j+4)-G_padded2(i+5,j+2));
        
        
        
        G_r2(i,j+1) = G_padded2(i+2,j+3);
        G_r2(i,j) = 1/8*(4*G_padded2(i+2,j+1)+4*G_padded2(i+2,j+3)+...
            5*G_padded2(i+2,j+2)+0.5*G_padded2(i,j+2)+0.5*G_padded2(i+4,j+2)-...
            G_padded2(i+2,j)-G_padded2(i+1,j+1)-G_padded2(i+3,j+1)-G_padded2(i+1,j+3)-...
            G_padded2(i+3,j+3)-G_padded2(i+2,j+4));
        G_r2(i+1,j+1) = 1/8*(4*G_padded2(i+2,j+3)+4*G_padded2(i+4,j+3)+...
            0.5*G_padded2(i+3,j+1)+0.5*G_padded2(i+3,j+5)+5*G_padded2(i+3,j+3)-...
            G_padded2(i+2,j+2)-G_padded2(i+2,j+4)-G_padded2(i+4,j+2)-G_padded2(i+4,j+4)-...
            G_padded2(i+1,j+3)-G_padded2(i+5,j+3));
        G_r2(i+1,j) = 1/8*(6*G_padded2(i+3,j+2)+2*G_padded2(i+2,j+1)+2*G_padded2(i+2,j+3)+...
            2*G_padded2(i+4,j+1)+2*G_padded2(i+4,j+3)-1.5*G_padded2(i+1,j+2)-...
            1.5*G_padded2(i+3,j)-1.5*G_padded2(i+3,j+4)-1.5*G_padded2(i+5,j+2));
        % G_r2(1:2:end,2:2:end) = G_padded2(3:2:end-3,4:2:end-2);
% G_r2(1:2:end,1:2:end) = 1/8*(4*G_padded2(3:2:end-3,2:2:end-4)+4*G_padded2(3:2:end-3,4:2:end-2)+...
%     5*G_padded2(3:2:end-3,3:2:end-3)+0.5*G_padded2(1:2:end-5,3:2:end-3)+0.5*G_padded2(5:2:end-1,3:2:end-3)-...
%     G_padded2(3:2:end-3,1:2:end-5)-G_padded2(2:2:end-4,2:2:end-4)-G_padded2(4:2:end-2,2:2:end-4)-G_padded2(2:2:end-4,4:2:end-2)-...
%     G_padded2(4:2:end-2,4:2:end-2)-G_padded2(3:2:end-3,5:2:end-1));
% G_r2(2:2:end,2:2:end) = 1/8*(4*G_padded2(3:2:end-3,4:2:end-2)+4*G_padded2(5:2:end-1,4:2:end-2)+...
%     0.5*G_padded2(4:2:end-2,2:2:end-4)+0.5*G_padded2(4:2:end-2,6:2:end)+5*G_padded2(4:2:end-2,4:2:end-2)-...
%     G_padded2(3:2:end-3,3:2:end-3)-G_padded2(3:2:end-3,5:2:end-1)-G_padded2(5:2:end-1,3:2:end-3)-G_padded2(5:2:end-1,5:2:end-1)-...
%     G_padded2(2:2:end-4,4:2:end-2)-G_padded2(6:2:end,4:2:end-2));
% G_r2(2:2:end,1:2:end) = 1/8*(6*G_padded2(4:2:end-2,3:2:end-3)+2*G_padded2(3:2:end-3,2:2:end-4)+2*G_padded2(3:2:end-3,4:2:end-2)+...
%     2*G_padded2(5:2:end-1,2:2:end-4)+2*G_padded2(5:2:end-1,4:2:end-2)-1.5*G_padded2(2:2:end-4,3:2:end-3)-...
%     1.5*G_padded2(4:2:end-2,1:2:end-5)-1.5*G_padded2(4:2:end-2,5:2:end-1)-1.5*G_padded2(6:2:end,3:2:end-3));
        
        
        G_b2(i+1,j) = G_padded2(i+3,j+2);
        G_b2(i,j) = 1/8*(4*G_padded2(i+1,j+2)+4*G_padded2(i+3,j+2)+5*G_padded2(i+2,j+2)+...
            0.5*G_padded2(i+2,j+4)+0.5*G_padded2(i+2,j)-...
            G_padded2(i,j+2)-G_padded2(i+4,j+2)-...
            G_padded2(i+1,j+1)-G_padded2(i+1,j+3)-G_padded2(i+3,j+3)-G_padded2(i+3,j+1));
        G_b2(i+1,j+1) = 1/8*(4*G_padded2(i+3,j+2)+4*G_padded2(i+3,j+4)+...
            5*G_padded2(i+3,j+3)+0.5*G_padded2(i+1,j+3)+0.5*G_padded2(i+5,j+3)-...
            G_padded2(i+2,j+2)-G_padded2(i+2,j+4)-G_padded2(i+4,j+4)-G_padded2(i+4,j+2)-...
            -G_padded2(i+3,j+1)-G_padded2(i+3,j+5));
        G_b2(i,j+1) = 1/8*(2*G_padded2(i+1,j+2)+2*G_padded2(i+1,j+4)+2*G_padded2(i+3,j+2)+2*G_padded2(i+3,j+4)+...
            6*G_padded2(i+2,j+3)-1.5*G_padded2(i,j+3)-1.5*G_padded2(i+2,j+5)-...
            1.5*G_padded2(i+2,j+1)-1.5*G_padded2(i+4,j+3));
        % G_b2(2:2:end,1:2:end) = G_padded2(4:2:end-2,3:2:end-3);
        % G_b2(1:2:end,1:2:end) = 1/8*(4*G_padded2(2:2:end-4,3:2:end-3)+4*G_padded2(4:2:end-2,3:2:end-3)+5*G_padded2(3:2:end-3,3:2:end-3)+...
        %     0.5*G_padded2(3:2:end-3,5:2:end-1)+0.5*G_padded2(3:2:end-3,1:2:end-5)-...
        %     G_padded2(1:2:end-5,3:2:end-3)-G_padded2(5:2:end-1,3:2:end-3)-...
        %     G_padded2(2:2:end-4,2:2:end-4)-G_padded2(2:2:end-4,4:2:end-2)-G_padded2(4:2:end-2,4:2:end-2)-G_padded2(4:2:end-2,2:2:end-4));
        % G_b2(2:2:end,2:2:end) = 1/8*(4*G_padded2(4:2:end-2,3:2:end-3)+4*G_padded2(4:2:end-2,5:2:end-1)+...
        %     5*G_padded2(4:2:end-2,4:2:end-2)+0.5*G_padded2(2:2:end-4,4:2:end-2)+0.5*G_padded2(6:2:end,4:2:end-2)-...
        %     G_padded2(3:2:end-3,3:2:end-3)-G_padded2(3:2:end-3,5:2:end-1)-G_padded2(5:2:end-1,5:2:end-1)-G_padded2(5:2:end-1,3:2:end-3)-...
        %     -G_padded2(4:2:end-2,2:2:end-4)-G_padded2(4:2:end-2,6:2:end));
        % G_b2(1:2:end,2:2:end) = 1/8*(2*G_padded2(2:2:end-4,3:2:end-3)+2*G_padded2(2:2:end-4,5:2:end-1)+2*G_padded2(4:2:end-2,3:2:end-3)+2*G_padded2(4:2:end-2,5:2:end-1)+...
        %     6*G_padded2(3:2:end-3,4:2:end-2)-1.5*G_padded2(1:2:end-5,4:2:end-2)-1.5*G_padded2(3:2:end-3,6:2:end)-...
        %     1.5*G_padded2(3:2:end-3,2:2:end-4)-1.5*G_padded2(5:2:end-1,4:2:end-2));
    end
end

% G_g2(1:2:end,1:2:end) = G_padded2(3:2:end-3,3:2:end-3);
% G_g2(2:2:end,2:2:end) = G_padded2(4:2:end-2,4:2:end-2);
% G_g2(1:2:end,2:2:end) = 1/8*(2*G_padded2(3:2:end-3,3:2:end-3)+2*G_padded2(3:2:end-3,5:2:end-1)+...
%     2*G_padded2(2:2:end-4,4:2:end-2)+2*G_padded2(4:2:end-2,4:2:end-2)+4*G_padded2(3:2:end-3,4:2:end-2)-...
%     G_padded2(3:2:end-3,2:2:end-4)-G_padded2(3:2:end-3,6:2:end)-G_padded2(1:2:end-5,4:2:end-2)-G_padded2(5:2:end-1,4:2:end-2));
% G_g2(2:2:end,1:2:end) = 1/8*(2*G_padded2(3:2:end-3,3:2:end-3)+2*G_padded2(4:2:end-2,2:2:end-4)+...
%     2*G_padded2(5:2:end-1,3:2:end-3)+2*G_padded2(4:2:end-2,4:2:end-2)+4*G_padded2(4:2:end-2,3:2:end-3)-...
%     G_padded2(2:2:end-4,3:2:end-3)-G_padded2(4:2:end-2,1:2:end-5)-G_padded2(4:2:end-2,5:2:end-1)-G_padded2(6:2:end,3:2:end-3));
% 
% 

% 
% 
% 

% 
G_mhc = zeros(col,row,3);
G_mhc(:,:,1) = G_r2;
G_mhc(:,:,2) = G_g2;
G_mhc(:,:,3) = G_b2;

figure(3);
imshow(G_mhc/255);
writeraw(G_mhc,"demosaicing_b.raw",3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Problem1(c)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = readraw('Toy.raw',400,560,3);
T_r = T(:,:,1);
T_g = T(:,:,2);
T_b = T(:,:,3);
T_rhis = zeros(1,256);
T_ghis = zeros(1,256);
T_bhis = zeros(1,256);
T_bucket_r=zeros(256,875,875);
T_bucket_g=zeros(256,875,875);
T_bucket_b=zeros(256,875,875);
for i = 1:560
    for j = 1:400
        T_rhis(1,T_r(i,j)+1) = T_rhis(1,T_r(i,j)+1)+1;
        T_ghis(1,T_g(i,j)+1) = T_ghis(1,T_g(i,j)+1)+1;
        T_bhis(1,T_b(i,j)+1) = T_bhis(1,T_b(i,j)+1)+1;
    end
end

for i =2:256
    T_rhis(1,i) = T_rhis(1,i)+T_rhis(1,i-1);
    T_ghis(1,i) = T_ghis(1,i)+T_ghis(1,i-1);
    T_bhis(1,i) = T_bhis(1,i)+T_bhis(1,i-1);
end

T_rhis  = 255.*T_rhis./(560*400);
T_ghis  = 255.*T_ghis./(560*400);
T_bhis  = 255.*T_bhis./(560*400);
figure(4);
plot(T_rhis);
figure(5);
plot(T_ghis);
figure(6);
plot(T_bhis);

for i = 1:560
    for j = 1:400
        T_r(i,j) = round(T_rhis(1,T_r(i,j)+1));
        T_g(i,j) = round(T_ghis(1,T_g(i,j)+1));
        T_b(i,j) = round(T_bhis(1,T_b(i,j)+1));
    end
end
T_rhis = zeros(1,256);
T_ghis = zeros(1,256);
T_bhis = zeros(1,256);
for i = 1:560
    for j = 1:400
        T_rhis(1,T_r(i,j)+1) = T_rhis(1,T_r(i,j)+1)+1;
        T_ghis(1,T_g(i,j)+1) = T_ghis(1,T_g(i,j)+1)+1;
        T_bhis(1,T_b(i,j)+1) = T_bhis(1,T_b(i,j)+1)+1;
    end
end
figure(4);
bar(T_rhis);
figure(5);
bar(T_ghis);
figure(6);
bar(T_bhis);

T_hm = zeros(560,400,3);
T_hm(:,:,1) = T_r;
T_hm(:,:,2) = T_g;
T_hm(:,:,3) = T_b;
figure(7);
imshow(T/255);
figure(8);
imshow(T_hm/255);
writeraw(T_hm,"histogram_c.raw",3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Problem2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Problem2(a)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = readraw("Corn_gray.raw",320,320,1);
C_n = readraw("Corn_noisy.raw",320,320,1);
figure(9);
imshow(C/255);
figure(10);
imshow(C_n/255);
Noise = zeros(1,510);
for i =1:320
    for j =1:320
        Noise(1,C_n(i,j)-C(i,j)+255) = Noise(1,C_n(i,j)-C(i,j)+255)+1;
    end
end
figure(11);
plot(Noise(191:271));
h1 = [0.1,0.1,0.1;0.1,0.2,0.1;0.1,0.1,0.1];
h_temp1 = [-0.5,-0.5,-0.5;0,0,0;-0.5,-0.5,-0.5];
h_temp2 = [-0.5,0,-0.5;-0.5,0,-0.5;-0.5,0,-0.5];

h2 = exp(h_temp1+h_temp2);
A = sum(sum(h2));
h2 = h2/A;
K = convolution(C_n,h1);
K_g = convolution(C_n,h2);
figure(12);
imshow(K/255);
writeraw(K,"basic_uniform.raw",1);
figure(13);
imshow(K_g/255);
writeraw(K_g,"basic_gaussain.raw",1);
psnr = PSNR(C,K);
psnr_g = PSNR(C,K_g);
PSNR_o = PSNR(C,C);
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
J = Bi_filter(C_n,100,100,1);
figure(14);
imshow(J/255);
psnr_b  =PSNR(C,J);
writeraw(J,"bilateral.raw",1);
% 
% 
J2 = nlm_filter(C_n,100,100,1,1,1);
figure(15);
imshow(J2/255);
psnr_n = PSNR(C,J2);
writeraw(J2,"nlm_filtering.raw",1);



















