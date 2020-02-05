%EE569 HOMEWORK ASSGINMENT #1
%Date: 2020/01/28
%Name: Shen Zhihong
%ID: 3645974217
%email: shenzhih@usc.edu


function psnr = PSNR(O_p,G_p)
%PSNR 此处显示有关此函数的摘要
%   此处显示详细说明
temp = 0;
for i  = 1: size(O_p,1)
    for j = 1:size(O_p,2)
        temp = temp+(O_p(i,j)-G_p(i,j))^2;
    end
end
temp = temp/(size(O_p,1)*size(O_p,2));
psnr = 10*log10(255^2/temp);
end

