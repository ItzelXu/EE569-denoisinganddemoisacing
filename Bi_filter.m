%EE569 HOMEWORK ASSGINMENT #1
%Date: 2020/01/28
%Name: Shen Zhihong
%ID: 3645974217
%email: shenzhih@usc.edu


function J = Bi_filter(I,sigma_c,sigma_s,filter_R)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
J = zeros(size(I,1),size(I,2));
K = zeros(size(I,1)+2*filter_R,size(I,2)+2*filter_R);
K(filter_R+1:filter_R+size(I,1),filter_R+1:filter_R+size(I,2)) = I;
weight = zeros(size(I,1),size(I,2),2*filter_R+1,2*filter_R+1);
for i = 1:size(I,1)
    for j = 1:size(I,2)
        temp1 = 0;
        temp2 = 0;
        for k = -filter_R:filter_R
            for l = -filter_R:filter_R
                weight(i,j,k+filter_R+1,l+filter_R+1) = exp(-((k^2+l^2)/(2*sigma_c^2))-((K(i+filter_R,j+filter_R)-K(i+k+filter_R,j+l+filter_R))^2/(2*sigma_s^2)));
                temp1 = temp1+K(i+k+filter_R,l+j+filter_R)*weight(i,j,k+filter_R+1,l+filter_R+1);
                temp2 = temp2+weight(i,j,k+filter_R+1,l+filter_R+1);
            end
        end
        J(i,j) = temp1/temp2;
    end
end
end

