
%EE569 HOMEWORK ASSGINMENT #1
%Date: 2020/01/28
%Name: Shen Zhihong
%ID: 3645974217
%email: shenzhih@usc.edu


function J = nlm_filter(I,h,a,filter_R,n1,n2)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
J = zeros(size(I,1),size(I,2));
K = zeros(size(I,1)+2*(filter_R+n1),size(I,2)+2*(filter_R+n2));
K(n1+filter_R+1:n1+filter_R+size(I,1),n2+filter_R+1:n2+filter_R+size(I,2)) = I;
weight = zeros(size(I,1),size(I,2),2*filter_R+1,2*filter_R+1);
%f(i,j,k,l)  i = i+2*filter_R   k = [i-filter_R,i+filter_R]
for i =1:size(I,1)
    for j =1:size(I,1)
        temp1 = 0;
        temp2 = 0;
        for k = -filter_R:filter_R
            for l = -filter_R:filter_R
                temp3 = 0;
                for n_1 = -n1:n1
                    for n_2 = -n2:n2
                        g_a = 1/(sqrt(2*pi)*a)*exp(-(n1^2+n2^2)/(2*a^2));
                        dis = (K(i+filter_R+n1-n_1,j+filter_R+n2-n_2)-K(i+filter_R+k+n1-n_1,j+filter_R+l+n2-n_2))^2;
                        temp3 = temp3+dis*g_a;
                    end
                end
                weight(i,j,filter_R+k+1,j+l+filter_R) = exp(-temp3^2/h^2);
                temp2 = temp2+weight(i,j,filter_R+k+1,j+l+filter_R);
                temp1 = temp1 +K(i+filter_R+n1+k,j+filter_R+n2+l)*weight(i,j,filter_R+k+1,j+l+filter_R);
            end
            J(i,j) = temp1/temp2;
        end
    end
end
end

