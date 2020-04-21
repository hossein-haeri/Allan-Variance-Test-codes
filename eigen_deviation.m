clear all, clc, close all

L = 1;
b0 = 1.5;
g = 10;
m0 = 2;
I0 = m0*L^2;
n = 100;
num_window_samples = 5;
t = linspace(0,6*pi,n);

m = normrnd(m0, 0.0, 1,n);
b = normrnd(b0, 0.5, 1,n)+2*sin(2*t);


% v = VideoWriter('video_mvg.avi');
% v.FrameRate = 10;
% open(v);
figure('position',[0 0 800 400])
hold on
A = [0 1; m0.*g./I0 -b0./I0];
B = [0; 1/I0];
win = [];
LAM = [];
E = [];
V = [];
w = 50;
for num_window_samples= 1:w
for i=1:n
    I = m(i)*L^2;
    A_tilde = [0 1; m(i).*g./I -b(i)./I];
    B_tilde = [0; 1/I];
    
    if i > num_window_samples
        b_hat = mean(b(i-num_window_samples:i));
        m_hat = mean(m(i-num_window_samples:i));
    else
        b_hat = mean(b(1:i));
        m_hat = mean(m(1:i));
    end
    
    
    I = m_hat*L^2;
    A_hat = [0 1; m_hat.*g./I -b_hat./I];
    B_hat = [0; 1/I];
    
    K = place(A_hat, B_hat, [-2, -4]);
    
    A_cl = A - B*K;
    A_cl_hat = A_hat - B_hat*K;
    
    lam = eig(A_cl);
    lam_hat = eig(A_cl_hat);
    
    e = norm(lam - lam_hat);
    E = [E e];
    
%     cla()
%     
% %     figure(1)
% %     LAM = [LAM real(lam)];
%     scatter([-2, -4], [0, 0],60,  '*', 'b');
%     scatter(real(lam_hat), imag(lam_hat),100, 'o','r','LineWidth',2);
%     scatter(real(lam), imag(lam),150, 'o', 'k','MarkerEdgeAlpha',0.8,'LineWidth',2);
%     xL = xlim;
%     yL = ylim;
%     line([0 0], yL, 'Color', [0.8,0.8,0.8],'LineWidth',2);  %x-axis
%     line(xL, [0 0], 'Color', [0.8,0.8,0.8],'LineWidth',2);  %y-axis
% 
%     xlim([-5,1])
%     xlabel('Re')
%     ylabel('Im')
% 
% 
%     grid on
%     legend('$\lambda^*$','$\hat{\lambda}$','$\lambda_{true}$','FontSize',16)  
%     pause(0.10)
%     frame = getframe(gcf);
%     writeVideo(v,frame);
    
end
v = var(E);
V = [V v];
end
% close(v);

plot(V, 'LineWidth',2)
xlabel('Window length M')
ylabel('Eigenvalue variance')
% figure(2)
% histogram(LAM(1,:),50)

