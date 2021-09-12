%OPTION PARAMETERS%
S = 50;
K = 50;
r = 0.01;
T=1;
sigma=0.25;
NSteps=24;
%Simulating with Halton random numbers%
N=1000;
G = zeros(1,N);
for i=1:N
    q = qrandstream('halton',NSteps,'Skip',1e3,'Leap',1e2);
    RandMat = qrand(q,i);
    z_RandMat = norminv(RandMat,0,1);
    dt = T/NSteps;
    C = (r-0.5*sigma^2)*dt + sigma*sqrt(dt)*z_RandMat;
    Paths = cumsum([log(S0)*ones(i,1),C],2);
    % Underlying asset price path
    SPaths = exp(Paths);
    Payoff = zeros(i,1);
    Payoff = max(0, mean(SPaths(:,2:(NSteps+1)),2)- K);
    G(i) = mean(exp(-r*T) * Payoff);
end 
figure(1)
plot(G,'-b')
grid minor
grid on
title('Convergence Diagram (Asian Option)')
ylabel('Option Price (Numerical)')
xlabel('No. of Simulations')
lgd=sprintf('$ %0.5f',G(end));
legend({lgd})
saveas(gcf,'Asian_MC_Halton','png')