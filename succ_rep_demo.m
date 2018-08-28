clear all

%Initialize parameters
buff_len = 51;  k = 10;  N = buff_len+2*k; Taustar_min = 1; Taustar_max = 10; 

%Create power-law growing Taustarlist and corresponding s
alpha = (Taustar_max/Taustar_min)^(1/buff_len)-1;
Taustarlist = Taustar_min * (1+alpha).^(-(k):(buff_len +(k) -1));
s = k./Taustarlist;


%Create DerivMatrix
DerivMatrix = zeros(N,N);
 for i = 2:N-1
   DerivMatrix(i,i-1) = -(s(i+1)-s(i))/(s(i)-s(i-1))/(s(i+1) - s(i-1));
   DerivMatrix(i,i) = ((s(i+1)-s(i))/(s(i)- s(i-1))-(s(i)-s(i-1))/(s(i+1)-s(i)))/(s(i+1) - s(i-1));
   DerivMatrix(i,i+1) = (s(i)-s(i-1))/(s(i+1)-s(i))/(s(i+1) - s(i-1));
 end

gamma = exp(-s);
 
%stim_seq = ['B', 'L', 'G', 'P', 'G', 'R'];
stim_seq = ['O', 'N', 'B', 'P', 'C', 'C', 'C', 'B', 'D', 'G'];
stimuli = unique(stim_seq);

target_stims1 = ['O', 'B'];
SR1 = compute_succ_rep(stim_seq, target_stims1, s);
SR_diff = DerivMatrix^k*SR1';
T1 = (-1)^k*s.^(k+1)'.*SR_diff/factorial(k)';

target_stims2 = ['O', 'C'];
SR2 = compute_succ_rep(stim_seq, target_stims2, s);
SR_diff = DerivMatrix^k*SR2';
T2 = (-1)^k*s.^(k+1)'.*SR_diff/factorial(k)';

%Plot the example

figure, subplot(3,1,1), plot([0:9],ones(10,1),'ok', 'MarkerSize', 12)
text([0:9],ones(10,1),stim_seq', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 10 );
set(gca,'xlim',[-1 11])
set(gca,'XTick',[0:9])
set(gca,'YTick',[])
set(gca,'YTickLabel',{[]})
xlabel('Time')
title('Sequence of simuli')
hold on, plot(find(stim_seq == target_stims1(2))-1,ones(length(find(stim_seq == target_stims1(2)))),...
    'or', 'MarkerSize', 12)
plot(find(stim_seq == target_stims2(2))-1,ones(length(find(stim_seq == target_stims2(2)))),...
    'ob', 'MarkerSize', 12)

subplot(3,1,2), plot(gamma(k+1:end-k),SR1(k+1:end-k),'r.')
hold on, plot(gamma(k+1:end-k),SR2(k+1:end-k),'b.')
title('Nodes of the successor representation for different values of \gamma')
xlabel('\gamma')
ylabel('Activation')
set(gca,'xlim',[0 0.45])
legend('M^{O,B}','M^{O,C}','location','northeast')

subplot(3,1,3), plot(1000+Taustarlist(k+1:end-k)*1000,T1(k+1:end-k,end),'r.')
hold on, plot(1000+Taustarlist(k+1:end-k)*1000,T2(k+1:end-k,end),'b.')
title('Nodes of the inverse transform, each firing for characteristic internal time')
set(gca,'xlim',[0 12000])
set(gca,'XTick',[1000:1000:10000])
set(gca,'XTickLabel',{[0:9]})
xlabel('Internal future-time representation')
ylabel('Activation')
legend('M^{O,B}','M^{O,C}','location','northeast')

set(gcf,'color','w')

%export_fig('succ_rep_sequence_demo_v1.pdf')
%saveas(gcf, 'succ_rep_sequence_demo_v1.fig')


function SR = compute_succ_rep(stim_seq, target_stims, s)
    for i = 1:2
        ind{i} = find(stim_seq == target_stims(i));
    end
    for l = 1:length(s)
        SR(l) = 0;
        for i=1:length(ind{2})
            SR(l) = SR(l) + exp(-s(l)*abs(ind{1,1}-ind{2}(i)));
        end
    end
end

 
        