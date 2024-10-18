%%NETWORK MODEL OF SPIKING NEURONS FOLLOWING THE Izhikevich MODEL
close all
clear
clf

% The number of excitatory neurons in the network. The mammalian cortex has about 4 times
%as many excitatory neurons as inhibitory ones.
Ne = 800;
%The number of inhibitory neurons in the network.
Ni = 200;

%Subpopulations of excitatory neurons:
p = 5;
f = 0.1;
n_esp = f*Ne; %number of neurons in each specific subpopulation
n_inesp = (1-f*p)*Ne; %number of neurons in the unspecific subpopulation

%Random numbers
re=rand(Ne, 1);
ri=rand(Ni, 1);

%This will set the value of 'a' for all excitatory neurons to 0.02 and the value of 'a'
%for inhibitory neurons to a random number between 0.02 and 0.1
a = [0.02*ones(Ne, 1); 0.02+0.08*ri];
%This will allow b to range from 0.2–0.25
b = [0.2*ones(Ne, 1); 0.25-0.05*ri];
%This will allow the spike reset membrane potential to range between -65 and -50
c = [-65+15*re.^2; -65*ones(Ni,1)];
%This will allow the recovery reset value to range between 2 and 8
d = [8-6*re.^2; 2*ones(Ni, 1)];

%S: weight matrix, holds info about the strength of connections: (it is a fully connected net)
%HEBBIAN RULE
win = 0.5;%2.1; %between excitatory of the same specific group.
wout = 1-(f*(win-1)/(1-f)); %between excitatory from different groups.
wns = 0.25;% between non selective;
S = ones(Ne+Ni,Ne+Ni); %square matrix of Ne+Ni *Ne+Ni

spi=0;
spj=0;
for i=1:(Ne)
    for j=1:(Ne)
        for k=1:5
            if i<=k*80
                spi=k;
                break
            else
                spi=6;
            end
            
        end
        for k=1:5
            if j<=k*80
                spj=k;
                break
            else
                spj=6;
            end
        end
        if spi==spj && spi<6
            S(i,j)=win*rand(1);
        elseif spi==spj && spi==6 
            S(i,j)=wns*rand(1);
        else
            S(i,j)=wout*rand(1);
        end
    end
end
for i=1:Ne
    for j=Ne+1:(Ne+Ni)
        %S(i,j)=wout*rand(1);
        S(i,j)=-rand(1);
    end
end
for j=1:Ne+Ni
    for i=Ne+1:(Ne+Ni)
        S(i,j)=-rand(1);
    end
end
    
% Very few elements of S will be exactly 0, so in this model almost every neuron has 
% synaptic contacts with all other neurons in the network.

%The initial values for v and u
v = -65*ones(Ne+Ni,1);
u = b.*v;

%Firings will be a two-column matrix. The first column will indicate the time that a
%neuron’s membrane potential crossed 30, and the second column will be a number
%between 1 and Ne+Ni that identifies which neuron fired at that time.
firings=[];

%stimulus:
stimulus1=1; %0-5; 0 means no stimulus.
N_range=[1,80; 81,160; 161,240; 241,320; 321,400; 401,800; 801,1000]; %index ranges for each populaton type.
time_s1= 2000:2500;%time when the stimulus is first applied.
stimulus2=0; %0-5; 0 means no stimulus.
time_s2= 6000:6500;%time when the stimulus is reapplied.
end_experiment=6500;
for t=1:end_experiment
    %Create some random input external to the network
    I=[5*randn(Ne, 1); 2*randn(Ni,1)]; %thalamic input
    %First Stimulus to a specific supbopulation
    if stimulus1 ==1
        if ismember(t,time_s1)==1
            I_ext=zeros(Ne+Ni,1);
            I_ext(N_range(stimulus1,1):N_range(stimulus1,2))=500;
            I=I+I_ext;%total input
        end       
    elseif stimulus1 ==2
        if ismember(t,time_s1)==1
            I_ext=zeros(Ne+Ni,1);
            I_ext(N_range(stimulus1,1):N_range(stimulus1,2))=500;
            I=I+I_ext;%total input
        end 
    elseif stimulus1 ==3
        if ismember(t,time_s1)==1
            I_ext=zeros(Ne+Ni,1);
            I_ext(N_range(stimulus1,1):N_range(stimulus1,2))=500;
            I=I+I_ext;%total input
        end 
    elseif stimulus1 ==4
        if ismember(t,time_s1)==1
            I_ext=zeros(Ne+Ni,1);
            I_ext(N_range(stimulus1,1):N_range(stimulus1,2))=500;
            I=I+I_ext;%total input
        end 
    elseif stimulus1 ==5
        if ismember(t,time_s1)==1
            I_ext=zeros(Ne+Ni,1);
            I_ext(N_range(stimulus1,1):N_range(stimulus1,2))=500;
            I=I+I_ext;%total input
        end 
    end
    %Second stimulus to a specific supbopulation
    if stimulus2~=0
        if ismember(t,time_s2)==1
            I_ext=zeros(Ne+Ni,1);
            I_ext(1:400)=500;
            I=I+I_ext;%total input
        end 
    end
       
            
    %Determine which neurons crossed threshold at the current time step t.
    fired=find(v>=30); % indices of spikeS
    %!!!fired=find((v+65)>=30); % indices of spikeS
    %Add the times of firing and the neuron number to firings.
    times= t*ones(1, length(fired));
    tn= horzcat(times',fired);
    if isempty(tn)
        firings = firings;
    else
        firings=cat(1,firings,tn);
    end
    %Reset the neurons that fired to the spike reset membrane potential and
    %recovery variable.
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    %Add to the input, I, for each neuron a value equal to the sum of the synaptic
    %strengths of all other neurons that fired in the last time step connected to that
    %neuron.
    I=I+sum(S(:,fired), 2);

    %Move the simulation forward using Euler’s method. Step=0.5 for numerical stability
    v=v+0.5*(0.04*v.^2+5*v+140-u+I);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I);
    u=u+a.*(b.*v-u);
    
end

%% Plot the raster plot of the network activity.
% We will randomly select 4#N neurons of each type to represent
U=unique(firings(:,2));
%First specific excitatory group: U(1-80)
E1_1 = firings(find(firings(:,2)==U(2)),:);
E1_2 = firings(find(firings(:,2)==U(33)),:);
E1_3 = firings(find(firings(:,2)==U(56)),:);
E1_4 = firings(find(firings(:,2)==U(78)),:);
%Second specific excitatory group: U(81-160)
E2_1 = firings(find(firings(:,2)==U(2+80)),:);
E2_2 = firings(find(firings(:,2)==U(33+80)),:);
E2_3 = firings(find(firings(:,2)==U(56+80)),:);
E2_4 = firings(find(firings(:,2)==U(78+80)),:);
%Third specific excitatory group:
E3_1 = firings(find(firings(:,2)==U(2+80*2)),:);
E3_2 = firings(find(firings(:,2)==U(33+80*2)),:);
E3_3 = firings(find(firings(:,2)==U(56+80*2)),:);
E3_4 = firings(find(firings(:,2)==U(78+80*2)),:);
%Forth specific excitatory group:
E4_1 = firings(find(firings(:,2)==U(2+80*3)),:);
E4_2 = firings(find(firings(:,2)==U(33+80*3)),:);
E4_3 = firings(find(firings(:,2)==U(56+80*3)),:);
E4_4 = firings(find(firings(:,2)==U(78+80*3)),:);
%Fifth specific excitatory group:
E5_1 = firings(find(firings(:,2)==U(2+80*4)),:);
E5_2 = firings(find(firings(:,2)==U(33+80*4)),:);
E5_3 = firings(find(firings(:,2)==U(56+80*4)),:);
E5_4 = firings(find(firings(:,2)==U(78+80*4)),:);
%Inspecific excitatory group:
Eis_1 = firings(find(firings(:,2)==U(402)),:);
Eis_2 = firings(find(firings(:,2)==U(466)),:);
Eis_3 = firings(find(firings(:,2)==U(492)),:);
Eis_4 = firings(find(firings(:,2)==U(502)),:);
Eis_5 = firings(find(firings(:,2)==U(567)),:);
Eis_6 = firings(find(firings(:,2)==U(589)),:);
Eis_7 = firings(find(firings(:,2)==U(607)),:);
Eis_8 = firings(find(firings(:,2)==U(652)),:);
Eis_9 = firings(find(firings(:,2)==U(698)),:);
Eis_10 = firings(find(firings(:,2)==U(706)),:);
Eis_11 = firings(find(firings(:,2)==U(755)),:);
Eis_12 = firings(find(firings(:,2)==U(783)),:);
Eis_13 = firings(find(firings(:,2)==U(532)),:);
Eis_14 = firings(find(firings(:,2)==U(576)),:);
Eis_15 = firings(find(firings(:,2)==U(675)),:);
Eis_16 = firings(find(firings(:,2)==U(723)),:);

%Inhibitory  group:
I_1 = firings(find(firings(:,2)==U(802)),:);
I_2 = firings(find(firings(:,2)==U(806)),:);
I_3 = firings(find(firings(:,2)==U(875)),:);
I_4 = firings(find(firings(:,2)==U(823)),:);
I_5 = firings(find(firings(:,2)==U(967)),:);
I_6 = firings(find(firings(:,2)==U(999)),:);
I_7 = firings(find(firings(:,2)==U(907)),:);
I_8 = firings(find(firings(:,2)==U(952)),:);
I_9 = firings(find(firings(:,2)==U(898)),:);
I_10 = firings(find(firings(:,2)==U(906)),:);
I_11 = firings(find(firings(:,2)==U(955)),:);
I_12 = firings(find(firings(:,2)==U(855)),:);

figure(1)
plot(E1_1(:,1),E1_1(:,2),'b.',E1_2(:,1),E1_2(:,2),'b.',E1_3(:,1),E1_3(:,2),'b.',E1_4(:,1),E1_4(:,2),'b.');
hold on
plot(E2_1(:,1),E2_1(:,2),'b.',E2_2(:,1),E2_2(:,2),'b.',E2_3(:,1),E2_3(:,2),'b.',E2_4(:,1),E2_4(:,2),'b.');
plot(E3_1(:,1),E3_1(:,2),'b.',E3_2(:,1),E3_2(:,2),'b.',E3_3(:,1),E3_3(:,2),'b.',E3_4(:,1),E3_4(:,2),'b.');
plot(E4_1(:,1),E4_1(:,2),'b.',E4_2(:,1),E4_2(:,2),'b.',E4_3(:,1),E4_3(:,2),'b.',E4_4(:,1),E4_4(:,2),'b.');
plot(Eis_1(:,1),Eis_1(:,2),'b.',Eis_2(:,1),Eis_2(:,2),'b.',Eis_3(:,1),Eis_3(:,2),'b.',Eis_4(:,1),Eis_4(:,2),'b.',Eis_5(:,1),Eis_5(:,2),'b.',Eis_6(:,1),Eis_6(:,2),'b.',Eis_7(:,1),Eis_7(:,2),'b.',Eis_8(:,1),Eis_8(:,2),'b.',Eis_9(:,1),Eis_9(:,2),'b.',Eis_10(:,1),Eis_10(:,2),'b.',Eis_11(:,1),Eis_11(:,2),'b.',Eis_12(:,1),Eis_12(:,2),'b.',Eis_13(:,1),Eis_13(:,2),'b.',Eis_14(:,1),Eis_14(:,2),'b.',Eis_15(:,1),Eis_15(:,2),'b.',Eis_16(:,1),Eis_16(:,2),'b.');
plot(I_1(:,1),I_1(:,2),'b.',I_2(:,1),I_2(:,2),'b.',I_3(:,1),I_3(:,2),'b.',I_4(:,1),I_4(:,2),'b.',I_5(:,1),I_5(:,2),'b.',I_6(:,1),I_6(:,2),'b.',I_7(:,1),I_7(:,2),'b.',I_8(:,1),I_8(:,2),'b.',I_9(:,1),I_9(:,2),'b.',I_10(:,1),I_10(:,2),'b.',I_11(:,1),I_11(:,2),'b.',I_12(:,1),I_12(:,2),'b.');
yline(80,'r')
yline(160,'r')
yline(240,'r')
yline(320,'r')
yline(400,'r')
yline(800,'r')

xlabel('time (ms)') 
ylabel('Neuron index') 
hold off

figure(2)
plot(firings(:,1), firings(:,2),'.');
xlabel('time (ms)') 
ylabel('Neuron index') 

%One can appreciate cortical-like asynchronous dynamics.






%% Calculation of firing rates for each type of neuron:
spiketimes = firings(:,1);
spiketimes = unique(spiketimes);

is_active=zeros(Ne+Ni,end_experiment);

for j=1:end_experiment
    if ismember(j, spiketimes)
        id=find(firings(:,1)==j);
        %for id=find(firings(:,1)==j)
        is_active(firings(id,2),j)=1;
        %end
    end
end

step=50;
%initializing firing rates variables.
n1=zeros(1,length(1:step:(end_experiment-step))); %#active neurons/10msec in the first specific exciting subpopulation of neurons 
n2=zeros(1,length(1:step:(end_experiment-step)));%#active neurons/10msec in the second specific exciting subpopulation of neurons 
n3=zeros(1,length(1:step:(end_experiment-step))); %#active neurons/10msec in the third specific exciting subpopulation of neurons 
n4=zeros(1,length(1:step:(end_experiment-step))); %#active neurons/10msec in the forth specific exciting subpopulation of neurons 
n5=zeros(1,length(1:step:(end_experiment-step))); %#active neurons/10msec in the fifth specific exciting subpopulation of neurons 
n6=zeros(1,length(1:step:(end_experiment-step)));%#active neurons/10msec in the non specific exciting subpopulation of neurons 
n7=zeros(1,length(1:step:(end_experiment-step))); %#active neurons/10msec in the inhibitor subpopulation of neurons 

i=1;
 for k=1:step:(end_experiment-step)
    n1(i)=(1/step)*sum(sum(is_active(1:80,k:(k+step))));
    n2(i)=(1/step)*sum(sum(is_active(81:160,k:(k+step))));
    n3(i)=(1/step)*sum(sum(is_active(161:240,k:(k+step))));
    n4(i)=(1/step)*sum(sum(is_active(241:320,k:(k+step))));
    n5(i)=(1/step)*sum(sum(is_active(321:400,k:(k+step))));
    n6(i)=(1/step)*sum(sum(is_active(401:800,k:(k+step))));
    n7(i)=(1/step)*sum(sum(is_active(801:1000,k:(k+step))));
    i=i+1;
 end

n = [n1./n_esp; n2./n_esp; n3./n_esp; n4./n_esp; n5./n_esp; n6./n_inesp; n7./Ni];

figure(3)
%Plot the firing rates for the different neuron types in the network.
time=1:step:end_experiment-step;
plot(time,n(1,:), time,n(2,:), time,n(3,:), time,n(4,:), time,n(5,:), time,n(6,:), time,n(7,:))
xlabel('time (ms)') 
ylabel('Firing rate') 
legend('Excitatory subpopulation 1', 'Excitatory subpopulation 2', 'Excitatory subpopulation 3','Excitatory subpopulation 4', 'Excitatory subpopulation 5','Excitatory inespecific subpopulation', 'Inhibitory neurons');


 
 
 
 
 