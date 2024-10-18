%%SPIKING NEURONS FOLLOWING THE Izhikevich MODEL
close all
clear all
clf
%These are some default parameter values
I=10;
spiketype=5;
spikeName={'Regular spiking','Intrinsically bursting behavior','Chattering','Fast spiking','Low-threshold spiking'};

RS=[0.02 0.2 -65 8]; %regular spiking
IB = [0.02, 0.2, -55, 4]; %intrinsically bursting behavior
CH = [0.02, 0.2, -50, 2]; %CH
FS = [0.1, 0.2, -65, 2]; %fast spiking
LTS = [0.1, 0.25, -65, 2];%Low-threshold spiking
abcd=[RS;IB;CH;FS;LTS]; %Rows: type of dynamics; Columns: [a,b,c,d]

a=abcd(spiketype,1);
b=abcd(spiketype,2);
c=abcd(spiketype,3);
d=abcd(spiketype,4);

%The initial values for v and u
v=-65;
u=b*v;

%Initialize the vector that will contain the membrane potential time series.
v_tot=zeros(1000, 1);

for t=1:1000
    %set v_tot at this time point to the current value of v
    v_tot(t)=v;
    %Reset v and u if v has crossed threshold. See Eq. 3 above.
    if (v>= 30)
    v=c;
    u=u+d;
    end;
    %Use Euler’s method to integrate Eqs. 1 and 2 from above. Here v is calculated in 
    %2 steps in order to keep the time step small (0.5 ms step in the line below).
    v=v+0.5*(0.04*v^2+5*v+140-u+I);
    v=v+0.5*(0.04*v^2+5*v+140-u+I);
    u=u+a*(b*v-u);
end;
%This line uses the function find to locate the indices of v_tot that hold elements with 
%values greater than or equal to 30 and then sets these elements to 30.
%This normalizes to heights of the action potential peaks to 30.
v_tot(find(v_tot >= 30))=30;
%Plot the neuron’s membrane potential.
time=1:1000;
figure(1)
plot(time(1:200), v_tot(1:200));
xlabel('Time (ms)', 'fontsize', 12);
ylabel('Membrane Potential (mV)', 'fontsize', 12);
title(spikeName(spiketype), 'fontsize', 14);
