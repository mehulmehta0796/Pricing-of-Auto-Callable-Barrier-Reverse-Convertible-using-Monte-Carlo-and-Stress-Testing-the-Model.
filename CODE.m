N=10^5; %No. of Simulation.
T1=0.50; %6th Month i.e. 13 October 2021
T2=0.75; %9th Month i.e. 13 January 2022
T3=1; %12th Month i.e. 11 April 2022
T4=1.25; %15th Month i.e. 13 July 2022

m=100; %No. of discretization.
dt=T4/m;
r=-0.00746; %Switzerland 1 Year Bond on 12th April 2021

S01=218.75; %Stock Price of Allianz SE on 12th April 2021(Stock 1)
S02=23.70; %Stock Price of AXA SA Bearer on 12th April 2021(Stock 2)
S03=157.10; %Stock Price of Hannover Rück SE on 12th April 2021(Stock 3)

K1=217.175; %Strike Price for Stock 1
K2=22.9075; %Strike Price for Stock 2
K3=156.575; %Strike Price for Stock 3

S1=S01*ones(m+1,N); 
S2=S02*ones(m+1,N);
S3=S03*ones(m+1,N);

B1=0.59*K1; %Barrier level(59%) for Stock 1
B2=0.59*K2; %Barrier level(59%) for Stock 2
B3=0.59*K3; %Barrier level(59%) for Stock 3

d1=0.0444; %Dividend Yield for Stock 1
d2=0.0617; %Dividend Yield for Stock 2
d3=0.0293; %Dividend Yield for Stock 3

CR1=4.6046; %Conversion Ratio for Stock 1
CR2=43.6538; %Conversion Ratio for Stock 2
CR3=6.3867; %Conversion Ratio for Stock 3

IV=1000; %Investment or Denomination
IA=0.00; %Interest Amount
PA=0.08; %Considered the Premium Amount as 8 percent

sig1=0.00780911*sqrt(252); %Standard Deviation for Stock 1
sig2=0.008202997*sqrt(252); %Standard Deviation for Stock 2
sig3=0.007596542*sqrt(252); %Standard Deviation for Stock 3

corr_mtrx=[1 0.874028799 0.798699928; 0.874028799 1 0.759168398; 0.798699928 0.759168398 1];
A=chol(corr_mtrx,'lower');

Z01=normrnd(0,1,m,N); %Creating a matrix of standard normal random variables. 
Z02=normrnd(0,1,m,N);
Z03=normrnd(0,1,m,N);

Z1=A(1,1)*Z01;
Z2=A(2,1)*Z01+A(2,2)*Z02;
Z3=A(3,1)*Z01+A(3,2)*Z02+A(3,3)*Z03; 

M1=zeros(m,N); %To check the extreme values of the stock 1. 
M2=zeros(m,N); %To check the extreme values of the stock 2.
M3=zeros(m,N); %To check the extreme values of the stock 3.

U1=unifrnd(0,1,m,N);
U2=unifrnd(0,1,m,N);
U3=unifrnd(0,1,m,N);

Barrier=zeros(m,N); %Initializing the Barrier Matrix
S_worst=zeros(1,N); %Capturing the Worst Performing Stock
Payoff=zeros(1,N); 

for n=1:N
    for i=1:m
        %Simulate the stock price using Geometric Brownian motion 
        S1(i+1,n)=S1(i,n).*exp((r-d1-0.5*sig1^2)*dt+sig1*sqrt(dt)*Z1(i,n));
        S2(i+1,n)=S2(i,n).*exp((r-d2-0.5*sig2^2)*dt+sig2*sqrt(dt)*Z2(i,n));
        S3(i+1,n)=S3(i,n).*exp((r-d3-0.5*sig3^2)*dt+sig3*sqrt(dt)*Z3(i,n));
        M1(i,n)=exp(0.5*(log(S1(i+1,n)*S1(i,n))-sqrt((log(S1(i+1,n)/S1(i,n)))^2-2*sig1^2*dt*log(U1(i,n)))));
        M2(i,n)=exp(0.5*(log(S2(i+1,n)*S2(i,n))-sqrt((log(S2(i+1,n)/S2(i,n)))^2-2*sig2^2*dt*log(U2(i,n)))));
        M3(i,n)=exp(0.5*(log(S3(i+1,n)*S3(i,n))-sqrt((log(S3(i+1,n)/S3(i,n)))^2-2*sig3^2*dt*log(U3(i,n)))));
        %Checking if the stock price was above the barrier.
        if M1(i,n)>B1 && M2(i,n)>B2 && M3(i,n)>B3 
            Barrier(i,n)=0; %Barrier is not achievied.
        else
            Barrier(i,n)=1; %Barrier is achieved.
        end    
    end
end

Check1=sum(Barrier); %Counting how many times the barrier has occured for every simulation.

CP1=exp(-r*0.25).*(0.25*PA*IV); %Premium Payment on 13 July 2021
CP2=exp(-r*T1).*(0.25*PA*IV); %Premium Payment on 13 October 2021
CP3=exp(-r*T2).*(0.25*PA*IV); %Premium Payment on 13 January 2022
CP4=exp(-r*T3).*(0.25*PA*IV); %Premium Payment on 11 April 2022
CP5=exp(-r*T4).*(0.25*PA*IV); %Premium Payment on 13 July 2022

S1_T1=S1(0.4*m+1,:); %Price of Stock S1 at T1 (6th month/15th month=0.4)
S1_T2=S1(0.6*m+1,:); %Price of Stock S1 at T2 (9/15=0.6)
S1_T3=S1(0.8*m+1,:); %Price of Stock S1 at T3 (12/15=0.8) 
S1_T4=S1(1*m+1,:); %Price of Stock S1 at T4 (15/15=1)

%Similarly calculated the price of S2 at T1, T2, T3 and T4.
S2_T1=S2(0.4*m+1,:);  
S2_T2=S2(0.6*m+1,:); 
S2_T3=S2(0.8*m+1,:); 
S2_T4=S2(1*m+1,:); 

%Similarly calculated the price of S3 at T1, T2, T3 and T4.
S3_T1=S3(0.4*m+1,:);  
S3_T2=S3(0.6*m+1,:); 
S3_T3=S3(0.8*m+1,:); 
S3_T4=S3(1*m+1,:);
        
for i=1:N
    %Calculating the Worst Performing Stock and dividing by their respective strike price
    S_worst(i)=min([S1_T4(i)/K1 ;S2_T4(i)/K2; S3_T4(i)/K3]);
    
    %Checking if the early redemption condition is met at T1.
    if (S1_T1(i)>=K1)&&(S2_T1(i)>=K2)&&(S3_T1(i)>=K3)
        Payoff(i)=(exp(-r*T1).*(IV))+CP1+CP2;
        
    %Checking if the early redemption condition is met at T2.
    elseif (S1_T2(i)>=K1)&&(S2_T2(i)>=K2)&&(S3_T2(i)>=K3)
        Payoff(i)=(exp(-r*T2).*(IV))+CP1+CP2+CP3; 
        
    %Checking if the early redemption condition is met at T3.
    elseif (S1_T3(i)>=K1)&&(S2_T3(i)>=K2)&&(S3_T3(i)>=K3)
        Payoff(i)=(exp(-r*T3).*(IV))+CP1+CP2+CP3+CP4;

    %Checking if the final redemption condition is met at T4.    
    %If no Barrier Event has occurred.
    elseif Check1(i) == 0 
            Payoff(i)=(exp(-r*T4).*(IV))+CP1+CP2+CP3+CP4+CP5;
    
    %If Barrier Event has occurred.    
    else 
        %The Final Level of each Underlying is at or above its Strike
        if S1_T4(i)>=K1 && S2_T4(i)>=K2 && S3_T4(i)>=K3 
            Payoff(i)=(exp(-r*T4).*(IV))+CP1+CP2+CP3+CP4+CP5;
       
        %The Final Level of at least one Underlying is below its Strike
        else
            Payoff(i)=(exp(-r*T4).*(S_worst(i)*IV))+CP1+CP2+CP3+CP4+CP5;
        end
    end
end

format longg
Payoff_Mean=mean(Payoff)
Payoff_Var=(var(Payoff))/N

%Using Antithetic Method to reduce the variance
T=1.25;
m=100;
dt=T/m;

S1p=S01*ones(m+1,N/2);
S1n=S01*ones(m+1,N/2);
S2p=S02*ones(m+1,N/2);
S2n=S02*ones(m+1,N/2);
S3p=S03*ones(m+1,N/2);
S3n=S03*ones(m+1,N/2);

M1p=zeros(m,N/2);
M1n=zeros(m,N/2);
M2p=zeros(m,N/2);
M2n=zeros(m,N/2);
M3p=zeros(m,N/2);
M3n=zeros(m,N/2);

U1p=unifrnd(0,1,m,N/2);
U1n=-U1p; %Introducing negative dependence between the pairs of replications.
U2p=unifrnd(0,1,m,N/2);
U2n=-U2p; %Introducing negative dependence between the pairs of replications.
U3p=unifrnd(0,1,m,N/2);
U3n=-U3p; %Introducing negative dependence between the pairs of replications.

Z_01=normrnd(0,1,m,N/2);
Z_02=normrnd(0,1,m,N/2);
Z_03=normrnd(0,1,m,N/2);

Z1p=A(1,1)*Z_01;
Z1n=-Z1p; %Introducing negative dependence between the pairs of replications
Z2p=A(2,1)*Z_01+A(2,2)*Z_02;
Z2n=-Z2p; %Introducing negative dependence between the pairs of replications
Z3p=A(3,1)*Z_01+A(3,2)*Z_02+A(3,3)*Z_03; 
Z3n=-Z3p; %Introducing negative dependence between the pairs of replications

Barrier_a1=zeros(m,N/2);
Barrier_a2=zeros(m,N/2);

Payoff_p=zeros(1,N/2);
Payoff_n=zeros(1,N/2);

S_worst_p=zeros(1,N/2);
S_worst_n=zeros(1,N/2);


for n=1:N/2
    for i=1:m
        S1p(i+1,n)=S1p(i,n).*exp((r-d1-0.5*sig1^2)*dt+sig1*sqrt(dt)*Z1p(i,n));
        S1n(i+1,n)=S1n(i,n).*exp((r-d1-0.5*sig1^2)*dt+sig1*sqrt(dt)*Z1n(i,n));
        S2p(i+1,n)=S2p(i,n).*exp((r-d2-0.5*sig2^2)*dt+sig2*sqrt(dt)*Z2p(i,n));
        S2n(i+1,n)=S2n(i,n).*exp((r-d2-0.5*sig2^2)*dt+sig2*sqrt(dt)*Z2n(i,n));
        S3p(i+1,n)=S3p(i,n).*exp((r-d3-0.5*sig3^2)*dt+sig3*sqrt(dt)*Z3p(i,n));
        S3n(i+1,n)=S3n(i,n).*exp((r-d3-0.5*sig3^2)*dt+sig3*sqrt(dt)*Z3n(i,n));
        M1p(i,n)=exp(0.5*(log(S1p(i+1,n)*S1p(i,n))-sqrt((log(S1p(i+1,n)/S1p(i,n)))^2-2*sig1^2*dt*log(U1p(i,n)))));
        M1n(i,n)=exp(0.5*(log(S1n(i+1,n)*S1n(i,n))-sqrt((log(S1n(i+1,n)/S1n(i,n)))^2-2*sig1^2*dt*log(U1n(i,n)))));
        M2p(i,n)=exp(0.5*(log(S2p(i+1,n)*S2p(i,n))-sqrt((log(S2p(i+1,n)/S2p(i,n)))^2-2*sig2^2*dt*log(U2p(i,n)))));
        M2n(i,n)=exp(0.5*(log(S2n(i+1,n)*S2n(i,n))-sqrt((log(S2n(i+1,n)/S2n(i,n)))^2-2*sig2^2*dt*log(U2n(i,n)))));
        M3p(i,n)=exp(0.5*(log(S3p(i+1,n)*S3p(i,n))-sqrt((log(S3p(i+1,n)/S3p(i,n)))^2-2*sig3^2*dt*log(U3p(i,n)))));
        M3n(i,n)=exp(0.5*(log(S3n(i+1,n)*S3n(i,n))-sqrt((log(S3n(i+1,n)/S3n(i,n)))^2-2*sig3^2*dt*log(U3n(i,n)))));
    end
    if M1p(i,n)>B1 && M2p(i,n)>B2 && M3p(i,n)>B3
            Barrier_a1(i,n)=0; %The underlying is above the Barrier level.
        else
            Barrier_a1(i,n)=1; 
    end
    if M1n(i,n)>B1 && M2n(i,n)>B2 && M3n(i,n)>B3
            Barrier_a2(i,n)=0; 
        else
            Barrier_a2(i,n)=1; 
    end    
end

Check1_p=sum(Barrier_a1); 
Check1_n=sum(Barrier_a2); 

S1_T1p=S1p(0.4*m+1,:); 
S1_T1n=S1n(0.4*m+1,:); 
S1_T2p=S1p(0.6*m+1,:); 
S1_T2n=S1n(0.6*m+1,:); 
S1_T3p=S1p(0.8*m+1,:); 
S1_T3n=S1n(0.8*m+1,:); 
S1_T4p=S1p(1*m+1,:); 
S1_T4n=S1n(1*m+1,:); 

S2_T1p=S2p(0.4*m+1,:); 
S2_T1n=S2n(0.4*m+1,:); 
S2_T2p=S2p(0.6*m+1,:); 
S2_T2n=S2n(0.6*m+1,:); 
S2_T3p=S2p(0.8*m+1,:); 
S2_T3n=S2n(0.8*m+1,:); 
S2_T4p=S2p(1*m+1,:); 
S2_T4n=S2n(1*m+1,:); 

S3_T1p=S3p(0.4*m+1,:); 
S3_T1n=S3n(0.4*m+1,:); 
S3_T2p=S3p(0.6*m+1,:); 
S3_T2n=S3n(0.6*m+1,:); 
S3_T3p=S3p(0.8*m+1,:); 
S3_T3n=S3n(0.8*m+1,:); 
S3_T4p=S3p(1*m+1,:); 
S3_T4n=S3n(1*m+1,:); 

for i=1:N/2
    %Calculating the Worst Performing Stock and dividing by their respective strike price
    S_worst_p(i)=min([S1_T4p(i)/K1 ;S2_T4p(i)/K2; S3_T4p(i)/K3]);
    
    %Checking if the early redemption condition is met at T1.
    if (S1_T1p(i)>=K1)&&(S2_T1p(i)>=K2)&&(S3_T1p(i)>=K3)
        Payoff_p(i)=(exp(-r*T1).*(IV))+CP1+CP2;
        
    %Checking if the early redemption condition is met at T2.
    elseif (S1_T2p(i)>=K1)&&(S2_T2p(i)>=K2)&&(S3_T2p(i)>=K3)
        Payoff_p(i)=(exp(-r*T2).*(IV))+CP1+CP2+CP3;
    
    %Checking if the early redemption condition is met at T3.
    elseif (S1_T3p(i)>=K1)&&(S2_T3p(i)>=K2)&&(S3_T3p(i)>=K3)
        Payoff_p(i)=(exp(-r*T3).*(IV))+CP1+CP2+CP3+CP4;
   
    %Checking if the final redemption condition is met at T4.    
    %If no Barrier Event has occurred.
    elseif Check1_p(i) == 0 
            Payoff_p(i)=(exp(-r*T4).*(IV))+CP1+CP2+CP3+CP4+CP5;
    
    %If Barrier Event has occurred.        
    else
        %The Final Level of each Underlying is at or above its Strike
        if S1_T4p(i)>=K1 && S2_T4p(i)>=K2 && S3_T4p(i)>=K3
            Payoff_p(i)=(exp(-r*T4).*(IV))+CP1+CP2+CP3+CP4+CP5;
        else
            %The Final Level of at least one Underlying is below its Strike
            Payoff_p(i)=(exp(-r*T4).*(S_worst_p(i)*IV))+CP1+CP2+CP3+CP4+CP5;
        end
    end
end

for i=1:N/2
    S_worst_n(i)=min([S1_T4n(i)/K1 ;S2_T4n(i)/K2; S3_T4n(i)/K3]);
    
    if (S1_T1n(i)>=K1)&&(S2_T1n(i)>=K2)&&(S3_T1n(i)>=K3)
        Payoff_n(i)=(exp(-r*T1).*(IV))+CP1+CP2;
        
    elseif (S1_T2n(i)>=K1)&&(S2_T2n(i)>=K2)&&(S3_T2n(i)>=K3)
        Payoff_n(i)=(exp(-r*T2).*(IV))+CP1+CP2+CP3;
        
    elseif (S1_T3n(i)>=K1)&&(S2_T3n(i)>=K2)&&(S3_T3n(i)>=K3)
        Payoff_n(i)=(exp(-r*T3).*(IV))+CP1+CP2+CP3+CP4;
    
    elseif Check1_n(i) == 0 
            Payoff_n(i)=(exp(-r*T4).*(IV))+CP1+CP2+CP3+CP4+CP5;

    else 
        if S1_T4n(i)>=K1 && S2_T4n(i)>=K2 && S3_T4n(i)>=K3 
            Payoff_n(i)=(exp(-r*T4).*(IV))+CP1+CP2+CP3+CP4+CP5;
        else
            Payoff_n(i)=(exp(-r*T4).*(S_worst_n(i)*IV))+CP1+CP2+CP3+CP4+CP5;
        end
    end
end

Payoff_Anti=(Payoff_p+Payoff_n)/2; %Average of the positve and negative Payoff.
Payoff_Anti_Mean=mean(Payoff_Anti) %Mean of the Payoff.
Payoff_Anti_Var=(var(Payoff_Anti))*2/N %Variance of the Payoff.