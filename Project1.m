%% Constants 
clear
clc
close all
A0 = 750000;
r = 0.03;

%% Question 1
timesPerYear = [1 2 4 12];
A = [];
k = 1;
figure(1)
hold on 
for n = timesPerYear
    for t = 0:5
            a = A0 * ( 1+ r/n) ^ (n*t);
            A(k,t+1) = a;
    end
    plot(0:5,A(k,:))
    k = k+1;
end
title("Compounded Interest Comparison");
xlabel("Time (Years)");
ylabel("Outstanding Balance (Dollars)");
legend("Compounded Anually","Compounded Semi-Anually", "Compounded Quarterly","Compounded Monthly","Location", 'northwest');
hold off
figure(2)
hold on
k = 1;
A1 = [];
for j = [ 4 12]
    for ti = 0:30
        a1 = A0 * ( 1+ r/j) ^ (j*ti);
        A1(k,ti+1) = a1;
    end
    plot(0:30, A1(k,:))
    k = k+1;
end
for tim = 0:30
    a2 = A0*exp(r*tim);
    A1(3, tim+1) = a2;
end
plot(0:30, A1(3,:))
title("Compounded Interest Comparison");
xlabel("Time (Years)");
ylabel("Outstanding Balance (Dollars)");
ylim([750000 2000000]);
legend("Compounded Quarterly","Compounded Monthly","Compounded Continuously","Location", 'northwest');
hold off

%% Question 2
%EQulibrium solutions at 
%A = 12p/r

%Testing Constants 
R = 0.05;
P = 100;
Atest = (12 * P / R ) - 1;
APrim = R * Atest - 12* P;
disp(APrim)

R1 = 0;
P1 = 0;

for R = linspace(0.001,0.1)
    R1 = R1 +1;
    for P = linspace(10,200)
        P1 = P1 + 1;
        Atest1 = (12 * P / R ) - 1;
        APrim1 = R * Atest1 - 12* P;
        Atest2 = (12 * P / R ) + 1;
        APrim2 = R * Atest2 - 12* P;
        if (APrim1 > 0) && (APrim2 < 0)
            Stabils(R1,P1) = 1;
            Stab1 = "Stable";
        elseif (APrim1 > 0) && (APrim2 > 0)
            Stabils(R1,P1) = 2;
            Stabl = "Semi";
        elseif (APrim1 < 0) && (APrim2 < 0)
            Stabils(R1,P1) = 2;
            Stabl = "Semi";
        else 
            Stabils(R1,P1) = 3;
            Stabl = "Un";
        end
    end
    P1 = 0;
end

% 1 is stable 
% 2 is semistable 
% 3 is unstable

[rows, cols, StableSolutions] = find(Stabils == 3);

% All the values of the solutions are going to be unstable solutions. We
% were able to determine this by checking the values above and below the
% stability solution for 100 terms of R and 100 terms of P. 
%The equilibria represent where you are only paying off the interest on the
%mortgages and not in fact paying off any of the initial loan 

%% Question 3

%see hand written work

%% Question 4
format bank

r = 0.03;
t = 10;
P10year = -A0*r*exp(r*t)/(12*(1-exp(r*t)));

r = 0.05;
t = 30;
P30year = -A0*r*exp(r*t)/(12*(1-exp(r*t)));

fprintf("Payment per month for a 10 year and 30 year mortgage respectivly \n")
disp(P10year)
disp(P30year)

%% Question 5
Total10Yr = P10year * 10 * 12;
Total30Yr = P30year * 30 * 12;

fprintf("Total Payments for 10 yr and 30 yr mortgages respectivly \n")
disp(Total10Yr)
disp(Total30Yr)

Interest10Year = Total10Yr - A0;
Interest30Year = Total30Yr - A0;

fprintf("Total Interests for 10 yr and 30 yr mortgages respectivly \n")
disp(Interest10Year)
disp(Interest30Year)

%% Question 6
A01 = 650000;
r = 0.03;
t = 10;
P10yearNew = -A01*r*exp(r*t)/(12*(1-exp(r*t)));

r = 0.05;
t = 30;
P30yearNew = -A01*r*exp(r*t)/(12*(1-exp(r*t)));

fprintf("Payment per month for the new 10 year and 30 year mortgage respectivly \n")
disp(P10yearNew)
disp(P30yearNew)

Interest10YearNew = (P10yearNew * 10 * 12) - A0;
Interest30YearNew = (P30yearNew * 30 * 12) - A0;

fprintf("Total Interests for 10 yr and 30 yr mortgages respectivly with downpayment \n")

disp(Interest10YearNew)
disp(Interest30YearNew)

fprintf("Difference in interest paid, ie. savings \n ")
disp(Interest10Year - Interest10YearNew)
disp(Interest30Year - Interest30YearNew)

%% Question 7 
%The advantages to having a 30 year mortgage is that you will be paying
%considerably less per month on the scale of almost half of what you would
%have to pay on a 10 year motgage. This however comes with its downside as
%it results in far more interest paid on the initial amount. Depending upon
%the initial loan amount and whether or not there is a downpayment this can
%be anywhere between 6 and nearly 250 times as much interest that would be
%paid on a 10 year mortgage. 
%This then translates to the 10 year mortgage as follows. The 10 year
%mortgage has a much higher monthly payment ($7234/month vs $4022/month
%with no downpayment and $6269/month vs $3486/month with a $100k
%downpayment) but subsequntly has a much lower total interest paid given
%that the total amount of money decreases faster and as such the interest
%has less time to accrew. The interest paid on a 10 year mortgage is only
%$118k versus $698k on a 30 year mortgage with no downpayment, which then
%decreases drastically when an initial downpayment of 100k is paid. At the
%new initial mortgage amount with the $100k downpayment the 10 year
%mortgage interest is only $2367 while the 30 year mortgage interest total
%remains at $505k. 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION 3.2 

%% Question 1
%Eulers method with the following:
h = 0.5;
r = 0.05;
p = 4000; 
A = zeros();
A(1) = 750000;

n = 62;

for i = 1:n
    Aprime = r*A(i) - 12*p;
    A(i+1) = A(i) + (Aprime * h);
end

%% Question 2
years = (0:n);
Atrue = (12 * p)/r + (A0 - (12*p)/r).*(exp(r.*(years)./2));

figure(4)
hold on 
title("step size of 0.5")
plot((1:n+1)/2,A, 'g')
plot((1:n+1)/2,Atrue, 'r')
yline(0)
xlabel("Time (Years)");
ylabel("Outstanding Balance (Dollars)");
legend("Eulers method values", "Truth values")
hold off

%% Question 3
h1 = 0.01;
A1 = zeros();
A1(1) = 750000;
n1 = 3100;
divisor = 1/h1;

for i = 1:n1
    Aprime = r*A1(i) - 12*p;
    A1(i+1) = A1(i) + (Aprime * h1);
end

years1 = (0:n1);
Atrue1 = (12 * p)/r + (A0 - (12*p)/r).*(exp(r.*(years1)./divisor));

figure(5)
hold on 
title("Step Size of 0.01")
plot((1:n1+1)/100,A1, 'g', 'LineWidth', 3)
plot((1:n1+1)/100,Atrue1, 'r')
xlabel("Time (Years)");
ylabel("Outstanding Balance (Dollars)");
yline(0)
legend("Eulers method values", "Truth values", 'Location', 'SW')
hold off

figure(6)
title("Errors")
Error = Atrue - A;
Error1 = Atrue1 - A1;
hold on 
plot(1:n+1, Error)
plot((1:n1+1)/50,Error1)
xlabel("Time (Years)");
ylabel("Error in Balance (Dollars)");
legend("Error in step size 0.5", "Error is step size 0.01", "Location","Southwest")
hold off

%Figure 6 outlines the difference in error or the difference between the
%Euler approximation and the truth values for each term. While the plot
%itself is somewhat harder to interpret because of the horizontal scale it
%is very informative on the vertical scale. It shows how at the calculated
%values where the mortgage is paid off the error using a smaller step size
%is less than $500 off while the error for the larger step size is more
%than $18k. This outlines the value of using a smaller step size on a
%numberical approximation such as Euler's method. 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION 3.2.2

%% Question 1
htest = 0.01;
Ptest = 4000;
Atest = zeros();
Atest(1) = 750000;

maxnums = 10;
while Atest(end) >= 0
    for timestep = 1:maxnums 
        if timestep <= 500
            rtest = 0.03;
        elseif timestep > 500
            rtest = 0.03 + 0.015* sqrt((timestep*htest)-5);
        end
        
            AprimeTest = rtest*Atest(timestep) - 12*Ptest;
            Atest(timestep+1) = Atest(timestep) + (AprimeTest * htest);
    end 
    maxnums = maxnums +1;
end
            
figure(7)
hold on 
plot((1:maxnums)/100, Atest)
yline(0);
hold off

timeToPay = length(Atest)* 0.01;
fprintf("The time to pay off the Mortgage with the new variable rate is %.2f years\n", timeToPay)
            
%% Question 2
htest = 0.01;
Ptest1 = 4500;
Atest1 = zeros();
Atest1(1) = 750000;

maxnums = 10;
while Atest1(end) >= 0
    for timestep = 1:maxnums 
        if timestep <= 500
            rtest = 0.03;
        elseif timestep > 500
            rtest = 0.03 + 0.015* sqrt((timestep*htest)-5);
        end
        
            AprimeTest1 = rtest*Atest1(timestep) - 12*Ptest1;
            Atest1(timestep+1) = Atest1(timestep) + (AprimeTest1 * htest);
    end 
    maxnums = maxnums +1;
end
            
figure(8)
hold on 
plot((1:maxnums)/100, Atest1)
yline(0);
hold off

timeToPay1 = length(Atest1)* 0.01;
fprintf("The time to pay off the Mortgage with the new variable rate is %.2f years\n", timeToPay1)

%% Question 3
VarInterest = Ptest*timeToPay*12 - A0;
VarInterest1 = Ptest1*timeToPay1*12 -A0;
fprintf("The interest paid on $4000/month will be $%.2d while the interest paid on the $4500/month will be $%.2d", VarInterest, VarInterest1)

%% Question 4
figure(9)
title("Variable Rate Payment Plans")
hold on 
plot((1:length(Atest))/100, Atest)
plot((1:length(Atest1))/100, Atest1)
yline(0)
xlabel("Time (Years)");
ylabel("Outstanding Balance (Dollars)");
hold off
legend("$4000/month", "$4500/month")
%When the monthly payment is $4000 the variable rate payment causes the
%length to pay off to increase by about 4.5 years but when the monthly
%payment is $4500 it will decrease the amount being paid by nearly 8 years.
%This means it will be advantagious for the borrower to take a variable
%rate if they can pay $4500 but not if they can only pay $4000. 
