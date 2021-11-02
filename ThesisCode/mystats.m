function [MainA,MainB,Interaction,SimpleAatB,SimpleBatA] = mystats(data,n_factors,lF)
 
% onlinestatbook.com/calculators/tukeycdf.html
% alternatively you can use a table of the studentized distribution: 
% http://cse.niaes.affrc.go.jp/miwa/probcalc/s-range/srng_tbl.html#fivepercent

% IF YOU DONT HAVE THE STATS TOOLBOX (for fcdf), CAN USE THIS WEBSITE
% http://davidmlane.com/hyperstat/F_table.html

% Attempt to do my stats myself for pairwise comparisons
 % works, if same nb of participants in each condition
%clear all
% close all 

% Need initially: save('filename.txt', 'data','-ascii')
 % data=load('data.txt'); % or save the data in a variable called 'data'

levels=[];Y=[];S=[];F1=[];F2=[];F3=[];F01=[];F02=[];F03=[];
n_cond=1;

% n_factors=input('How many factors?  ');
for i=1:n_factors
    levels(i)=lF(i);
    n_cond=n_cond*levels(i);
end

[n_participants,~]=size(data);
for i=1:n_participants
    Y=cat(2,Y,data(i,:));
    S=cat(2,S,zeros(1,n_cond)+i);
end


if n_factors==1
    F01=[1:levels];
    names={'factor1'};
    for i=1:n_participants
        F1=cat(2,F1,F01);
    end
    mik_anova1 
end


if n_factors==2
    for i=1:levels(1)
        F01=cat(2,F01,zeros(1,levels(2))+i);
        F02=cat(2,F02,[1:levels(2)]);
    end
    names={'factor1','factor2'};
    for i=1:n_participants
        F1=cat(2,F1,F01);
        F2=cat(2,F2,F02);
    end
    [MainA,MainB,Interaction,SimpleAatB,SimpleBatA] = mik_anova2(F1,F2,S,Y);
end

if n_factors==3
    for i=1:levels(1)
        F01=cat(2,F01,zeros(1,levels(2)*levels(3))+i);
        for j=1:levels(2)
            F02=cat(2,F02,zeros(1,levels(3))+j);
            F03=cat(2,F03,[1:levels(3)]);
        end
    end
    names={'factor1','factor2','factor3'};
    for i=1:n_participants
        F1=cat(2,F1,F01);
        F2=cat(2,F2,F02);
        F3=cat(2,F3,F03);
    end
    mik_anova3 
end

    % I actually realised that there are already some programs anova in the
    % stats toolbox of Matlab, but those ones, I've created are better
    % since they give the results of ANOVA, the simple main effects (up to
    % 3-way simple main effects) and the pairwise comparisons (up to 3-way)
    % for Scheffe and Tukey. This program only works for WITHIN-SUBJECTS!!
    
