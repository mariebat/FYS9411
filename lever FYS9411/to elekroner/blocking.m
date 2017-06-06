clear 
close all 

filenameIN = 'energies_6el_w1'; 
suffixIN = '.txt'; 
%filenameIN = char(A); 
%suffixIN = char(B);
filename = strcat(filenameIN,suffixIN); 

ID = fopen(filename,'r');
formatspec = '%f'; 
data = fscanf(ID,formatspec); 

nBlocks = 10000; 
%minBlockSize = 100; 
%maxBlockSize = length(data)/100; 
%blockStepLength = round((maxBlockSize-minBlockSize+1)/nBlocks);

n = length(data); 

%nBlocks_vec = [n,floor(n/4),floor(n/2)];

[mean1, var1] = calculateMeanAndVar(data);
std = sqrt(var1/(n-1));

block_sizes = linspace(0,0,nBlocks); 
block_means = linspace(0,0,nBlocks); 
block_vars = linspace(0,0,nBlocks); 
block_std = linspace(0,0,nBlocks); 

count = 1;
for i=1:n
    if mod(n,i) < 200
        
        BS = i;
        startPoint = 1; 
        endPoint = BS; 
        
        %block_results(1,i) = blockSize;
        block_sizes(count) = BS;
        
        %newnBlocks = round(n/blockSize);
        NB = n/BS;  
        %meanvec = linspace(0,0,newnBlocks);
        meanvec = linspace(0,0,NB); 
        
        for j = 1:NB
            if endPoint <= length(data)
            tempvec = data(startPoint:endPoint);
            %tempvec = data(j*BS:(j+1)*BS);
            meanvec(j) = sum(tempvec)/BS;

            startPoint = endPoint;
            endPoint = endPoint + BS;
            end
        end
        
        
        [mean,var] = calculateMeanAndVar(nonzeros(meanvec));
        
        block_means(count) = mean;
        block_vars(count) = var;
        block_std(count) = sqrt(var/((n/BS) -1.0));
        count = count+1;
    disp(i);
    end
end 

%block_sizes = nonzeros(block_Sizes); 
%block_vars = nonzeros(block_Vars); 
%block_means = nonzeros(block_Means); 
%block_std = nonzeros(block_Std); 

newvar = max(block_vars);

figure(01)
semilogx(block_sizes(1:count),block_vars(1:count))
xlabel('Block sizes')
ylabel('Variance')

figure(02)
semilogx(block_sizes(1:count),block_std(1:count))
%axis([0 10^6 0 1e-2])
xlabel('Block size')
ylabel('\sigma')
title('Standard deviation')


x = linspace(1,count,count); 

figure(03) 
plot(x,block_std(1:count))
xlabel('Number of block transformations')
ylabel('standard deviation') 

fileout1 = 'blocking_results_'; 
fileout2 = filenameIN; 
fileout3 = '.mat'; 

fileout = strcat(fileout1,fileout2,fileout3); 

save(fileout, 'block_sizes','block_means','block_vars','block_std'); 
disp('Saved successfully!'); 