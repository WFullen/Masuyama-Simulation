function [estBias, estRmse, ntrialVec] = p01Masuyama(ntrials,csvfile,r,plotsOn)

% p01Masuyama(ntrials, csvfile, r, plotsOn): this function creates a
% simulation of Masuyama's method; we select a random point on the entire stand
% (that is, we find a point in the range [-r, 750 + 2] for both the x and y
% coordinate; then, we select an input radius r > 0 and record the trees total
% basal area (TBA) within this circle that is generated on the plot stand.
% Note: for this method, each tree has the same probability of selection
% (that is, the probability is equal to the area of the sampled circle over
% the entire area of the stand - (pi*r*r)/( 750 + 2*r )

% INPUTS: ntrials, the number of trials for which we will run our
% simulation (this will give us a good guage of what values the bias / RMSE
% start to converge to
%         csvfile, a common separated values (CSV) file that contains the
% x and y coordinates of the trees, as well as the TBA of each tree
%         r, the radius of the plots for which we will generate our
%         simulated circles and collect the TBA
%         plotsOn, a True/False statement that determines whether or not we
%         display the RMSE and Bias of Estimate plots

% This is the true TBA (we will compare our simulation estimates against
% this)
t = 311.906;

% The following instantiates the three vectors that we will be calculating
% in the simulation; the while loop constructs a vector of length ~ 50, for
% which we will calculate the bias and RMSE estimates from our simulation

k=1;
while ceil(ntrials/k) > 50
    k=k+1;
end

nSample=k;

ntrialVec= nSample:nSample:ntrials; 

estBias = zeros(1, length(ntrialVec));
estRmse = zeros(1, length(ntrialVec));

% This reads in the csv file of the trees data, and parses out three
% things: the x-coordinates, y-coordinates, and the total basal area

trees = csvread(csvfile, 1, 0);
x = trees(:,end-1);
y = trees(:,end);
ba = trees(:,3);

% This instantiates the sample area and total area values - as stated in
% the initial description, all of these values (called pi_i's) are equal.
% That is, each tree has equal probability of being selected

sample_area = pi*r^2;
total_area = (750 + 2*r)^2;
pi_i = sample_area/total_area;

t_hat = zeros(1, ntrials);

counter = 1;

% Below are the steps that are conducted in this simulation:
%   - Calculate two random points between [-r , 750 + r] (as you can see,
%   these points can be outside the 750 x 750 square)
%   - Calculate the distance between our random points and every tree in
%   the uploaded csv file
%   - Sum the Basal Area if the tree is within or on the circle (this will
%   become our total basal area
%   - If the ntrials that the "for" loop is on is divisible by the nSample
%   metric described above, we will plot it (therefore, we will calculate
%   the Bias and RMSE as set out by their respective formulas)

for i = 1:ntrials
    
    a = (750+2*r)*rand() - r;
    b = (750+2*r)*rand() - r;
    
    distance = ( (a - x).^2 + ( b - y).^2 );
    
    y_i = sum( ba( distance <= r^2 ) );

    t_hat(i) = (1/pi_i)*y_i;
    
    if mod(i, nSample) == 0

        estBias(counter) = (100/t)*( (sum(t_hat(1:i)) / i ) - t);
        estRmse(counter) = (100/t)*sqrt( var(t_hat(1:i)) );
        counter = counter + 1;
            
    end
end

% If the plotsOn is 'on', display the plots of the bias and RMSE estimates
% for the Masuyama Method

if strcmp(plotsOn, 'on')
    plot(ntrialVec, estBias, '*b-');
    grid on
    xlabel('Number of Trials');
    ylabel('Percentage Bias of Estimate');
    title(sprintf('Percent Bias of TBA Estimate Using Masuyamas Method \n ntrials = %5.0e, r = %d', ntrials, r));
    figure()
    
    plot(ntrialVec, estRmse, '*b-');
    grid on
    xlabel('Number of Trials');
    ylabel('Percentage RMSE of Estimate');
    title(sprintf('Percent RMSE of TBA Estimate Using Masuyamas Method \n ntrials = %5.0e, r = %d', ntrials, r));
    
end
    
end
