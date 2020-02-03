%Amir Kalmoni 
%100987101
%ELEC 4700 Assignment 1

clearvars
clearvars -GLOBAL
close all
format short


set(0, 'DefaultFigureWindowStyle', 'docked')


qo = 1.60217653e-19;              % electron charge
hb = 1.054571596e-34;             % Dirac constant
h = hb * 2 * pi;                  % Planck constant
m_0 = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;               % Boltzmann constant
eps0 = 8.854187817e-12;           % vacuum permittivity
mu0 = 1.2566370614e-6;            % vacuum permeability
c = 299792458;                    % speed of light
g = 9.80665;                      % metres (32.1740 ft) per s²
am = 1.66053892e-27;

effectmass = 0.26 * m_0;
thermvel = sqrt((kb * 300) / effectmass)
standardev = thermvel/(sqrt(2));    %Standard deviation for x and y velocities
timestep = 7.5 * 10 ^ -15;           % time-step value for iteration

width = 200 * 10 ^ -9; % x-boundarie of semiconductor crystal
length = 100 * 10 ^ -9; % y-boundaries of semiconductor crystal

numofelec = 75; %Number of Electrons

temparray = zeros(1, 1000);             %Array of temperature
tmps = (1:1:1000);

%%Question 1 - Assigning random posaitions to the electrons within the
%%boundaries
pos_of_X = rand(1, numofelec) .* width;     %x position
pos_of_Y = rand(1, numofelec) .* length;    %y position


isinbox = true;
while isinbox == true
   inbox = ((pos_of_X <= (1.15 * width/2) & (pos_of_X >= (0.85 * width/2))) & ((pos_of_Y < (length/3)) | pos_of_Y >= (2*length/3)));
   if (sum(inbox) > 0)
       pos_of_X(inbox) = rand(1, sum(inbox)) .* width;
       pos_of_Y(inbox) = rand(1, sum(inbox)) .* length;
   else 
       isinbox = false;
       
   end 
       
end

%%Question 2 - Assigning a Random Velocity to the particles, following a
%Maxwell-Boltsmann distribution
xvel = randn(1, numofelec) .* standardev;
yvel = randn(1, numofelec) .* standardev;
vrms = sqrt((xvel .^ 2) + (yvel .^ 2));
vrmsarray = zeros(1, 1000);

%Calculating the probability of scattering 
pscat = 1 - (exp((-1 * timestep) / (0.2 * 10 ^ -12)));
temp_K = 300;

is2 = zeros(1, numofelec);

numcol = 0;
tdiff = 0;
sumtdiff = 0;

boundtimestep = 0; %Boundary conditions for Question 3

%%Main Iteration loop 1000 times

for i = 1:200
    
    is = pscat > rand(1,numofelec);
    
    xvel(is) = randn .* standardev;
    yvel(is) = randn .* standardev;


    if is(1) ~= is2(1)
        numcol = numcol + 1;
        sumtdfiff = sumtdiff + tdiff;
        tdiff = 0;
    else
        tdiff = tdiff + 1;
       
    end
        
    %%This section models the electron motion specified in the third part
    %of question 1
    pos_of_X(pos_of_X >= width) = pos_of_X(pos_of_X >= width) - width;
    pos_of_X(pos_of_X <= 0) = pos_of_X(pos_of_X <= 0) + width;
    
    ylg = (pos_of_Y >= length);
    ylg1 = (pos_of_Y <= 0);
    
    yvel(ylg) = -yvel(ylg);
    yvel(ylg1) = -yvel(ylg1);
    
    pos_of_XPrev = pos_of_X;
    pos_of_YPrev = pos_of_Y;
   
    %%Question 3 - Enhancements
    %This section deals with the behaviours around the boundaries in the
    %semiconductor. If boundtimestep is equal to 0, the boundaries will be
    %specular, if the boundtimestep is equal to 1, the boundaries are
    %diffusive.
    
    inbox = (((pos_of_X < (1.15 * width/2)) & (pos_of_X > (0.85 * width/2))) & ((pos_of_Y < (length/3)) | pos_of_Y > (2*length/3)));
    
    if ((boundtimestep == 0) && (sum(inbox) >= 1))
        
        if ((pos_of_XPrev < (1.15 * width/2)) & (pos_of_XPrev > (0.85 * (width/2)) & (sum(inbox) >= 1)))
            if (pos_of_Y(inbox) > (2*length/3))
                pos_of_Y(inbox) = pos_of_Y(inbox) - (2 * (pos_of_Y(inbox) - (2*length/3)));
            elseif (pos_of_Y(inbox) < (length/3))
                pos_of_Y(inbox) = pos_of_Y(inbox) + (2 * ((length/3) - pos_of_Y(inbox)));
            end
            yvel(inbox) = -yvel(inbox);
            pos_of_Y(inbox) = pos_of_Y(inbox) + (yvel(inbox) .* timestep);
            pos_of_X(inbox) = pos_of_X(inbox) + (xvel(inbox) .* timestep);
        else
            xvel(inbox) = -xvel(inbox);
            pos_of_Y(inbox) = pos_of_Y(inbox) + (yvel(inbox) .* timestep);
            pos_of_X(inbox) = pos_of_X(inbox) + (xvel(inbox) .* timestep);
        end
        
        vrms = sqrt((xvel .^ 2) + (yvel .^ 2));
        
    elseif ((boundtimestep == 1) && (sum(inbox) >=1))
        
        if ((pos_of_XPrev < (1.15 * width/2)) & (pos_of_XPrev > (0.85 * (width/2)) & (yvel(inbox) > 0)))
            pos_of_Y(inbox) = pos_of_Y(inbox) - (2 * (pos_of_Y(inbox) - (2*length/3)));
            xvel(inbox) = randn .* standardev;
            yvel(inbox) = -1. * (abs(randn .* standardev));
        elseif ((pos_of_XPrev < (1.15 * width/2)) & (pos_of_XPrev > (0.85 * (width/2)) & (yvel(inbox) < 0)))
            pos_of_Y(inbox) = pos_of_Y(inbox) + (2 * ((length/3) - pos_of_Y(inbox)));
            xvel(inbox) = randn .* standardev;
            yvel(inbox) = (abs(randn .* standardev));
        elseif (xvel(inbox) > 0)
            pos_of_X(inbox) = pos_of_X(inbox) - (2 * (pos_of_X(inbox) - (0.85*width/2)));
            xvel(inbox) = -1 .* abs(randn .* standardev);
            yvel(inbox) = randn .* standardev;
        else
            pos_of_X(inbox) = pos_of_X(inbox) + ((2 *(1.15*width/2)) - pos_of_X(inbox));
            xvel(inbox) = abs(randn .* standardev);
            yvel(inbox) = abs(randn .* standardev);
        end
        
        vrms = sqrt((xvel .^ 2) + (yvel .^ 2));
    end
    
    pos_of_X = pos_of_XPrev + (xvel .* timestep);
    pos_of_Y = pos_of_YPrev + (yvel .* timestep);
    
    
    vrms = sqrt((xvel .^ 2) + (yvel .^ 2));
    temp_K = (sqrt(2)*(mean(vrms) ^ 2) * effectmass) / kb;
    temparray(1, i) = temp_K;
    is2 = is;
    
    figure(1);
    plot(pos_of_X, pos_of_Y, 'o');
    xlabel("X axis of Semiconductor Crystal");
    ylabel("Y axis of Semiconductor Crystal");
    title(["Average Temperature = " num2str(temp_K)]);
    
    xlim([0 width]);
    ylim([0 length]);
    hold on
    
    %%This section of code draws the lines of the enhancement
    
    line([0.85*width/2 0.85*width/2], [length 2*length/3]);
    line([1.15*width/2 1.15*width/2], [length 2*length/3]);
    line([0.85*width/2 1.15*width/2], [length length]);
    line([0.85*width/2 1.15*width/2], [2*length/3 2*length/3]);
    
    line([0.85*width/2 0.85*width/2], [0 length/3]);
    line([1.15*width/2 1.15*width/2], [0 length/3]);
    line([0.85*width/2 1.15*width/2], [0 0]);
    line([0.85*width/2 1.15*width/2], [length/3 length/3]);    
    
end

%%plot of temperature versus time
figure (2)
plot(tmps, temparray);
xlabel('time');
ylabel('Average Temperature');
title('Temperature vs Time');
hold on

%  meanfreepath and meanfreetime Calculations
meanfreetime = (sumtdiff * timestep)/numcol;
meanfreepath = mean(vrms) * meanfreetime;

fprintf("The Mean Free Time is = %12f", meanfreetime);
fprintf("The Mean Free Path is = %12f", meanfreepath);

%%Question 3.3 and 3.4

[xgr, ygr] = meshgrid(0:(width/10):width, 0:(length/10):length);
elecmatrx = zeros(11, 11);
temprmatrx = zeros(11, 11);
numelec = 0;
totvelocity = 0;

for e = 1:10
    minx = xgr(1, e);
    maxx = xgr(1, e+1);
    for a = 1:10
        miny = ygr(a, 1);
        maxy = ygr(a+1, 1);
        for kk = 1:numofelec
            if((pos_of_X(kk) > minx) && (pos_of_X(kk) < maxx) && ((pos_of_Y(kk) > miny) && pos_of_Y(kk) < maxy))
                numelec = numelec + 1;
                elecmatrx(e, a) = elecmatrx(e, a) + 1;
                totvelocity = totvelocity + sqrt((xvel(kk) .^ 2) + (yvel(kk) .^ 2));
                temprmatrx(e, a) = ((sqrt(2)*(totvelocity/numelec) ^ 2) * effectmass) / kb;
            end
        end
        totvelocity = 0;
        numelec = 0;
    end
end

%%Question 3.3 and 3.4 - Creating the plots of the histogram from question
%2, and the electron density map, along with the temperature map.

figure(3); hist(vrms, 10); 
title('Histogram of Thermal Velocities');

figure(4); surf(elecmatrx);
title('Electron Density Map');

figure(5); surf(temprmatrx);
title('Temperature Map');
