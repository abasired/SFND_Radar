clear all
close all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Defined Range and Velocity of target
% 
% Define the target's initial position and velocity. 
% Note : Velocity remains contant
 
R = 110;
v = 20;


%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the design requirements.


%Operating carrier frequency of Radar 
fc= 77e9;                   %carrier freq
c= 3e8;                     %speed of light
range_res = 1;              %range resolution
R_max = 200;                %max range
B = c/(2*range_res);        %bandwidth of the chirp
Tchirp = 5.5*(2*R_max)/c;   %duration of single chirp
slope = B/Tchirp;           %slope of chirp


                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v*t(i);
    
    
    %For each time sample we need update the transmitted and
    %received signal. 
    td(i) = 2*r_t(i)/c;
    delay = t(i) - td(i);
    Tx(i) = cos(2*pi*(fc*t(i) + (slope*(t(i)^2)/2)));
    Rx(i) = cos(2*pi*(fc*delay + (slope*(delay^2)/2)));
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);
    
end

%% RANGE MEASUREMENT


%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr,Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(Mix,Nr)/Nr;

% Take the absolute value of FFT output
signal_fft = abs(signal_fft);


% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft  = signal_fft(1:Nr/2);

%L = Tchirp*B;

% plot FFT output 
figure ('Name',' range FFT')
subplot(2,1,1)
%Fs = Nr/Tchirp;
f = B*(0:(Nr/2 - 1))/Nr;
%f = Fs*(0:(L-1))/L;
plot(f/1e6,signal_fft) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (MHz)')
ylabel('|P1(f)|')
%axis ([0 100 0 1]);

%plotting the range
subplot(2,1,2)
R = c*Tchirp*f/(2*B);
plot(R,signal_fft)
title('Range computation')
axis ([0 200 0 0.5]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map
%Select the number of Training Cells in both the dimensions.
Tr = 15;
Td = 4;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;


% offset the threshold by SNR value in dB
offset = 15;

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);
signal_cfar = zeros(Nr/2,Nd);
Threshold_cfar = zeros(Nr/2,Nd);

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing

   % CFAR
   
%%

for i = Tr+Gr+1:Nr/2-(Gr+Tr)
    for j = Td+Gd+1:Nd-(Gd+Td)
        noise_level = 0;
        num_tr_cells = 0;
        signal = db2pow(RDM(i,j));
        %CUT (i,j)
        for p = i-(Tr+Gr):i+Tr+Gr
            for q = j-(Td+Gd):j+Td+Gd
                %Training/Gaurd cells (p,q)
                if (abs(i-p)>Gr || abs(j-q)>Gd)
                    %Training cells (p,q)
                    num_tr_cells = num_tr_cells + 1; 
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        Threshold_cfar(i,j) = pow2db(noise_level/num_tr_cells) + offset;
        if signal > db2pow(Threshold_cfar(i,j))
            signal_cfar(i,j) = 1;
        end
            
    end
end
%%
%                     
% for i = 1: Tr+Gr
%     %Top left corner
%     for j = 1:Td+Gd
%         %CUT (i,j)
%         noise_level = 0;
%         num_tr_cells = 0;
%         signal = db2pow(RDM(i,j));
%         
%         for p = 1:i+Tr+Gr
%             for q = 1:j+Td+Gd
%                 %Training/Gaurd cells (p,q)
%                 if (abs(i-p)>Gr || abs(j-q)>Gd)
%                     %Training cells (p,q)
%                     num_tr_cells = num_tr_cells + 1;
%                     noise_level = noise_level + db2pow(RDM(p,q));
%                 end
%             end
%         end
%         Threshold_cfar(i,j) = pow2db(noise_level/num_tr_cells) + offset;
%         if signal > db2pow(Threshold_cfar(i,j))
%             signal_cfar(i,j) = 1;
%         end
%     end
%     %Top right corner
%     for j = Nd-(Gd+Td)+1:Nd
%         %CUT (i,j)
%         noise_level = 0;
%         num_tr_cells = 0;
%         signal = db2pow(RDM(i,j));
%         for p = 0:i+Tr+Gr
%             for q = j-(Td+Gd):j
%                 %Training/Gaurd cells (p,q)
%                 if (abs(i-p)>Gr || abs(j-q)>Gd)
%                     %Training cells (p,q)
%                     num_tr_cells = num_tr_cells + 1;
%                     noise_level = noise_level + db2pow(RDM(p,q));
%                 end
%             end
%         end
%         Threshold_cfar(i,j) = pow2db(noise_level/num_tr_cells) + offset;
%         if signal > db2pow(Threshold_cfar(i,j))
%             signal_cfar(i,j) = 1;
%         end
%     end
% end
% 
% for i = Nr/2 - (Gr+Tr)+1: Nr/2
%     %Bottom left corner
%     for j = 1:Td+Gd
%         %CUT (i,j)
%         noise_level = 0;
%         num_tr_cells = 0;
%         signal = db2pow(RDM(i,j));
%         for p = i-(Tr+Gr):i
%             for q = 1:j+Td+Gd
%                 %Training/Gaurd cells (p,q)
%                 if (abs(i-p)>Gr || abs(j-q)>Gd)
%                     %Training cells (p,q)
%                     num_tr_cells = num_tr_cells + 1;
%                     noise_level = noise_level + db2pow(RDM(p,q));
%                 end
%             end
%         end
%         Threshold_cfar(i,j) = pow2db(noise_level/num_tr_cells) + offset;
%         if signal > db2pow(Threshold_cfar(i,j))
%             signal_cfar(i,j) = 1;
%         end
%     end
%     %Bottom right corner
%     for j = Nd-(Gd+Td)+1:Nd
%         %CUT (i,j)
%         noise_level = 0;
%         num_tr_cells = 0;
%         signal = db2pow(RDM(i,j));

%         for p = i-(Tr+Gr):i
%             for q = j-(Td+Gd):j
%                 %Training/Gaurd cells (p,q)
%                 if (abs(i-p)>Gr || abs(j-q)>Gd)
%                     %Training cells (p,q)
%                     num_tr_cells = num_tr_cells + 1;
%                     noise_level = noise_level + db2pow(RDM(p,q));
%                 end
%             end
%         end
%         Threshold_cfar(i,j) = pow2db(noise_level/num_tr_cells) + offset;
%         if signal > db2pow(Threshold_cfar(i,j))
%             signal_cfar(i,j) = 1;
%         end
%     end
% end
                        
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;


 
 