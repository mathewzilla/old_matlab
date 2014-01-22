% AR_train_temp    Store templates for later classification. Modified for
%   publication in Evans et al 2014
%   'State of the art in whisker-based object localisation'.

%   out = AR_train_temp(data_train) for DATA an array.
%   Outputs an array of templates. Not very clever! This is basically for
%   consistency and timing.
%   Old code executed all kinds of preprocessing, but now that is handled
%   by AR_preprocessing.m
%   An option for template generation method is included, e.g. for
%   averaging over train data instead of taking a single sample
%

function templates = AR_train_temp(data_train,method)

% Included option for template 'method'
% 1 = low pass filtered
% 2 = average template
% 3 = Velocity
% 4 = Velocity of single trial
% 5 = Frequency template

[nt nn] = size(data_train); 


%if nnd==1; d = repmat(d,[1,nn]); end

%% Low pass filter
if method == 1; 
    f = 1000;                % Sample rate. Should probably have this
    % as input
    
    % Low pass filter the samples
    fNorm = 20/(f/2);        % Low pass cutoff at 20Hz. ADD OPTION
    [b a] = butter(10, fNorm, 'low');
    
    for i = 1:nn;
        %   clear x;
        clear y;
        for j = 2;  % FIXED TO SINGLE TRIAL
            %   y(1,:) = data{j,i};
            y(1,:) = filter(b,a, data_train{j,i});
        end
        
        templates(i,:) = y;
        
    end
end

%% Average template
if method == 2;
    
    for i = 1:nn;
        clear x;
        clear y;
        for j = 1:3;  %:nt;
            x(j,:) = data_train{j,i};
        end
        
        for k = 1:length(x);
            y(k) = mean(x(:,k));
        end
        
        templates(i,:) = y;
        
    end
end

%% First derivative of average template REDUNDANT?
if method == 3;
    

    for i = 1:nn;
        clear x;
        clear y;
        for j = 1:3;  %:nt;
            x(j,:) = data_train{j,i};
        end
        
        for k = 1:length(x);
            y(k) = mean(x(:,k));
        end
        
        y = diff(y);
        
        templates(i,:) = y;
        
    end
end

%% First derivative of single trial REDUNDANT?
if method == 4;
    
    for i = 1:nn;
        clear x;
        clear y;
        for j = 2;  %:nt;
            x = data_train{j,i};
        end
        
        y = diff(x);
        
        if smoo;
            y = smooth(y);
        end
        templates(i,:) = y;
        
    end
end

%% Frequency spectrum template
if method == 5;

    for i = 1:nn;
        clear x;
        clear y;
        for j = 2;
            x = abs(real(fft(data_train{j,i}(501:750)))); 
            % Only looks at the offset here, must have been a one-off test 
            % to see what sticks
        end
        
        x(1) = 0;
%         if smoo;
            x = smooth(smooth(x));
%         end
        templates(i,:) = x;
    end
end
