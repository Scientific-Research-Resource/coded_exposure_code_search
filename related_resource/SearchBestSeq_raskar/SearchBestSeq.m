
% Copyright 2010 Mitsubishi Electric Research Laboratories All Rights Reserved.
% 
% Permission to use, copy and modify this software and its documentation without fee for educational, research and non-profit purposes, is hereby granted, provided that the above copyright notice and the following three paragraphs appear in all copies.
% 
% To request permission to incorporate this software into commercial products contact: Vice President of Marketing and Business Development; Mitsubishi Electric Research Laboratories (MERL), 201 Broadway, Cambridge, MA 02139 or .
% 
% IN NO EVENT SHALL MERL BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF MERL HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
% 
% MERL SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND MERL HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS OR MODIFICATIONS.



function [bNoPad,sminimum] = SearchBestSeq(n,k,objectwidth, max_exhuasive_search_times, max_random_search_times)
% Params:
%	n = code length
%	k = number of ones in the code
%	objectwidth = zero padding for the code. 
%	max_exhuasive_search_times: maximum exhuasive search times
%	max_random_search_times: maximum random search times
% Return:
%	bNoPad: searched code
%	sminimum: spectral minimum
% 
% zzh note: 
%	1. the param 'objectwidth` should be given a relative % large number, like n*8, or the result might  be incorrect
%	2. if the total number of possible solutions < max_exhuasive_search_times, use
%	exhausive searching by matlab; otherwise, use random searching by
%	calling 'SearchBestSeqRandom.exe' (maximum random searching times =
%	max_random_search_times)

if nargin<5
	max_random_search_times = 5*10^6;
end
if nargin<4
	max_exhuasive_search_times = 5*10^5; % the param 'objectwidth` should be given a relative large number
end
if nargin<3
	objectwidth = min(16*n, 2^20); % the param 'objectwidth` should be given a relative large number
end

OPTIMIZE_F = 0;

if(OPTIMIZE_F)
    sminimum = Inf;
else
    sminimum = -1;
end



tt = nchoosek(n-2,k-2);

fprintf('Total permutations = %d\n',tt);


tic

if(tt < max_exhuasive_search_times)

    fprintf('Doing Permutation. Total number = %d\n',tt);

    idx = [2:n-1]';
    % find all permutations
    C = nchoosek(idx,k-2);

    nn = size(C,1);

    if(nn~=tt)
        disp('error cannot come here');
        keyboard
    end



    for i = 1:nn

        bNoPad = zeros(n,1);
        bNoPad(1) = 1;
        bNoPad(end) = 1;

        id = C(i,:);
        bNoPad(id) = 1;



        sumbnopad = sum(bNoPad(:));
        if(sumbnopad~=k)
            disp('error cannot come here');
            keyboard
        end

        E = bNoPad/sumbnopad;
        E = E(:);

        if(OPTIMIZE_F)
            A = motionblurmatrix(E,objectwidth);
            A = full(A');
            fminimum = sqrt(mean(diag(inv(A'*A))));
            clear A
        else
            f = abs(fft([E;zeros(objectwidth,1)]));
            fminimum = min(f(:));
        end



        if(OPTIMIZE_F)
            if(fminimum < sminimum)   %min FFT, we got an improvement in mask
                sminimum = fminimum ;
                bNoPadSaved = bNoPad;
                fprintf('n=%d, k=%d, minFofAprimeA=%f, seq=%s\n',n,k,sminimum,char(bNoPad+48))
            end

        else


            %disp(sprintf('%s',char(bNoPad+48)))

            if(sminimum < fminimum)   %min FFT, we got an improvement in mask
                sminimum = fminimum ;
                bNoPadSaved = bNoPad;
                fprintf('n=%d, k=%d, minFFT=%f, seq=%s\n',n,k,sminimum,char(bNoPad+48))
            end
        end

    end


    bNoPad = bNoPadSaved;

else

    %maximum numbers of search
    startnum = 1
    endnum = max_random_search_times
    jump = 1

    tic
    str = sprintf('SearchBestSeqRandom %d %d %d %d',n,k,objectwidth+n,endnum)
    system(str);
    toc





    str = sprintf('BestSeq_n=%dk=%dT=%d.txt',n,k,objectwidth+n);
    x = load(str);


    a = dec2bin(x(end,4));
    bNoPad = double(a)-48;
    bNoPad = bNoPad(:);
    sminimum = x(end,3);



end


toc













