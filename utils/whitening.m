function [xnew] = whitening(x, Fs, varargin) 
%   SPECTRAL WHITENING
%
%   This function generates flat Fourier spectrum for a given signal (which
%   is originally not white) either for the full range of 0 Hz to the
%   Nyquist frequency or for a user defined frequency band. This operation
%   tends to sharpen signal, as well as the noise. The whitening process is
%   often used for ambient vibration data before stacking waveforms for
%   cross-correlation. The whitening process is often used for ambient
%   vibration data before stacking waveforms for cross-correlation. The
%   process is simple as Fourier transforming the signal after applying
%   Hann window, then normalizing its magnitude, and then inverse Fourier
%   transforming it.
%
%   Syntax:    
%       xnew = whitening(x, Fs, freq, []) for the full range of 0 Hz to the
%              Nyquist frequency
%       xnew = whitening(x, Fs, freq, [0.1, 20]) for 0.1 Hz to 20 Hz as 
%              an example
%
%    Input: 
%          x = input signal (x must be a row vector)
%         Fs = sampling rate (e.g., 200)
%       freq = frequency limit for whitening in Hz (e.g., freq, [] or 
%              freq, [0.1, 20])
%
%   Output: 
%       xnew = spectrally whitened signal for the full range of 0 Hz to the
%              Nyquist frequency or within a user defined frequency band
%
%   Example: See demo.m file 
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 1.0 $  $Date: 2017/12/08 09:00:00 $
%
%% Input parser 
p = inputParser();
p.CaseSensitive = false;
p.addOptional('freq',[]); % prescribed or not prescribed frequencies (empty or vector)
p.parse(varargin{:});

% shorthen the variable name
freq = p.Results.freq;

% check if freq is empty
if isempty(freq)
%     warning('freq is specified as empty. The default value is taken as 0 Hz to Nyquist');
    freq = [0, Fs/2];
end
L = numel(x);

% L-point symmetric Hann window in the column vector W
W = hann(L); 

% Multiply input signal, x, with Hann window and take FFT in
% order to make sure that the ends of signal in frequency domain match up
% while keeping everything reasonably smooth; this greatly diminishes
% spectral leakage
xf = fft(W.*x,L); 

% Magnitude vector   
mag = abs(xf); 

% Phase vector
phase = unwrap(angle(xf)); 

% Frequency vector
f = (0:L-1)*(Fs/L);

% Find lower and upper indices in f vector
[~,index(1)] = min(abs(f-freq(1))); 
[~,index(2)] = min(abs(f-freq(2))); 

mag(index(1):index(2)) = 1.0;
xnew = ifft((mag.*exp(1i*phase)),L,'symmetric');

return
