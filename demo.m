clear all

addpath(genpath('f_math'));
addpath(genpath('f_em_doa_mimo'));
addpath(genpath('f_tensor_esprit'));
addpath(genpath('f_somp'));
addpath(genpath('f_figure'));

%% System Info
MIMO_info.Nt = 64; % transmitting antennas
MIMO_info.P = 16; % number of training frames
MIMO_info.Nr = 32; % receiving antennas
MIMO_info.Q = 6; % number of combiner output

MIMO_info.K_0 = 256; % total number of subcarriers
MIMO_info.K = 16; % number of training subcarriers
MIMO_info.K_select = 0 : floor(MIMO_info.K_0/MIMO_info.K) : floor(MIMO_info.K_0/MIMO_info.K)*(MIMO_info.K-1); % subcarriers for training: 0, K_0/K, ..., (K-1)K_0/K

MIMO_info.f_c = 60*1e9; % carrier frequency 60GHz
MIMO_info.f_s = 1.76*1e9; % bandwidth 1.76GHz, 256 IFFT, 16 subbands used for training, corresponding Ts = 145.5 ns

% =========================================================================
% Choice 1:
% Identity matrix, do not support for P>Nr or Q > Nr
% It generally have better performance, as
% 1. for tensor esprit: it only works when the vandermonde structure is retained
% 2. for cs methods: it makes the correlation of the dictionary small
% However, why it has a higher CRB than random matrices?
% -------------------------------------------------------------------------
MIMO_info.F = zeros(MIMO_info.Nt,MIMO_info.P); % precoding matrix
MIMO_info.F(1:MIMO_info.P,:) = 1./sqrt(MIMO_info.Nt).*eye(MIMO_info.P); % then each antenna transmits signal power 1/N_t, and SNR defined as 1/sigma^2 or 1/N_t/sigma^2?
MIMO_info.W = zeros(MIMO_info.Nr,MIMO_info.Q); % combining matrix
MIMO_info.W(1:MIMO_info.Q,:) = 1./sqrt(MIMO_info.Nr).*eye(MIMO_info.Q);
% =========================================================================
% Choice 2:
% Random matrix
% -------------------------------------------------------------------------
% MIMO_info.F = 1./sqrt(MIMO_info.Nt).*exp(1i * 2 * pi/4 * randi(4,MIMO_info.Nt,MIMO_info.P));
% MIMO_info.W = 1./sqrt(MIMO_info.Nr).*exp(1i * 2 * pi/4 * randi(4,MIMO_info.Nr,MIMO_info.Q));
% =========================================================================
%% Channel built
Channel_info.L = 4; % number of path
Channel_info.SNR = 10; % \dB

Channel_info.alpha = 1/sqrt(2)*(   randn(Channel_info.L,1) + 1i*randn(Channel_info.L,1)   ); % path gain: CN(0,1)
Channel_info.tau = 1e-7 * rand(Channel_info.L,1);
Channel_info.phi = -pi/2 + pi/256 + pi/128*[40;60;80;110] + pi/200;
Channel_info.theta = -pi/2 + pi/256 + pi/128*[106;68;35;82] + pi/200;

[Yn,Ytrue,H] = channel_build(MIMO_info,Channel_info);
clf
stem3(   sin(Channel_info.theta)/2,sin(Channel_info.phi)/2,zeros(Channel_info.L,1),'*r'   ); hold on
stem3(   sin(Channel_info.theta)/2,sin(Channel_info.phi)/2,abs(Channel_info.alpha),'*r'   )
hold on;

%% Estimation
[Channel_est,Y_est] = em_offgrid_dualwideband(Yn,MIMO_info,...
                                        'show_info',true);


%% Evaluation
Hmmse = (   Channel_est.H(:)-H(:)   )'*(   Channel_est.H(:)-H(:)   )/(   H(:)'*H(:)   )
Ymmse = (   Y_est(:)-Ytrue(:)   )'*(   Y_est(:)-Ytrue(:)   )/(   Ytrue(:)'*Ytrue(:)   )
