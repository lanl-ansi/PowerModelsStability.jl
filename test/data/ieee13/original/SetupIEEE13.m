% Changhong Zhao, 01/12/2017, convert IEEE 13 bus test feeder to a model
% that is convenient for an SDP formulation of OPF

jay = sqrt(-1);

Zbase = 1;    % Ohm
Vbase = 4.16;  % in kV.
Sbase = ((Vbase*1e3)^2)/Zbase/1e3; %in kVA
Ybase = 1/Zbase; % Siemens (Ohm^-1)

Nnode = 11;
Nline = Nnode - 1;

BusNames = [{'sourcebus'} {'632'} {'633'} {'634'} {'645'} {'646'} {'671'}  {'684'} {'611'} {'652'} {'675'}];
% ignore no-load bus 680 (always has the same voltage as 671)
% merge bus 692 into 671 (consolidate the loads) since they are connected
% by a closed switch
% sourcebus = 650
BusPhases = [1 1 1;
    1 1 1;
    1 1 1;
    1 1 1;
    0 1 1;
    0 1 1;
    1 1 1;
    %1 1 1;
    1 0 1;
    0 0 1;
    1 0 0;
    %1 1 1;
    1 1 1;
    ];


overline_V = 1;
V_sub = [overline_V  ;
    overline_V*exp(jay*(-120)/180*pi);
    overline_V*exp(jay*120/180*pi)];
%substation bus index
substation=1;

% real power loads in kW
% bus order consistent with network Y matrix and BusNames
P_Delta_l = [0 0 0;   % 650
       0 0 0; %632
       0 0 0;  % 633
       0 0 0;  % 634
       0 0 0; % 645
       0 230 0; % 646
       385 385 385 + 170; % 671
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 0; % 611
       0 0 0; % 652
       %0 0 170;  % 692
       0 0 0; % 675
       ];

% reactive power loads in kVar
Q_Delta_l = [0 0 0;   % 650
       0 0 0; %632
       0 0 0;  % 633
       0 0 0;  % 634
       0 0 0; % 645
       0 132 0; % 646
       220 220 220+151; % 671
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 0; % 611
       0 0 0; % 652
       %0 0 151;  % 692
       0 0 0; % 675
       ];

DeltaLoadPhases = [0 0 0;   % 650
       0 0 0; %632
       0 0 0;  % 633
       0 0 0;  % 634
       0 0 0; % 645
       0 1 0; % 646
       1 1 1; % 671
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 0; % 611
       0 0 0; % 652
       %0 0 170;  % 692
       0 0 0; % 675
       ];


P_Y_l = [0 0 0;   % 650
       8.5 33 58.5; %632    % include half distributed load on the connected line
       0 0 0;  % 633
       160 120 120;  % 634
       0 170 0; % 645
       0 0 0; % 646
       8.5 33 58.5; % 671   % include half distributed load on the connected line
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 170; % 611
       128 0 0; % 652
       %0 0 0;  % 692
       485 68 290; % 675
       ];


Q_Y_l = [0 0 0;   % 650
       5 19 34; %632    % include half distributed load on the connected line
       0 0 0;  % 633
       110 90 90;  % 634
       0 125 0; % 645
       0 0 0; % 646
       5 19 34; % 671   % include half distributed load on the connected line
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 80; % 611
       86 0 0; % 652
       %0 0 0;  % 692
       190 60 212; % 675
       ];

YLoadPhases = [0 0 0;   % 650
       1 1 1; %632
       0 0 0;  % 633
       1 1 1;  % 634
       0 1 0; % 645
       0 0 0; % 646
       1 1 1; % 671
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 1; % 611
       1 0 0; % 652
       %0 0 0;  % 692
       1 1 1; % 675
       ];


% reactive injections from capacitors, kVAr
Q_Y_cap = [0 0 0;   % 650
       0 0 0; %632
       0 0 0;  % 633
       0 0 0;  % 634
       0 0 0; % 645
       0 0 0; % 646
       0 0 0; % 671
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 100; % 611
       0 0 0; % 652
       %0 0 0;  % 692
       200 200 200; % 675
       ];

CapPhases = [0 0 0;   % 650
       0 0 0; %632
       0 0 0;  % 633
       0 0 0;  % 634
       0 0 0; % 645
       0 0 0; % 646
       0 0 0; % 671
       %0 0 0;  % 680
       0 0 0;% 684
       0 0 1; % 611
       0 0 0; % 652
       %0 0 0;  % 692
       1 1 1; % 675
       ];

q_max_cap = sum(sum(Q_Y_cap)); % total (reference) reactive injection from capacitors

% substation power tracking
sub_p_ref = 2800;            % substation reference real power to track, kW (80.78% * overline_sub_p)
sub_q_ref = 1000;           % substation reference reactive power to track, kVAR (80.88% * overline_sub_q - q_injection_all_capacitors)
% substation total maximum power (as given by the original data)
overline_sub_p = sum(sum(P_Y_l)) + sum(sum(P_Delta_l));
overline_sub_q = sum(sum(Q_Y_l)) + sum(sum(Q_Delta_l));

% Line configurations----------------------------------------------------
% How to convert to DSS format for R X: take R(Z601(1,1)) = 0.3465 as example. let
% R_DSS be the value in ~rmatrix. Then 0.3465*1850/5280 = R_DSS*1.85 where
% 1850 is the line_length in this code (also in IEEE doc) and 1.85 is the
% length value in the New Line object in .dss file
% How to convert to DSS format for Y: take Y601(1,1) = 6.2998e-6
% siemens/mile as example.
% 6.2998e-6*1230ft/(5280ft/mile) = C_DSS e-9 * 2pi60 * 1.230 kft where
% 1.230 is the value that appears in the New Line Lenth in .dss
Z601 = [0.3465+ 1.0179i 0.1560+0.5017i 0.1580+0.4236i;
       0.1560+0.5017i  0.3375+1.0478i  0.1535+0.3849i  ;
       0.1580+0.4236i  0.1535+0.3849i  0.3414+1.0348i]/Zbase;
Y601 = jay*[6.2998 -1.9958 -1.2595;
    -1.9958 5.9597 -0.7417;
    -1.2595 -0.7417 5.6386]*(10^-6)/Ybase;

Z602 = [0.7526 + 1.1814i  0.1580 + 0.4236i 0.1560 + 0.5017i;
        0.1580 + 0.4236i  0.7475 + 1.1983i 0.1535 + 0.3849i;
        0.1560 + 0.5017i  0.1535 + 0.3849i 0.7436 + 1.2112i]/Zbase;
Y602 = jay*[5.6990  -1.0817  -1.6905;
            -1.0817 5.1795 -0.6588;
            -1.6905 -0.6588 5.4246]*(10^-6)/Ybase;

Z603 = [0.0 + 0.0i  0.0 + 0.0i 0.0 + 0.0i;
        0.0 + 0.0i  1.3294 + 1.3471i 0.2066 + 0.4591i;
        0.0 + 0.0i  0.2066 + 0.4591i 1.3238 + 1.3569i]/Zbase;
Y603 = jay*[0.0 0.0 0.0;
       0.0 4.7097  -0.8999;
       0.0 -0.8999 4.6658]*(10^-6)/Ybase;

Z604 = [1.3238 + 1.3569i  0.0 + 0.0i 0.2066 + 0.4591i;
        0.0 + 0.0i       0.0 + 0.0i   0.0 + 0.0i;
        0.2066 + 0.4591i  0.0 + 0.0i 1.3294 + 1.3471i]/Zbase;
Y604 = jay*[4.6658 0.0 -0.8999;
             0.0  0.0  -0.0;
            -0.8999 0.0 4.7097]*(10^-6)/Ybase;


Z605 = [0.0 + 0.0i      0.0 + 0.0i   0.0 + 0.0i;
        0.0 + 0.0i       0.0 + 0.0i   0.0 + 0.0i;
        0.0 + 0.0i        0.0 + 0.0i  1.3292 + 1.3475i]/Zbase;
Y605 = jay*[0.0  0.0  0.0;
            0.0  0.0  0.0;
            0.0  0.0  4.5193]*(10^-6)/Ybase;

Z606 = [0.7982 + 0.4463i      0.3192 + 0.0328i   0.2849 - 0.0143i;
        0.3192 + 0.0328i       0.7891 + 0.4041i   0.3192 + 0.0328i;
        0.2849 - 0.0143i        0.3192 + 0.0328i  0.7982 + 0.4463i]/Zbase;
Y606 = jay*[96.8897  0.0  0.0;
            0.0   96.8897  0.0;
            0.0  0.0  96.8897]*(10^-6)/Ybase;


Z607 = [1.3425 + 0.5124i      0.0 + 0.0i   0.0 + 0.0i;
        0.0 + 0.0i       0.0 + 0.0i   0.0 + 0.0i;
        0.0 + 0.0i        0.0 + 0.0i  0.0 + 0.0i]/Zbase;
Y607 = jay*[88.9912  0.0  0.0;
            0.0  0.0  0.0;
            0.0  0.0  0.0]*(10^-6)/Ybase;


% % calculate transformer XFM-1 equivalent series impedance
% % using 633, 634 voltages (in V), 633--634 currents (in A), and loss
% % in (W)
% % solve an optimization to get the best estimate of Z_XFM-1, in Ohm
% V633 = [1.018*exp(jay*(-2.56)/180*pi);
%       1.040*exp(jay*(-121.77)/180*pi);
%       1.015*exp(jay*117.82/180*pi)]*4160;   % in Volts
% V634 = [0.994*exp(jay*(-3.23)/180*pi);
%       1.022*exp(jay*(-122.22)/180*pi);
%       0.996*exp(jay*117.34/180*pi)]*4160;   % in Volts
% I633634 = [81.33*exp(jay*(-37.74)/180*pi);
%       61.12*exp(jay*(-159.09)/180*pi);
%       62.71*exp(jay*80.47/180*pi)];            % In Ampere
% loss633634 = [2513;
%               1420;
%               1494];               % in W
%
% % call cvx to solve for Z_XFM
% cvx_begin
%
% cvx_solver sedumi
%
% variable Z_XFM(3,3) complex symmetric
%
% minimize(  1e-2 * sum(  pow_abs( real(diag( Z_XFM * (I633634*I633634') )) - loss633634 , 2) ) +  sum( pow_abs( V633 - Z_XFM*I633634 - V634 , 2) ) );
% % minimize(    sum( pow_abs( V633 - Z_XFM*I633634 - V634 , 2) ) );
%
% cvx_end
% The result of the optimization above is the following:
Z_XFM = [0.3842 + 1.1952i  0.0 + 0.0i  0.0 + 0.0i;
         0.0 + 0.0i    0.3870 + 1.1785i  0.0+0.0i;
         0.0 + 0.0i    0.0 + 0.0i     0.3871 + 1.2061i]/Zbase;   % this is the ohm value of Xformer impedance. Consider Z_XFM as ohm/mile value, then the length of the equivalent line should be 1 mile = 5280 ft
Y_XFM = zeros(3,3)/Ybase;  % ignore shunt admittance of transformer

% Switch modeled as zeros impedance line
Z_switch = zeros(3,3)/Zbase;
Y_switch = zeros(3,3)/Ybase;

% arrange all configurations in a row
Zconfig = [Z601 Z602 Z603 Z604 Z605 Z606 Z607 Z_XFM Z_switch];
Yconfig = [Y601 Y602 Y603 Y604 Y605 Y606 Y607 Y_XFM Y_switch];
%------------------------------------------------------------------------


% Convert the configurations above to DSS format
Zconfig_dss = Zconfig*1000/5280/1.000;
Yconfig_dss = Yconfig*1000/5280/1.000/(2*pi*60)/(1e-9);

%------------------------------------------------------------------------
% from and to bus index
line_from = [2 2 3 5 1 8 2 7 8 7];
line_to =   [5 3 4 6 2 10 7 8 9 11];
% line configurations
% line_config = [3 2 8 3 1 7 1 4 5 6];
line_config = [3 2 8 3 1 6 1 4 5 6];
% line length
line_length = [500 500 5280 300 2000 800 2000 300 300 500];
% line phases
LinePhases = [0 1 1;
    1 1 1;
    1 1 1;
    0 1 1;
    1 1 1;
    1 0 0;
    1 1 1;
    1 0 1;
    %1 1 1;
    %1 1 1;
    0 0 1;
    1 1 1;
    ];
% mile = 5280 feet
convfm = (1/5280);


% line series impedance matrices
Zseries = zeros(3,3*Nline);   % all series impedance matrices are arranged in a row
for ll = 1:Nline
    Zseries(:,3*(ll-1)+1:3*ll) = Zconfig(:, 3*(line_config(ll)-1)+1 : 3*line_config(ll))*line_length(ll)*convfm;    % per unit based on Zbase, already scaled in Zconfig
end


% Line shunt admittance aggregated at terminal buses (each side 0.5)
% Self shunt if there is any
Yshunt = zeros(3,3*Nnode);
for ll = 1:Nline
    terminal_shunt_thisline = 0.5 * Yconfig(:, 3*(line_config(ll)-1)+1 : 3*line_config(ll))*line_length(ll)*convfm;    % per unit based on Ybase, already scaled in Yconfig
    Yshunt(:,3*(line_from(ll)-1)+1 : 3*line_from(ll)) = Yshunt(:,3*(line_from(ll)-1)+1 : 3*line_from(ll)) + terminal_shunt_thisline;
    Yshunt(:,3*(line_to(ll)-1)+1 : 3*line_to(ll)) = Yshunt(:,3*(line_to(ll)-1)+1 : 3*line_to(ll)) + terminal_shunt_thisline;
end


Ynet = zeros(Nnode*3);
for ll = 1:Nline
    ff = line_from(ll);
    tt = line_to(ll);
    Ynet(3*ff-2:3*ff,3*tt-2:3*tt) = -pinv(Zseries(:,3*ll-2:3*ll));
    Ynet(3*tt-2:3*tt,3*ff-2:3*ff) = -pinv(Zseries(:,3*ll-2:3*ll));
    Ynet(3*ff-2:3*ff,3*ff-2:3*ff) = Ynet(3*ff-2:3*ff,3*ff-2:3*ff) + pinv(Zseries(:,3*ll-2:3*ll)) + Yshunt(:,3*ll-2:3*ll);
    Ynet(3*tt-2:3*tt,3*tt-2:3*tt) = Ynet(3*tt-2:3*tt,3*tt-2:3*tt) + pinv(Zseries(:,3*ll-2:3*ll)) + Yshunt(:,3*ll-2:3*ll);
end