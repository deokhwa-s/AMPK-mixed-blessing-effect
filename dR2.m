function y = dR2(t, x, param, varargin)
%-----------------------------------------------------------------------
%Function for the system of ODEs for the mTORC1/AMPK/ULK1/DEPTOR Pathway
%-----------------------------------------------------------------------
%Input arguments:
%   t (float)    : Time (Not directly used within the function but set to solve the ODE)
%   x (float)    : 20 by 1 column vector containing concentrations of network species   
%   param (float): 58 by 1 column vector containing parameters of the system
%   
%   An extra parameter can be passed on to set total SIRT1 manually
%       e.g. dR2(t, x, param) -> SIRT1_total = 375 (Default)
%            dR2(t, x, param, 100) -> SIRT1_total = 100
%-----------------------------------------------------------------------
%Output:
%   y (float)    : 20 by 1 column vector containing first order derivatives of the output signal
%                  (concentration of each component)
%-----------------------------------------------------------------------

    %Define each network component to the corresponding vector element
    IR            = x(1);
    pIR           = x(2);
    IRS           = x(3);
    pIRS          = x(4);
    iIRS          = x(5);
    AKT           = x(6);
    pAKT          = x(7);
    mTORC1        = x(8);
    pmTORC1       = x(9);
    mTORC2        = x(10);
    pmTORC2       = x(11);
    mTORC1_DEPTOR = x(12);
    mTORC2_DEPTOR = x(13);
    DEPTOR        = x(14);
    pDEPTOR       = x(15);
    AMPK          = x(16);
    pAMPK         = x(17);
    SIRT1         = x(18);
    ULK1          = x(19);
    pULK1         = x(20);

    % Set total SIRT1 depending on number of input parameters
    % SIRT1_total is assumed to be constant
    if nargin == 3
        SIRT1_total = 375; % Default total SIRT1
    else
        SIRT1_total = varargin{1}; % Set total SIRT1 to user input
    end

    %Set the parameters for the rate equations
    V_IR                        = param(1);         % Rate of activation of IR
    Km_IR                       = param(2);         % MM constant for the activation of IR
    V_pIR                       = param(3);         % Rate of deactivation of IR
    Km_pIR                      = param(4);         % MM constant for the deactivation of IR
    K_IRS_by_pIR                = param(5);         % Rate of activation of IRS via pIR
    Km_IRS_by_pIR               = param(6);         % MM constant for the activation of IRS via pIR
    V_pIRS                      = param(7);         % Rate of deactivation of IRS
    Km_pIRS                     = param(8);         % MM constant for the deactivation of IRS
    K_AKT_by_pIRS               = param(9);         % Rate of activation of AKT via pIRS
    Km_AKT_by_pIRS              = param(10);        % MM constant for the activation of AKT via pIRS
    K_AKT_by_pmTORC2            = param(11);        % Rate of activation of AKT via pmTORC2
    Km_AKT_by_pmTORC2           = param(12);        % MM constant for the activation of AKT via pmTORC2
    V_pAKT                      = param(13);        % Rate of deactivation of AKT
    Km_pAKT                     = param(14);        % MM constant for the deactivation of AKT
    K_mTORC1_by_pAKT            = param(15);        % Rate of activation of mTORC1 via pAKT
    Km_mTORC1_by_pAKT           = param(16);        % MM constant for the activation of mTORC1 via pAKT
    K_pmTORC1                   = param(17);        % Rate of background deactivation of mTORC1
    K_pmTORC1_by_pAMPK          = param(18);        % Rate of deactivation of mTORC1 via pAMPK
    Km_pmTORC1_by_pAMPK         = param(19);        % MM constant for the deactivation of mTORC1
    K_pmTORC1_by_pULK1          = param(20);        % Rate of deactivation of mTORC1 via pULK1
    Km_pmTORC1_by_pULK1         = param(21);        % MM constant for the deactivation of mTORC1 via pULK1
    K_mTORC2_by_pIRS            = param(22);        % Rate of activation of mTORC2 via pIRS
    Km_mTORC2_by_pIRS           = param(23);        % MM constant for the activation of mTORC2 via pIRS
    K_mTORC2_by_pAMPK           = param(24);        % Rate of activation of mTORC2 via pAMPK
    Km_mTORC2_by_pAMPK          = param(25);        % MM constant for the activation of mTORC2 via pAMPK
    V_pmTORC2                   = param(26);        % Rate of deactivation of mTORC2
    Km_pmTORC2                  = param(27);        % MM constant for the deactivation of mTORC2
    K_DEPTOR_by_pmTORC1         = param(28);        % Rate of activation of DEPTOR via pmTORC1
    Km_DEPTOR_by_pmTORC1        = param(29);        % MM constant for the activation of DEPTOR via pmTORC1
    K_DEPTOR_by_pmTORC2         = param(30);        % Rate of activation of DEPTOR via pmTORC2
    Km_DEPTOR_by_pmTORC2        = param(31);        % MM constant for the activation of DEPTOR via pmTORC2
    V_pDEPTOR                   = param(32);        % Rate of deactivation of DEPTOR
    Km_pDEPTOR                  = param(33);        % MM constant for the deactivation of DEPTOR
    K_mTORC1_DEPTOR_form        = param(34);        % Rate of formation of the mTORC1-DEPTOR complex
    K_mTORC1_DEPTOR_diss        = param(35);        % Rate of dissociation of the mTORC1-DEPTOR complex
    K_mTORC2_DEPTOR_form        = param(36);        % Rate of formation of the mTORC2-DEPTOR complex
    K_mTORC2_DEPTOR_diss        = param(37);        % Rate of dissociation of the mTORC2-DEPTOR complex
    K_IRS_to_iIRS               = param(38);        % Rate of inactivation of IRS
    Km_IRS_to_iIRS              = param(39);        % MM constant for the inactivation of IRS
    V_iIRS                      = param(40);        % Rate of activation of IRS from iIRS
    Km_iIRS                     = param(41);        % MM constant for the activation of IRS from iIRS
    K_AMPK                      = param(42);        % Rate of background activation of AMPK
    K_AMPK_by_SIRT1             = param(43);        % Rate of activation of AMPK via SIRT1
    Km_AMPK                     = param(44);        % MM constant for the activation of AMPK
    K_pAMPK                     = param(45);        % Rate of background deactivation of AMPK
    K_pAMPK_by_pULK1            = param(46);        % Rate of deactivation of AMPK via pULK1
    K_pAMPK_by_pmTORC1          = param(47);        % Rate of deactivation of AMPK via pmTORC1
    Km_pAMPK                    = param(48);        % MM constant for the deactivation of AMPK
    K_SIRT1                     = param(49);        % Rate of background activation of SIRT1
    K_SIRT1_by_pAMPK            = param(50);        % Rate of activation of SIRT1 via pAMPK
    Km_SIRT1                    = param(51);        % MM constant for the activation of SIRT1
    K_SIRT1_diss                = param(52);        % Rate of dissociation of SIRT1
    K_ULK1                      = param(53);        % Rate of background activation of ULK1
    K_ULK1_by_pAMPK             = param(54);        % Rate of activation of ULK1 via pAMPK
    Km_ULK1                     = param(55);        % MM constant for the activation of ULK1
    K_pULK1                     = param(56);        % Rate of background deactivation of pULK1
    K_pULK1_by_pmTORC1          = param(57);        % Rate of deactivation of pULK1 via pmTORC1
    Km_pULK1                    = param(58);        % MM constant for the deactivation of ULK1







    %Set the rate equation for each reaction in the system
    r1  = (V_IR*IR)/(Km_IR+IR); % IR...pIR
    r2  = (V_pIR*pIR)/(Km_pIR+pIR); % pIR...IR
    r3  = (K_IRS_by_pIR*pIR*IRS)/(Km_IRS_by_pIR+IRS); % IRS...pIRS
    r4  = (V_pIRS*pIRS)/(Km_pIRS+pIRS); % pIRS...IRS
    r5  = (K_AKT_by_pIRS*pIRS*AKT)/(Km_AKT_by_pIRS+AKT)+(K_AKT_by_pmTORC2*pmTORC2*AKT)/(Km_AKT_by_pmTORC2+AKT); % AKT...pAKT
    r6  = (V_pAKT*pAKT)/(Km_pAKT+pAKT); % pAKT...AKT
    r7  = (K_mTORC1_by_pAKT*pAKT*mTORC1)/(Km_mTORC1_by_pAKT+mTORC1); % mTORC1...pmTORC1
    r8  = (K_pmTORC1 + K_pmTORC1_by_pAMPK*pAMPK)*pmTORC1/(Km_pmTORC1_by_pAMPK+pmTORC1)+(K_pmTORC1_by_pULK1*pULK1)*pmTORC1/(Km_pmTORC1_by_pULK1+pmTORC1); % pmTORC1....mTORC1
    r9  = (K_mTORC2_by_pIRS*pIRS*mTORC2)/(Km_mTORC2_by_pIRS+mTORC2)+(K_mTORC2_by_pAMPK*pAMPK*mTORC2)/(Km_mTORC2_by_pAMPK+mTORC2); % mTORC2...pmTORC2
    r10 = (V_pmTORC2*pmTORC2)/(Km_pmTORC2+pmTORC2); % pmTORC2...mTORC2
    r11 = (K_DEPTOR_by_pmTORC1*pmTORC1*DEPTOR)/(Km_DEPTOR_by_pmTORC1+DEPTOR)+(K_DEPTOR_by_pmTORC2*pmTORC2*DEPTOR)/(Km_DEPTOR_by_pmTORC2+DEPTOR); % DEPTOR...pDEPTOR
    r12 = (V_pDEPTOR*pDEPTOR)/(Km_pDEPTOR+pDEPTOR); % pDEPTOR...DEPTOR
    r13 = (K_mTORC1_DEPTOR_form*mTORC1*DEPTOR)-(K_mTORC1_DEPTOR_diss*mTORC1_DEPTOR); % mTORC1+DEPTOR...mTORC1_DEPTOR
    r14 = (K_mTORC2_DEPTOR_form*mTORC2*DEPTOR)-(K_mTORC2_DEPTOR_diss*mTORC2_DEPTOR); % mTORC2+DEPTOR...mTORC2_DEPTOR 
    r15 = (K_IRS_to_iIRS*pmTORC1*IRS)/(Km_IRS_to_iIRS+IRS); % IRS...iIRS
    r16 = (V_iIRS*iIRS)/(Km_iIRS+iIRS); % iIRS...IRS
    r17 = (K_AMPK+K_AMPK_by_SIRT1 * SIRT1)*AMPK/(Km_AMPK + AMPK); % AMPK...pAMPK
    r18 = (K_pAMPK+K_pAMPK_by_pULK1 * pULK1 + K_pAMPK_by_pmTORC1 * pmTORC1)*pAMPK/(Km_pAMPK + pAMPK); % pAMPK...AMPK
    r19 = (K_SIRT1+K_SIRT1_by_pAMPK * pAMPK)*(SIRT1_total - SIRT1) /(Km_SIRT1 + SIRT1_total - SIRT1) - K_SIRT1_diss * SIRT1; % pAMPK...SIRT
    r20 = (K_ULK1+K_ULK1_by_pAMPK * pAMPK)*ULK1/(Km_ULK1 + ULK1); % ULK1...pULK1
    r21 = (K_pULK1+K_pULK1_by_pmTORC1 * pmTORC1)*pULK1/(Km_pULK1 + pULK1); % pULK1...ULK1


    %Set the ODEs using the rate equations defined above. 
    dIR            = r2-r1;
    dpIR           = r1-r2;
    dIRS           = r4+r16-r3-r15;
    dpIRS          = r3-r4;
    diIRS          = r15-r16;
    dAKT           = r6-r5;
    dPAKT          = r5-r6;
    dmTORC1        = r8-r7-r13;
    dpmTORC1       = r7-r8;
    dmTORC2        = r10-r9-r14;
    dpmTORC2       = r9-r10;
    dmTORC1_DEPTOR = r13;
    dmTORC2_DEPTOR = r14;
    dDEPTOR        = r12-r11-r13-r14;
    dpDEPTOR       = r11-r12;
    dAMPK          = r18-r17; 
    dpAMPK         = r17-r18; 
    dSIRT1         = r19; 
    dULK1          = r21-r20; 
    dpULK1         = r20-r21; 


    %Output the system of ODEs as a column vector 
    y=[
    dIR;
    dpIR;
    dIRS;
    dpIRS;
    diIRS;
    dAKT;
    dPAKT;
    dmTORC1;
    dpmTORC1;
    dmTORC2;
    dpmTORC2;
    dmTORC1_DEPTOR;
    dmTORC2_DEPTOR;
    dDEPTOR;
    dpDEPTOR;
    dAMPK; 
    dpAMPK; 
    dSIRT1; 
    dULK1; 
    dpULK1; 
    ];
end









