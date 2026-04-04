function dCells = cyto(~, C, parameters, parameters_logfit)

  %tumor parameters
  a = parameters_logfit(1);
  b = parameters_logfit(2);

  %fitting parameters
  Kappa = parameters(1);
  k = parameters(2);
  kc = parameters(3);
  kp2 = parameters(4);
  kp3 = parameters(5);
  l = parameters(6);
  s = parameters(7);

  %cell populations
  T = C(1);
  CART = C(2);
  TC = C(3);

  %cytotoxicity control
  
  if T < 1
    T = 0;
    FC = 0; 
  elseif CART < 1
    CART = 0;
    FC = 0;
  else
    FC = kc * (((CART/T)^l) / ((s^l) + ((CART/T)^l))) * T;
  end

  %CART control

  if CART < 1
    CART = 0;
    FP = 0;
  else
    FP = kp2 * (FC^2)/((k^2)+(FC^2)) * CART;
  end

  %TC control

  if TC < 1
    TC = 0;
  end

  %total effectors
  
  CP = CART + TC;

  if CP < 1
      CP = 1;
  end

  %system of odes (as cells)

   %tumor
   dTumor = a*T*(1-(b*T)) - FC;                           % Tumor cells

   %Em1
   dCART = FP*log(Kappa/CP) + -1*kp3*CART;                % CART cells
   dTC = -1*kp3*TC;                                       % T cells

   dCells = [dTumor dCART dTC (FC-C(4)) (FP-C(5))]';       % Output
  
end