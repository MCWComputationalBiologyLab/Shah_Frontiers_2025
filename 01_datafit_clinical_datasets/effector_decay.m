function deff = effector_decay(t, effector, parameters)

kd = parameters(1);

% system of ODEs (as cells)

deff = -1*kd*effector; % cell death

end