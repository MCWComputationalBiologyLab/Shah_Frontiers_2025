function dtumor = logistic_growth(t, tumor, parameters)

a = parameters(1);
b = parameters(2);

% system of ODEs (as cells)

dtumor = a*tumor*(1-(b*tumor)); % logistic tumor growth

end