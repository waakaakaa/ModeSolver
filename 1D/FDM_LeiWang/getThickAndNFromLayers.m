function [ thickness, neff ] = getThickAndNFromLayers( layers )

thickness = zeros(1,length(layers));
neff = zeros(1,length(layers));

for layer = 1:length(layers)
    thickness(1,layer) = layers(layer,1)*1e-6;
    neff(1,layer) = layers(layer,2);
end


end