function [ thickness, neff, width ] = getThickAndNAndWidFromLayers( layers )

thickness = zeros(1,length(layers));
neff = zeros(1,length(layers));
width = zeros(1,length(layers));

for layer = 1:length(layers)
    thickness(1,layer) = layers(layer,1)*1e-6;
    neff(1,layer) = layers(layer,2);
    width(1,layer) = layers(layer,3)*1e-6;
end


end