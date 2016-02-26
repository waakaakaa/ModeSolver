function [ layers ] = getLayers( t )
layers = [[3.0,    3.4, 8];...  %% substrate
    [t,    3.44, 8];...  %% slab
    [1-t,  3.44, 3];...     %% ridge
    [0.5,  1.00, 8]];... %% air
end

