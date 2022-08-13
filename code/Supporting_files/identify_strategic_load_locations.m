buses_1_hop = [];
buses_2_hop = [];
buses_3_hop = [];
buses_4_hop = [];
load_buses_near_high_MW_gens = [];
load_buses_near_low_MW_gens  = [];
load_buses_near_high_inertia_gens = [];
load_buses_near_low_inertia_gens  = [];

for i_gen = 1:length(generator_locations)
    % Gives all 1 hop buses only
    buses_1_hop = [buses_1_hop; nearest(G,generator_locations(i_gen),1)];
    % Gives all 2 hops buses only
    buses_2_hop = [buses_2_hop; setdiff(nearest(G,generator_locations(i_gen),2),nearest(G,generator_locations(i_gen),1))];
    % Gives all 3 hops buses only
    buses_3_hop = [buses_3_hop; setdiff(nearest(G,generator_locations(i_gen),3),nearest(G,generator_locations(i_gen),2))];
    % Gives all 4 hops buses only
    buses_4_hop = [buses_4_hop; setdiff(nearest(G,generator_locations(i_gen),4),nearest(G,generator_locations(i_gen),3))];
end

buses_1_hop=unique(buses_1_hop);
buses_2_hop=unique(buses_2_hop);
buses_3_hop=unique(buses_3_hop);
buses_4_hop=unique(buses_4_hop);

buses_1_hop = intersect(buses_1_hop, load_locations);
buses_2_hop = intersect(setdiff(buses_2_hop,buses_1_hop), load_locations);
buses_3_hop = intersect(setdiff(buses_3_hop,[buses_2_hop;buses_1_hop]), load_locations);
buses_4_hop = intersect(setdiff(buses_4_hop,[buses_3_hop;buses_2_hop;buses_1_hop]), load_locations);

for i_gen = 1:5
    load_buses_near_high_MW_gens = [load_buses_near_high_MW_gens; nearest(G,generator_size_ordering(end-i_gen),2)];
    load_buses_near_low_MW_gens  = [load_buses_near_low_MW_gens; nearest(G,generator_size_ordering(i_gen),2)];
    
    load_buses_near_high_inertia_gens = [load_buses_near_high_inertia_gens; nearest(G,inertia_ordering(end-i_gen),2)];
    load_buses_near_low_inertia_gens  = [load_buses_near_low_inertia_gens; nearest(G,inertia_ordering(i_gen),2)];
end

load_buses_near_high_MW_gens = intersect(load_buses_near_high_MW_gens(1:min(15,length(load_buses_near_high_MW_gens))), load_locations);
load_buses_near_low_MW_gens  = intersect(load_buses_near_low_MW_gens, load_locations);

load_buses_near_high_inertia_gens = intersect(load_buses_near_high_inertia_gens(1:min(15,length(load_buses_near_high_inertia_gens))), load_locations);
load_buses_near_low_inertia_gens  = intersect(load_buses_near_low_inertia_gens(1:min(25,length(load_buses_near_low_inertia_gens))), load_locations);
% -------------------------------------------------------------------------