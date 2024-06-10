# Heterogenous Membranes with Metabolite Consumption

## Run 1
Set the consumption to the wrong particle type, so there actually is no consumption here.
Possibly accidentally overwritten in faulty file start.

## Run2 
Consumption now activated between metabolites and `active_vertices`

## Run3
Maybe it is not such a good idea to have consumption of active vertices happen at the same range as the cutoff of their interactions.
Vary the consumption range to be less than interaction range.

# Run 4
Made a mistake, the group `inert_vertices` should have been consuming.
Altered code, and rerun different column heights and interaction probabilities.