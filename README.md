# Ecological recoveries in at-risk metapopulations
This project folder contains the code and analyses for a simulation model on understanding the recovery dynamics in heavily disturbed metapopulations governed by local density dependence, dispersal, spatial habitat networks, and disturbance regimes.

## Walking through an example

For illustration, we show two example recovery regimes for metapopulations that varied only in network topology when disturbance was localized and uneven. We first examine a metapopulation with linear topology composed of 16 identical patches, a dispersal rate of 1%, and low stochasticity (0.1%). This lower dispersal limits the strength of the “rescue regime” from dominating patterns in recovery. As shown in the below figure (upper panels), the linear network has a slowed and modular recovery pattern where recovery cascades outwards from the last remaining patch (e.g., California sea otters – see Case Studies in main paper). In contrast, the dendritic network recovered sub-regions within the whole network as dispersal more quickly spread rescue effects to neighboring patches (below figure – lower panels).

![example model](https://github.com/klwilson23/SpatialWatershed/blob/master/Figures/example%20landscape%20results.jpeg?raw=true)

## Highlights

MA few common spatial recovery regimes emerged from the interplay between local productivity, dispersal, and spatial disturbances (Table 1; Figure S16). Overall, clustering analyses of our model results found evidence for five common outcomes (Figure S16): (1) resilient recovery – fast metapopulation recoveries, (2) slow recovery – metapopulation recoveries were slow but no other metrics were affected, (3) hidden collapses – metapopulation recovered (or were recovering) but 60% of patches remained unoccupied, (4) lost capacity – slow recovery rates, patch occupancy was ≤60%, long-term relative production was ≤80% pre-disturbance, and the risk of not recovering was relatively high, and (5) critical risk – metapopulation abundances remained ≤10% of pre-disturbance. Resilient recoveries were most common when disturbance was uniformly applied across the network and local patches had homogenous demographic parameters (Table S2, Figures S14 and S17). Local disturbances generated poorer recovery patterns (defined as any non-resilient recovery), and the frequency of these outcomes were exacerbated when patches varied in their per-capita productivity. Locally even disturbances tended to reduce recovery rates by 16% compared to uniform disturbances. The poorest recovery outcomes emerged when disturbance was locally uneven. In these scenarios, recovery was slowed by 42% (from 0.82∙yr-1 to 0.47∙yr-1 on average), relative production was often reduced by 28%, long-term patch occupancy was reduced by 59%, and the risk of non-recovery was increased by almost seven-fold (from 2.3% to 18%). These dynamics illustrate some of the context dependence underlying how local and regional processes interact to shape metapopulation recoveries. 

![Some common outcomes can be shown here](https://github.com/klwilson23/SpatialWatershed/blob/master/Figures/surprising%20outcomes.jpeg?raw=true)
