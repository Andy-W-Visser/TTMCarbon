# TTMCarbon
These MATLAB programes illustrate the use of the transport matrix OCIM in estimating the residence times of injected carbon onto the ocean interior. 
The matrix itself is has been developed by Tim deVries at the University of Sanata Barbra aand can be downloaded at
https://tdevries.eri.ucsb.edu/models-and-data-products/.
Any use of this product should be referenced to 
DeVries, T. (2014), The oceanic anthropogenic CO2 sink: Storage, air-sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.

Within the code, the TMM is imported as load('CTL.mat')

PassageTimes : calculates the first and last passage time from and to the sea floor.

PassageTimes1000 : calculates the first and last passage time from and to a depth of 1000 m.

krillinterp001 : calculates the sequestered carbon and residence time associated with krill poop in the Southern Ocean. Data for distribution is inthe .csv file. Model used in Cavan EL, Mackay N, Hill SL, Atkinson A, Belcher A, Visser AW. Antarctic krill sequester similar amounts of carbon to key coastal blue carbon habitats. Nature Communications. 2024 Sep 8;15(1):7842.

TransportMatrixAppSeaQuester_exported : the code to generate a MatLab application
