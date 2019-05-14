source("make networks.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("popDynFn.R")
source("Metapop function.R")

library(mvtnorm)
library(marima)
library(diagram)

tiff("linear example.tiff",res=800,units="in",compression="lzw",height=6,width=6)
metaPop(Npatches=16,
        networkType="linear",
        patchDistance=1,
        Nburnin=50,
        NyrsPost=100,
        omega=0.01,
        m=1,
        alpha=2,
        metaK=1600,
        cv=1e-2,
        DistScenario="uniform",
        magnitude_of_decline=0.9,
        lagTime=1,
        prodType="Beverton-Holt",
        rho.time=1e-5,
        rho.dist=1e5,
        compensationLag=25,
        dataWeighting=0.1,
        alphaVariable=FALSE,
        kVariable=FALSE)
dev.off()

tiff("dendritic example.tiff",res=800,units="in",compression="lzw",height=6,width=6)
metaPop(networkType = "dendritic")
dev.off()

tiff("star example.tiff",res=800,units="in",compression="lzw",height=6,width=6)
metaPop(networkType = "star")
dev.off()

tiff("complex example.tiff",res=800,units="in",compression="lzw",height=6,width=6)
metaPop(networkType = "complex")
dev.off()
