require(LocalControlStrategy)

#input the radon data of Krstic & Obenchain (2016).
data(radon)
# outcome:    lcanmort ...continuous
# treatment:  lnradon  ...continuous Exposure (natural log, Winsorized to exceed -3.0)

# Define Cluster Hierarchy for UNSUPERVISED, nonparametric analyses...
xvars  <- c("obesity","over65","cursmoke")
hclobj <- LCcluster(radon, xvars)   # defalut clustering method = "ward.D"
hclobj  
plot(hclobj)

# Save Local Control basic parameter settings to an environment that will be updated...
e <- LCsetup(hclobj, radon, lnradon, lcanmort)
ls.str(e)

# Compute and Save LRC distributions for a range of K = Number of Clusters...
mort010 <- lrcagg( 10, e)
mort050 <- lrcagg( 50, e)   # Average Cluster Size: ~58 US Counties
plot(mort050, show="ecdf", e)
mort100 <- lrcagg(100, e)
mort200 <- lrcagg(200, e)   # Average Cluster Size: Only ~14 US Counties
plot(mort200, show="ecdf", e)

# "Sensitivity Analysis" Summary...
LCcompare(e) 
# LTD Distribution for 50 Clusters appears to Optimize Variance-Bias Trade-Offs... 

# Compute and Save IV distributions for 2 values of K = Numbers of Clusters...
iv050  <- ivadj(mort050)
plot(iv050)  # graphical display

iv200  <- ivadj(mort200)
plot(iv200)  # graphical display

# Confirm: Does the Observed LRC Distribution for 50 clusters really differ
#          from its Random Permutation counterpart?
system.time( conf050 <- confirm(mort050) )   # Simulation takes ~6 seconds. 
conf050
plot(conf050)
# Due to the presence of MANY ties in LRC distributions, the above p.value
# from ks.test() is badly BIASED downwards!
warnings()

# Simulate crude pmax.value for Kolmogorov-Smirnov D-statistic... 
system.time( ksd050 <- KSperm(conf050) )    # Simulation takes ~12 seconds.   
ksd050
plot(ksd050)
warnings()
