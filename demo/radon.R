require(LocalControlStrategy)

# Use the radon data of Krstic & Obenchain (2016).
data(radon)
# outcome:    lcanmort ...continuous
# treatment:  lnradon  ...continuous Exposure (natural log, Winsorized to exceed -3.0)

# Define Cluster Hierarchy for UNSUPERVISED, nonparametric analyses...
xvars  <- c("obesity","over65","cursmoke")
hclobj <- LCcluster(radon, xvars)   # defalut clustering method = "ward.D"
hclobj  
plot(hclobj)

# Save Local Control basic parameter settings to an environment that will be Updated...
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

# Save and plot IV distributions for 2 values of K = Number of Clusters...
iv050  <- ivadj(mort050)
plot(iv050)  # graphical display
iv200  <- ivadj(mort200)
plot(iv200)  # graphical display

# Confirm: Does the Observed LRC distribution for 50 clusters truly differ from
# the Random Permutation NULL distribution assuming x_Covariates are Ignorable?
system.time( conf050 <- confirm(mort050) )   # Simulation takes ~6 seconds. 
conf050
plot(conf050)

# Simulate maximum p-value for observed Kolmogorov-Smirnov D-statistic... 
system.time( ksd050 <- KSperm(conf050) )    # Simulation takes ~12 seconds.   
ksd050
plot(ksd050)
