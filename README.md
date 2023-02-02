# bsb-CMR

This is data and code which contributed to the following paper:

SPORADIC, WARM-WATER RECRUITMENT PULSES AND AGGREGATION-BASED FISHING DRIVE MULTIDECADAL BOOM AND BUST CYCLES IN A TEMPERATE AGGREGATE SPAWNER 

Erica T. Jarvis Mason(1)*, Thomas V. Riecke(2), Lyall F. Bellquist(1,3), Daniel J. Pondella II(4), Brice X. Semmens(1)

(1)Scripps Institution of Oceanography, University of California San Diego, California, USA
(2)Department of Ecosystem and Conservation Sciences, University of Montana, USA
(3)California Oceans Program, The Nature Conservancy, California, USA
(4)Vantuna Research Group, Occidental College, California, USA

R-script files:
1) Double-tag_TagRetentionModel_Masonetal.R (Tag retention model for double-tagging data (2 identical tags) with individual recovery times) Includes the following sections:
Explore relationship between the ratio of fish with 2 vs 1 tag r --------
Simulate data and define model ------------------------------------------
Test model --------------------------------------------------------------
Load and format Kelp Bass double tag data -------------------------------
Tag retention model (tag loss as a function of age of tag) --------------
Run tag retention model on double-tag data ------------------------------
Plot cumulative and non-cumulative tag retention over time --------------

2) Age_Growth_Masonetal.R (Age and Growth Model Fits (Traditional Von Bertalanffy and Francis parameterization)) Includes the following sections:

3) CMR_wGrowth_Masonetal.R (Capture-mark-reencounter (CMR) data simulation and model testing) Includes the following sections:
Simulate CMR data--------------------------------------------------------
Define CMR model in JAGS-------------------------------------------------
Set up and run model in JAGS---------------------------------------------

4) 60s_Final_Masonetal.R (Capture-mark-reencounter (CMR) model for modeling the 1960s CDFW CMR data)

5) 90s_Final_Masonetal.R (Capture-mark-reencounter (CMR) model for modeling the 1990s CDFW CMR data)

6) 10s_Final_Masonetal.R (Capture-mark-reencounter (CMR) model for modeling the 2010s Bellquist CMR data)


