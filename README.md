# bsb-CMR

This is data and code that contributed to the following paper:

SPORADIC, WARM-WATER RECRUITMENT PULSES AND AGGREGATION-BASED FISHING DRIVE MULTIDECADAL BOOM AND BUST CYCLES IN A TEMPERATE AGGREGATE SPAWNER (tentative title) 

*Erica T. Jarvis Mason<sup>1</sup>, Thomas V. Riecke<sup>2</sup>, Lyall F. Bellquist<sup>1,3</sup>, Daniel J. Pondella II<sup>4</sup>, Brice X. Semmens<sup>1</sup>*

<sup>1</sup>Scripps Institution of Oceanography, University of California San Diego, California, USA  
<sup>2</sup>Department of Ecosystem and Conservation Sciences, University of Montana, USA  
<sup>3</sup>California Oceans Program, The Nature Conservancy, California, USA  
<sup>4</sup>Vantuna Research Group, Occidental College, California, USA

**R-script files:**

1. Double-tag_TagRetentionModel_Masonetal.R (Tag retention model for double-tagging data with individual recovery times) Includes the following sections:
   - Explore relationship between the ratio of fish with 2 vs 1 tag
   - Simulate data and define model
   - Test model
   - Load and format Kelp Bass double tag data
   - Tag retention model (tag loss as a function of age of tag)
   - Run tag retention model on double-tag data
   - Plot cumulative and non-cumulative tag retention over time

2. Age_Growth_Masonetal.R (Age and Growth Model Fits using the Traditional Von Bertalanffy and Francis parameterizations)

3. CMR_wGrowth_Masonetal.R (Capture-mark-reencounter data simulation and model testing) Includes the following sections:
   - Simulate CMR data
   - Define CMR model in JAGS
   - Set up and run model in JAGS

4. 60s_Final_Masonetal.R (Capture-mark-reencounter model for modeling the 1960s CDFW CMR data)

5. 90s_Final_Masonetal.R (Capture-mark-reencounter model for modeling the 1990s CDFW CMR data)

6. 10s_Final_Masonetal.R (Capture-mark-reencounter model for modeling the 2010s L. Bellquist CMR data)


