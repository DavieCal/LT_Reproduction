---
output:
  word_document:
    fig_caption: yes
    pandoc_args: --smart
    reference_docx: format.docx
  html_document:
    fig_caption: yes
csl: canadian-journal-of-fisheries-and-aquatic-sciences.csl
---

Title: Lake trout (*Savelinus namaycush*) reproductive behaviour in a northern lake. SPatial temporal clustering

Authors: David Thomas Callaghan^1^*, Paul James Blanchfield^1,2^, Matt M. Guzzo^1^, and Peter A. Cott^3^

^1^Department of Biological Sciences, University of Manitoba, 50 Sifton Road, Winnipeg, MB R3T 2N2, Canada
^2^ Freshwater Institute, Fisheries and Oceans Canada, 501 University Crescent, Winnipeg, MB, R3T 2N6, Canada 
^3^ Cumulative Impact Monitoring Program - Environment and Natural Resources, Government of the Northwest Territories, Box 1320, Yellowknife, Northwest Territories X1A 2L9, Canada.

*Corresponding Author Email: callagdt@myumanitoba.ca

Dave Dave Dave

#Abstract

This thesis investigates lake trout (*Salvelinus namaycush*) reproduction northern lakes. Lake trout have a broad distribution across Canada's North, yet most studies that describe reproductive habitat are from the southern extent of their range. I first assessed whether lake trout spawning habitat, typically characterized as wave-swept shoals with clean cobble that face predominant wind directions, is similar for a northern lake. Specifically, I examined a dozen sites around Alexie Lake, Northwest Territories, to test if physical habitat and wind exposure were important determinants of spawning site use and embryonic survival. Spawning occurred in ~2 m water depth, on 3–15 cm diameter clean substrate found on the leading edge of shoals that ended in a rock crib rising abruptly in nearshore regions around the lake. Wind direction was predominantly from the west, although it was highly variable within and among spawning seasons. I found evidence of lake trout spawning at each site examined, but not limited to shoals facing a predominate wind direction. High variation in embryonic survival (2–83%) from incubation trays was observed among spawning sites, suggesting a large gradient in habitat quality exists within a given lake. Modelled wind exposure did not predict embryonic survival, nor did physical characteristics that may influence interstitial water flow on spawning shoals. 
  I also provide a detailed description of the lake trout mating system at the whole lake scale in a typical northern boreal lake. Using an acoustic telemetry monitoring system and a novel spatial temporal clustering analysis, I was able to quantify lake trout spawning movements and behaviours over the course of an entire spawning season. Lake trout were found to cluster on spawning shoals virtually around the entire nearshore region of Alexie Lake, as well as around several islands, which appears to further confirm previous findings that subtable spawning habitat is abundant in Alexie Lake. Males arrived earlier than females and spent longer durations on spawning shoals over the course of the spawning season. Males formed >4 times as many spawing clusters and visited more sites than females. Spawning clusters predominantly were formed at night but were also observed during daylight hours, especially during the peak spawning season. I found males may exert more energy than females during the spawning season, with males showing higher activity rates and longer periods spent on spawning shoals than females, in spite of similar daily travel distances between sexes. Overall, females performed more linear movements over the course of the spawning season suggesting a searching behaviour, while males were less persitent and more random in there movements. Our findings challenge the conventional role of wind as a predominant predictor of lake trout spawning site quality. We propose that the unpredictable nature of wind and abundance of suitable habitat may favour lake-wide spawning by lake trout as a bet-hedging strategy in northern lakes with limited fetch.

##Introduction  

  The evolution of animal mating systems is mainly influenced by sexual selection [@andersson_sexual_1994], parental care [@trivers_parental_1972] and the spatial temporal distribution of resources and mates [@emlen_ecology_1977]. In many species, males optimize their fitness by mating with multiple females. By contrast, the optimal mating rate of females is limited by the production of progeny per mating event [@bateman_intra_1948]. This high energetic investment in gametes by females relative to males often results in conventional sex roles, whereby females provide parental care and males compete for access to females [@kokko_parental_2008]. Because parental care is uncommon in fishes [only found in 21% of families; @blumer_bibliography_1982], their mating systems are ultimately shaped by the distribution of resources necessary for each sex to ensure successful reproduction [@emlen_ecology_1977].  
  
Salmonine fishes (salmons, trouts and chars) are typically observed in a site-based competitive mating system where males compete for access to females, which is thought to be the limiting resource [@gross_sunfish_1984], while females compete for territories to establish and prepare their spawning sites [@fleming_pattern_1998]. Females construct a nest (a series of pits termed a redd) and deposit eggs that are fertilized externally by one or more males. Limited numbers of suitable spawning sites and variability in spawning site quality results in females becoming very selective when choosing their nest sites [@blanchfield_relative_2005; @esteve_observations_2005; @degaudemar_sexual_1998]. Males will fight for proximity to females, with body size and exagerated body shape being the main factors in establishing dominance hierarchies [@fleming_breeding_1994; @quinn_effects_1994]. Paternity tests have shown that the closest male to the female generally fertilizes the greatest proportion of eggs [@blanchfield_breeding_2003; @mjolnerod_mate_1998], thus proximity to females increases fertilization success. As a result, male reproductive success is maximized quantitatively, by mating with as many females as possible; whereas female reproductive success is maximized qualitatively, by choosing high quality nest sites and males  [@degaudemar_sexual_1998].  

  The reproductive behaviours of lake trout (*Salvelinus namaycush*) sharply contrasts the typical salmonine mating system described above [@gunn_spawning_1995; @esteve_lake_2008; @muir_lake_2012]. Spawning usually takes place in lakes on waveswept shoals, and unlike all other salmonines, lake trout do not provide any parental care [i.e., no redd is constructed by the female; @martin_lake_1980]. Eggs are spawned directly onto clean substrate, where they fall into interstitial spaces and incubate for several months before emerging [@royce_breeding_1951]. Lake trout also do not displayovert male-male agonistic behaviour [@royce_breeding_1951; @gunn_spawning_1995; @esteve_lake_2008], a predominant behavioural characteristic of the Salmoninae subfamily [@esteve_observations_2005]. Further, females do not show territorial behaviour (i.e. redd defence) or obvious mate selection [@esteve_lake_2008]. However, some similarities do occur in mating behaviour between lake trout and other salmonines, including males arriving earlier on breeding grounds and staying longer than females each year [@miller_observations_1948; @royce_breeding_1951; @martin_lake_1980; @muir_lake_2012]. In the absence of a site-based competitive mating system and parental care, mating system theory predicts that female reproductive success is driven by habitat quality and mate selection, whereas male reproductive success should be driven by spawning frequency [@degaudemar_sexual_1998]. 
  
  Lake trout (*Salvelinus namaycush*) are a long lived, iteroparous species that typically spawn in lakes during the fall [@gunn_spawning_1995]. Migration from the offshore summer refuge onto nearshore spawning shoals generally coincides with surface water temperatures declining to 12&deg; C or lower [@redick_review_1967]. Preferred spawning habitat is selected along exposed shorelines off points, islands or on mid-lake shoals containing clean substrate including pebble and cobble mainly 3-15 cm in diameter [@martin_lake_1980; See Chapter 2]. Spawning predominately occurs during night time [@gunn_spawning_1995] but has been observed during the day [@esteve_lake_2008; @muir_lake_2012; @binder_new_2015]. Males often constitute 60–85% of annual spawning populations,,, resulting in highly skewed sex ratios on spawning shoals. This imbalance is thought to be a product of earlier maturation, early arrival and longer duration on spawning grounds and the increased likelihood of males spawning each year [@miller_observations_1948; @eschmeyer_reproduction_1955; @martin_reproduction_1957].

  Lake trout reproductive behaviour is relatively understudied when compared to other salmonines, and what literature exists is largely from the southern extent of it's geographic range [@muir_lake_2012]. In general, information on the timing of movement onto spawning shoals and the degree of movements among spawning sites for males and females at seasonal and daily scales remains an important knowledge gap. Furthermore, few studies have investigated reproductive strategies of lake trout [@esteve_lake_2008; @muir_lake_2012]. 
  
The objectives of this chapter is three-fold: (i) to determine when and where lake trout spawn in a typical northern boreal lake; (ii) to describe sex-specific timing and movements of lake trout over the duration of a spawning season; and (iii) to test whether male and female lake trout employ a bet-hedging strategy by spawning on multiple sites. Using data collected from a whole-lake acoustic telemetry array, I will present the movements of five male and six female lake trout over the course of the 2013 spawning season. A spatial temporal clustering analysis was used to determine when and where lake trout formed spawning clusters. These data reveal new insights into the reproductive strategies of male and female lake trout and further our knowledge on spawning site use for managing this iconic Canadian fish species.
  
##Methods  

###Study Site

  The study was conducted in Alexie Lake (62&deg; 40' N, 114&deg; 05' W), located approximately 30 km northeast of Yellowknife, Northwest Territories, Canada (Fig. 3.1). Alexie Lake is a medium size (area: 402 ha, maximum depth: 32 m, mean depth: 11.7 m) oligotrophic lake, possessing 35 islands and a shorelength of 28 km. Stratifcation occurs during the summer until September when the lake becomes isothermal [@healey_experimental_1973; @cott_diel_2015]. Alexie Lake contains four large bodied species; lake trout, northern pike (*Esox lucius*), lake whitefish (*Coregonus clupeaformis*) and burbot (*Lota lota*) as well as numerous prey species [see @cott_food_2011]. Alexie Lake is a designated research lake that is closed to fishing.

<!--![**Fig. 3.1: ** Alexie Lake is located 30 km Northeast of Yellowknife, NWT. The bathymetric contours (blue lines) are at 10 m increments.](alexie_fig.tif)-->

###Lake water temperature

  Water temperature in Alexie Lake was recorded hourly on HOBO&reg; Pendant&trade; temperature loggers (Onset&reg; Computer Corporation, Borne, MA) set in the deepest basin at 0.5 m, at 1 m depth intervals from 1 m to 20 m and at 25 m and 30 m below the water surface. Mean daily temperature for each depth was calculated followed by a spline interpolation to obtain temperatures for every 0.1 m depth interval from lake surface to bottom. Interpolated daily temperature profiles were used to estimate water temperatures occupied by lake trout impanted with pressure-sensing telemetry transmitters [see@guzzo_resource_; @plumb_performance_2009], and determine the role of water temperature in the annual timing of spawning
```{r load_temp, echo=FALSE, cache=TRUE}
load("/Users/DavieCal/Documents/R-Directory/Masters Thesis/Reproductive Strategies 2/temp.RData")
```

###Lake bathymetry

  A high-resolution hydroacoustic survey was conducted in June 2012 (Milne Technologies, Keene,ON, Canada) producing a detailed bathymetric raster (2.5 x 2.5 m). Hydroacoustic data were collected using a 120 kHz Simrad EK60 7.0&deg; split-beam echo sounder system following a systematic parallel survey design with transects 25 m apart [for details see @cott_diel_2015]. These data were used to determine the bathymetric depth for each lake trout telemetry position (see below). This in conjunction with the pressure sensor reading (see below) allows the determination of fish depth from the lake bottom. 

###Fish acoustic telemetry
```{r fish_attributes, echo=FALSE}
#lake trout fork length and weight range attributes
attributes<-read.csv("~/Documents/Masters Thesis/Data/Tagged Fish/Fish_attributes.csv")
#tags of interest
tag<-c("LT-31","LT-32","LT-33","LT-34","LT-35","LT-37","LT-38","LT-40","LT-41","LT-43","LT-44")
attributes<-attributes[attributes$Fish.ID %in% tag,]

male_fl<-summary(attributes$FL..mm.[attributes$Sex=="M"])
male_wt<-summary(attributes$Wt..g.[attributes$Sex=="M"])

female_fl<-summary(attributes$FL..mm.[attributes$Sex=="F"])
female_wt<-summary(attributes$Wt..g.[attributes$Sex=="F"])
```
  Five male (fork length range: `r male_fl[["Min."]]`--`r male_fl[["Max."]]` mm, weight range: `r male_wt[["Min."]]`--`r male_wt[["Max."]]` g) and six female lake trout (fork length range: `r female_fl[["Min."]]`--`r female_fl[["Max."]]` mm, weight range: `r female_wt[["Min."]]`--`r female_wt[["Max."]]` g) were studied during the 2013 spawning season (approximately the month of September). All lake trout were captured by angling with a barbless lure in the spring (June 12--13, 2013) when surface water temperature was <15&deg; C. Once captured lake trout were brought to shore in holding containers, lightly anesthetized in a solution of 90 mg L^-1^ Tricaine Methanesulfonate buffered with sodium bicarbonate and surgically implanted with a coded acoustic transmitter with pressure (for depth) and accelerometer (for activity) sensors that  that randomly emitted an acoustic signal every 80--160 s (V13AP-1L, tail beat algorithm; VEMCO, Ltd., Bedford, NS. Canada). See Blanchfield et al. [-@blanchfield_relative_2005] for details on surgical procedures. Transmitters were 13 mm in diameter, 44 mm in length, and weighed 6 g in water. Each transmitter's depth sensor  was individually calibrated at 4 m depth intervals from the surface to lake bottom in Alexie Lake, and were accurate to &plusmn; 1.7 m. 

  The spatial positions of lake trout implanted with acoustic transmitters were monitored using a VEMCO Positioning System (VPS; VEMCO Ltd.) --- an fine-scale positioning system used for tracking aquatic animals. The VPS uses hyperbolic positioning, also known as multilateration, which measures the time difference on arrival (TDOA) of a signal from a transmitter at three or more time-synchronized receivers [@smith_understanding_2013]. An array of 72 underwater omni-directional acoustic receivers (VEMCO, VR2W) with a mean distance of 201.4 m (range: 120.4m--311.9m) between receivers, was deployed in Alexie Lake to monitor lake trout during spawning behaviour 2013 (see Appendix for further details). 

###Acoustic telemetry data filter 

  Prior to analyses telemetry positions were filtered to ensure that only the highest quality data were used to examine lake trout spawning behaviour. First, positions with successive timestamps less than the minimum interval (80 s) for each individual fish were removed as these positions were assumed to be false positions. Second, lake trout positions greater than 2.5 m (size of one raster grid) outside of the lakeshore were removed. The additional 2.5 m buffer outside of the lakeshore was applied because lake trout are known to spawn in the nearshore area of Alexie Lake (see Chapter 2) which is also a region of increased telemetry error (see Appendix), thus the additional buffur retains valuable spawning data. Third, lake trout depth data had to be greater than or within 2 m of the lake bottom to encompass pressure sensor error (&plusmn;1.7 m). Fourth, only successive positions less than realistic lake trout swimming speeds [<1.2 m s^-1^; @font_behavioural_2015] were included. Finally, I conducted stationary tag trials using VEMCO V13-1P tags with random signal transmission between 10 s to 30 s on known nearshore spawning shoals. I calculated the twice the distance root mean squared (2DRMS) for each estimated stationary tag position, relative to their known location in order to develop a relationship between measured error and hyperbolic positioning error (HPE), a unitless estimate of positioning error provided for each spatial position estimated by the VPS (See Appendix for further details). Using the relationship between 2DRMS and HPE for stationary tag data, I determined that removing all data with an HPE >106 would result in 95% confidence that our data has measurment error less than 27 m, which was determined to be sufficient for our cluster analyses (see Appendix). After filtering, 87% of the original dataset was retained.

###Movement metrics
 
The spawning season was defined as the period between August 30 (day after the 15&deg; C isocline broke; see below)-- September 30 (estimated end date; but see below). Using movement data collected during the spawning season, I calculated for each day, the (1) mean distance travelled (m), (2) mean acceleration (m s^-2^) for energy expenditure or activity of males and females and (3) daily persistence index (PI) or mean cosine of turning angles (unitless) to determine the linearity of movements. The PI ranges from -1 to 1, where a PI = 1 corresponds to linear or directed movement, PI = 0 corresponds to uncorrelated and therefore more tortuous movement, and PI = -1 corresponds to movements that oscillate back and forth [ @laidre_females_2012]. 

###Spatial temporal clustering analysis

Spatial temporal cluster analysis was performed using the density-based clustering algorithm, Spatial Temporal Density-Based Spatial Clustering of Applications with Noise [ST-DBSCAN; @birant_st_2007] to determine areas in Alexie Lake where individual lake trout remained for significant periods of time during the spawning season (i.e. during spawning behaviour, feeding behaviour, resting, etc.). The spatial temporal clustering algorithm requires data points with three required fields: x coordinate, y coordinate and timestamp. Additionally four parameters are required to create clusters: epsilon distance---the Euclidean distance parameter for spatial attributes; epsilon time---the Euclidean distance parameter for time attributes; minimum points-–-the minimum number of points required to create a cluster; and $\Delta$$\epsilon$-–-the acceptable change in time between the current cluster mean and the potential new addition to the cluster. Model parameters were estimated using a simple heuristic described by Ester et al. [-@ester_density_1996] and Birant and Kut [-@birant_st_2007] with my data (See Appendix for further details on parameter estimation). Parameters estimates were as follows: epsilon distance = 41 m; epsilon time = 3254 s; minimum points = 9; and $\Delta$$\epsilon$ = 2 SD. Details on the refinement of these data to determine whether clusters represent spawning activity are described below.  

###Spawning activity determined by cluster filtering

```{r load_data, echo=FALSE,include=FALSE, cache=FALSE}
load("/Users/DavieCal/Documents/R-Directory/Masters Thesis/Reproductive Strategies 2/fishspawn.RData")

load("/Users/DavieCal/Documents/R-Directory/Masters Thesis/Reproductive Strategies 2/SpawnClust.RData")
## calculate the density for all points to attain bandwidth for KDE 

#bandwidth estimaters, bgridsize = 27 because 2DRMS is 27
#f.bw<-hscv(fishspawn$SHOREDIST[fishspawn$SEX=="F"],bgridsize=27)
#m.bw<-hscv(fishspawn$SHOREDIST[fishspawn$SEX=="M"],bgridsize=27)


require(mixtools)
fit<-normalmixEM(fishspawn$SHOREDIST[fishspawn$SEX=="F" & !is.na(fishspawn$SHOREDIST)],lambda=c(0.10,0.30,0.30,0.20,0.10),mu=c(30,80,140,190,250),sigma=c(5,15,15,25,25))
peaks_dist_F<-fit$mu


fit<-normalmixEM(fishspawn$SHOREDIST[fishspawn$SEX=="M" & !is.na(fishspawn$SHOREDIST)],lambda=c(0.2,0.2,0.30,0.20,0.1),mu=c(15,60,140,210,300),sigma=c(5,20,20,40,40))
peaks_dist_M<-fit$mu

fit<-normalmixEM(fishspawn$DEPTH[fishspawn$SEX=="F" & !is.na(fishspawn$DEPTH)],lambda=c(0.40,0.60),mu=c(2,12),sigma=c(2,5))
peaks_depth_F<-fit$mu

fit<-normalmixEM(fishspawn$DEPTH[fishspawn$SEX=="M" & !is.na(fishspawn$DEPTH)],lambda=c(0.40,0.60),mu=c(1,12),sigma=c(1,3))
peaks_depth_M<-fit$mu

fit<-normalmixEM(fishspawn$BOTTOMDEPTH[fishspawn$SEX=="F" & !is.na(fishspawn$BOTTOMDEPTH)],lambda=c(0.2,0.6,0.2),mu=c(5,19,30),sigma=c(5,10,5))
peaks_cont_F<-fit$mu


fit<-normalmixEM(fishspawn$BOTTOMDEPTH[fishspawn$SEX=="M" & !is.na(fishspawn$BOTTOMDEPTH)],lambda=c(0.3,0.7),mu=c(1,17),sigma=c(1,4))
peaks_cont_M<-fit$mu

```

The spatial temporal clustering algorithm is non-biased towards the behaviour of the fish during the cluster formation. This means that when individual fish cluster, it could represent a number of different behaviours, such as spawning, feeding or resting. Therefore, filtering was required to identify potential spawning clusters. Because lake trout spawn nearshore at shallow depths (~2 m) during the month of September in Alexie Lake (see Chapter 2; personal observation), I used distance to shore and depth contour occupied as a filter for identifying potential spawning clusters. This was done as follows. I plotted kernel density estimates (KDE) of distance to nearest shore occupied by male and female lake trout (see Fig. 3.4) which indicated that males remain closer to the shore during the spawning season. It is well documented that males come onto the spawning shoals first and remain longer than females [@miller_observations_1948; @royce_breeding_1951; @martin_lake_1980; @muir_lake_2012], so I selected a 50 m distance to shore criteria based on the visual inspection of the male distance to nearest shore KDE. I also plotted the KDE of contour depth occupied by male and female lake trout (see Fig. 3.3) and found a bimodal distribution for males. The shallower depth contour KDE mode is likely potential spawning males on spawning shoals and therefore I included a criteria of <4 m depth contour. In order to be considered a potential spawning cluster, >50% of the detections within a given cluster must have met both criteria (distance to shore <50 m and depth contour <4m).
  
Sites where potential spawning clusters formed were considered spatially unique when the cluster centers were >70 m apart. This spatial criteria was based on mean distance between the center of known spawning shoals along a 780 m stretch of an Alexie Lake shoreline (site 1; see Chapter 2). Male lake trout often aggregate on spawning shoals but show very little behavioural activity other than when slowly swimming along the spawning shoals edge, however; when one or more females approach a male aggregation, a frenzy of activity ensues [personal observation; esteve_lake_2008; @muir_lake_2012; @binder_new_2015]. To distinguish between these low and high activity periods on spawning shoals, I categorized spawning clusters as low or high activity clusters based on whether the maximum acceleration recorded during a specific cluster was greater or less than 1.3 m s^-2^. This criteria was determined from the mean acceleration of all clusters (spawning and non-spawning; *n* = 4622) + 2 standard deviations (SD). 
  
###Spawning cluster metrics

Spawning clusters were used to determine the timing of lake trout spawning, where the  onset and conclusion of the spawning season were defined as the earliest and latest occurrence of a spawning cluster, respectively. Peak spawning was determined as the earliest and latest date where the daily number of spawning clusters was equal to at least half the maximum number of spawning clusters recorded in a day. To determine if lake trout were spawning in during the day or night, I used the *crepuscule* function from the maptools R package [@maptools_2015] to define periods of day, night, dawn and dusk for Alexie Lake, NWT, during the spawning season. Spawning cluster metrics were calculated for each individual fish and summarized into mean and standard deviation for: (1) frequency of spawning cluster formation, (2) total cumulative duration spent in spawning clusters, (3) duration per cluster formation, (4) total spatially unique sites visited, (5) the distance and (6) duration between subsequent spawning clusters, (7) frequency of cluster formations with low activity (site visits) prior to first high activity cluster (spawn event). Minimum convex polygons (MCP) were calculated using the *mcp* function from the adehabitatHR R package [@adehabitat_2006], to determine the spatial spread of lake trout spawning clusters.
  
###Statistical analysis

All analysis were completed in R V.3.2.1 [@R_2015].
Kernel density estimates (KDE) were computed for distance to shore, depth contour occupied, and depth occupied by lake trout using the *density* function in R. The resulting KDE's were multi-modal, therefore mode means of the KDE's were determined using the *normalmixEM* function, a Expectation-Maximization (EM) algorithm for mixtures of normal distributions, from the *mixtools* R package [@mixtools_2009].
 Mean times and standard deviations of all and high activity spawning clusters were calculated with the *circadian.mean* and *circadian.sd* functions from the *psych* R package [@psych_2015] which calculates the circular mean of circadian data. Using spawning cluster formation start times for all spawning clusters and high activity spawning clusters, I tested whether lake trout spawn during the night or day using a Pearson's Chi-squared test. Using the *loess* function in R, I applied a LOESS (local polynomial regression) curve fitting smoother to the frequency of lake trout spawning by binned one hour segements for all spawning clusters (both high and low activity) and high activity spawning clusers to show the general spawning time trends. A non-parametric Wilcoxon rank sum test was used to test differences in, daily displacement PI, acceleration and site visits prior to spawning between sexes. A log transformation was applied to cluster formation frequency, total duration, mean duration, site visits, distance and duration between consequtive clusters, and minimum convex polygon area data to meet parametric assumptions prior to testing diffferences between sexes, activity state (high or low) and the interaction of sex and activity state using a type III two-way ANOVA.

##Results  

###Spatial distribution of lake trout

  The migration of both male and female lake trout from the deeper offshore to shallower nearshore regions of Alexie Lake corresponded to the breakdown of the 15&deg;C isocline (August 29, 2013; Fig. 2). A clear bimodal depth distribution was apparent for both males and females during the spawning season, with males distributed around mean depths of `r round(peaks_depth_M[1])` m  and `r round(peaks_depth_M[2])` m and females distributed around mean depths `r round(peaks_depth_F[1])` m and `r round(peaks_depth_F[2])` m. Telemetry-tagged male lake trout occupied areas of Alexie Lake where mean lake depths were `r round(peaks_cont_M[1])` m and `r round(peaks_cont_M[2])` m while females occupied mean lake depths of `r round(peaks_cont_F[1])` m, `r round(peaks_cont_F[2])` m, and `r round(peaks_cont_F[3])` m (Fig. 3).


```{r figure_2, echo=FALSE, message=FALSE, fig.height=7, fig.width=6, fig.cap="**Fig. 2:** Depth distribution of males (A) and females (B) during the spawning season in Alexie Lake, NWT. The black line represents the 15&deg; C isocline. Points were given transparency value of 5%, therefore the darker the colour the higher the density of points. No data is available between September 25-27, 2013, while the VPS was removed from the lake for download"}


#        add depth distibution on end of plots      #

par(mfrow=c(2,1))
par(mar=c(0,4.1,4.1,2.1))
#5.1 4.1 4.1 2.1
#males
plot(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15" & fish$DATETIME<="2013-10-15",],
     main="",
     xlab="",
     ylab="",
     type="n",
     bty="n",
     axes=FALSE,
     ylim=c(32,0))

axis(2,las=1)
box()

points(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15"  & fish$DATETIME<="2013-10-15" & fish$SEX=="M",],
       col=rgb(0,0,1,0.05),pch=15)
lines(depth15$Date,depth15$Depth,lwd=4)
text(as.numeric(as.POSIXct("2013-08-15", tz="MST")),0.5,"(A)", font=2)

#females
par(mar=c(5.1,4.1,0,2.1))
plot(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15" & fish$DATETIME<="2013-10-15",],
     main="",
     xlab="",
     ylab="",
     type="n",
     bty="n",
     axes=FALSE,
     ylim=c(32,0))

axis.POSIXct(1,at=c(as.POSIXct("2013-08-15", tz="MST"),as.POSIXct("2013-09-01",tz="MST") ,as.POSIXct("2013-09-15", tz="MST"),as.POSIXct("2013-10-01", tz="MST"),as.POSIXct("2013-10-15", tz="MST")),
             format="%d %b")
axis(2,las=1)
box()

points(DEPTH~DATETIME,data=fish[fish$DATETIME>="2013-08-15"  & fish$DATETIME<="2013-10-15" & fish$SEX=="F",],
       col=rgb(1,0,0,0.05),pch=15)
lines(depth15$Date,depth15$Depth,lwd=4)
text(as.numeric(as.POSIXct("2013-08-15", tz="MST")),0.5,"(B)", font=2)
mtext("Depth (m)", side=2, line=3, adj=1.2)


```


```{r figure_3, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="**Fig. 3:** Kernel density estimate of bottom contour depth occupied by male (blue) and female (red) lake trout in Alexie Lake, NWT, between August 30 and September 30, 2013."}

#Depth Contour

## calculate the density - don't plot yet
densFemale_Spawn <- density(fishspawn$BOTTOMDEPTH[fishspawn$SEX=="F"],from=0)
densMale_Spawn<- density(fishspawn$BOTTOMDEPTH[fishspawn$SEX=="M"],from=0)

max_cont_F<-densFemale_Spawn$x[densFemale_Spawn$y==max(densFemale_Spawn$y)]
max_cont_M<-densMale_Spawn$x[densMale_Spawn$y==max(densMale_Spawn$y)]



## calculate the range of the graph
xlim <- range(0,32)
ylim <- range(0,0.1)
#pick the colours
FemaleCol <- rgb(1,0,0,0.5)
MaleCol <- rgb(0,0,1,0.5)
## plot the carrots and set up most of the plot parameters
plot(densMale_Spawn, xlim = xlim, ylim = ylim, 
     xlab = 'Depth Contour',
     main = "",
     type="n",
     bty="n",
     axes=FALSE
)

#put our density plots in
polygon(c(densMale_Spawn$x,rev(densMale_Spawn$x)),c(densMale_Spawn$y,rep(0,512)), density = -1, col = MaleCol, border = "blue")
polygon(c(densFemale_Spawn$x,rev(densFemale_Spawn$x)),c(densFemale_Spawn$y,rep(0,512)), density = -1, col = FemaleCol, border = "red")


# Draw our own x-axis
axis(1,pos=0)
axis(1,at=xlim,labels=c("",""),pos=0,lwd.ticks=0)

axis(2,las=1,pos=0)


## add a legend in the corner
legend("topright",c('Male', "Female"),
       fill = c(MaleCol,FemaleCol), bty = 'n',
       border = NA)

#abline(v=4)
```

  The spatial distribution of lake trout in Alexie Lake differed between sexes during the spawning season. Males had a higher density of detections in the nearshore region than females (Fig. 4 and 5). Males exhibited a clear peaked mode at a mean distance from shore of `r round(peaks_dist_M)[1]` m, and a second distinct mode at a mean distance from shore of `r round(peaks_dist_M)[3]` m (Fig. 4). Female lake trout exhibited three major modes at mean distances from shore of `r round(peaks_dist_F)[1]` m, `r round(peaks_dist_F)[2]` m, and `r round(peaks_dist_F)[3]` m, but none showed the same high density as males (Fig. 4). 
  

```{r figure_4, echo=FALSE, message=FALSE, fig.cap="**Fig. 4:** Kernel density estimate (KDE) of  distance to nearest shore occupied by male (blue) and female (red) lake trout in Alexie Lake, NWT, between August 30 and September 30, 2013."}
#Nearshore
densFemale_Spawn <- density(fishspawn$SHOREDIST[fishspawn$SEX=="F"],from=0)
densMale_Spawn<- density(fishspawn$SHOREDIST[fishspawn$SEX=="M"],from=0)
max_dist_F<-densFemale_Spawn$x[densFemale_Spawn$y==max(densFemale_Spawn$y)]
max_dist_M<-densMale_Spawn$x[densMale_Spawn$y==max(densMale_Spawn$y)]

## calculate the range of the graph
xlim <- range(0,450)
ylim <- range(0,0.012)
#pick the colours
FemaleCol <- rgb(1,0,0,0.5)
MaleCol <- rgb(0,0,1,0.5)
## plot the carrots and set up most of the plot parameters
plot(densMale_Spawn, xlim = xlim, ylim = ylim, 
     xlab = 'Distance to shore',
     main="",
     type="n",
     bty="n",
     axes=FALSE
)

#put our density plots in
polygon(c(densMale_Spawn$x,rev(densMale_Spawn$x)),c(densMale_Spawn$y,rep(0,512)), density = -1, col = MaleCol, border = "blue")
polygon(c(densFemale_Spawn$x,rev(densFemale_Spawn$x)),c(densFemale_Spawn$y,rep(0,512)), density = -1, col = FemaleCol, border = "red")


# Draw our own x-axis
axis(1,pos=0)
axis(1,at=xlim,labels=c("",""),pos=0,lwd.ticks=0)

axis(2,las=1,pos=0)


## add a legend in the corner
legend("topright",c('Male', "Female"),
       fill = c(MaleCol,FemaleCol), bty = 'n',
       border = NA)
```

  
  
```{r figure_5, echo=FALSE, message=FALSE, fig.cap="**Fig. 5:** Male (blue), and female (red) lake trout positions during the 2013 spawning season in Alexie Lake, NWT.  Points were given transparency value of 10%, therefore the darker the colour the higher the density of points"}

require(maptools)

plot(shore_outline)
points(UTM.Y~UTM.X,data=fishspawn[fishspawn$SEX=="F",],col=rgb(1,0,0,0.1),pch=19,cex=0.1)
points(UTM.Y~UTM.X,data=fishspawn[fishspawn$SEX=="M",],col=rgb(0,0,1,0.1),pch=19,cex=0.1)

```


```{r cluster_dist, echo=FALSE, message=FALSE}
require(adehabitatHR)
xy<-SpatialPoints(cbind(spawn$CLUST.X,spawn$CLUST.Y))
all_mcp<-mcp(xy,percent=100,unout="km2")

```

###Spatial distribution of spawning clusters

A total of `r dim(spawn)[1]` lake trout spawning clusters were estimated using the spawning cluster analysis, of which `r dim(spawn[spawn$SPAWN==1,])[1]` (`r round(dim(spawn[spawn$SPAWN==1,])[1]/dim(spawn)[1]*100)`%) were classified as high activity clusters and more likely to have spawning behaviour taken place (Fig. 6). Spawning clusters were found throughout the nearshore region of Alexie Lake, covering an area of `r round(all_mcp$area,digits=2)` km^2^ (based on MCP).

```{r Figure 6, echo=FALSE,message=FALSE,fig.cap="**Fig. 6:** Lake trout spawning clusters in Alexie Lake, NWT during the 2013 spawning season. The smaller red circles represent the center of high activity clusters and larger blue circles represent the center of low activity clusters"}
#require(adehabitatHR)
plot(shore_outline)
points(CLUST.Y~CLUST.X,data=spawn[spawn$SPAWN==0,],pch=19,col=rgb(0,0,1,1))

#points(CLUST.Y~CLUST.X,data=spawn[spawn$SPAWN==1,],pch=19,col=rgb(1,0,0,0.5))
points(CLUST.Y~CLUST.X,data=spawn[spawn$SPAWN==1,],pch=19,cex=0.5,col=rgb(1,0,0,1))

```



```{r cluster data, echo=FALSE, message=FALSE, cache=FALSE}

peakstart=min(spawn.num$DATE[spawn.num$FREQ>=(max(spawn.num$FREQ)/2)])
peakend=max(spawn.num$DATE[spawn.num$FREQ>=(max(spawn.num$FREQ)/2)])
spawnstart=min(spawn.num$DATE)
spawnend=max(spawn.num$DATE)
peakdays=spawn.num$DATE[spawn.num$FREQ==max(spawn.num$FREQ)]
malepeak=spawn.num$DATE[spawn.num$FREQ_MALE==max(spawn.num$FREQ_MALE)]
femalepeak=spawn.num$DATE[spawn.num$FREQ_FEMALE==max(spawn.num$FREQ_FEMALE)]
maledates=spawn.num$DATE[which(spawn.num$FREQ_MALE>0)]
femaledates=spawn.num$DATE[which(spawn.num$FREQ_FEMALE>0)]

mdates<-unique(spawn$DATE[spawn$SEX=="M"])
mdates<-sort(mdates)
mdates<-data.frame(DATE=mdates,
                   DIFF=c(NA,difftime(mdates[2:length(mdates)],mdates[1:length(mdates)-1],units="days")))
mdates<-mdates[-1,]
#assume clusters during download
mdays<-difftime(as.POSIXct("2013-09-24",tz="MST"),mdates$DATE[1],units="days")



fdates<-unique(spawn$DATE[spawn$SEX=="F"])
fdates<-sort(fdates)
fdates<-data.frame(DATE=fdates,
                   DIFF=c(NA,difftime(fdates[2:length(fdates)],fdates[1:length(fdates)-1],units="days")))
fdates<-fdates[-1,]
fdays<-difftime(as.POSIXct("2013-09-16",tz="MST"),as.POSIXct("2013-09-08",tz="MST"),units="days")


mdays_range<-range(tapply(spawn$DATE[spawn$SEX=="M"],spawn$TRANSMITTER[spawn$SEX=="M"], function(x) length(unique(x))),na.rm=TRUE)

fdays_range<-range(tapply(spawn$DATE[spawn$SEX=="F"],spawn$TRANSMITTER[spawn$SEX=="F"], function(x) length(unique(x))),na.rm=TRUE)

```


###3.3.3 Timing of spawning  


  All five males formed spawning clusters, while five out six females investigated formed spawning clusters; LT-32 never formed a spawning cluster in 2013 and was designated a non-spawning female (see Fig. 7). Lake trout showed a wide range in timing of spawning cluster formation, with males spending `r mdays_range[1]`--`r mdays_range[2]` days on spawning shoals compared to only `r fdays_range[1]`--`r fdays_range[2]` days by females. Most individuals (80%) were present during peak spawning, but two individuals, LT-31 (male) and LT-38 (female) only formed spawning clusters early in the season prior to peak spawn and did not return during or post peak spawn. 

```{r figure_7, echo=FALSE,warning=FALSE,fig.cap="**Fig. 7:** Spawning cluster formation by individual lake trout in Alexie Lake, NWT. Squares represent one or more spawning cluster formed on the date for males (blue) and females (red) and vertical black lines represent the onset and conclusion of peak spawning."}
spawn$SEX<-as.factor(spawn$SEX)
palette(c("red","blue"))
plot(spawn$DATE,spawn$TRANSMITTER,pch=15, cex=0.7,axes=FALSE,xlim=c(min(dateseq),max(dateseq)),ylab="Fish ID",xlab="",col=spawn$SEX)
axis(2,at=c(1:length(tag)),labels=tag,las=1)
axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d")
lines(c(peakstart-12*60*60,peakstart-12*60*60),c(0.5,11.5))
lines(c(peakend+12*60*60,peakend+12*60*60),c(0.5,11.5))
box(lwd=2)
```

  Based on spawning clusters, the onset of lake trout spawning in 2013 began `r format(spawnstart,"%B %e")` coinciding with a surface temperature of `r round(temp$Temp[temp$Date==spawnstart & temp$Depth==0], digits=1)`&deg;C (Fig. 8). Spawning occurred over a `r as.numeric(difftime(spawnend,spawnstart,units="days"))` d period and concluded on `r format(spawnend,"%B %e")` when surface temperature was `r round(temp$Temp[temp$Date==spawnend & temp$Depth==0], digits=1)`&deg;C. Peak spawn began on `r format(peakstart,"%B %e")`, and lasted `r as.numeric(difftime(peakend,peakstart,units="days"))` d, concluding on `r format(peakend,"%B %e")`. Surface temperature was `r round(temp$Temp[temp$Date==peakstart & temp$Depth==0], digits=1)`&deg;C at the start of peak spawn and dropped `r round(temp$Temp[temp$Date==peakstart & temp$Depth==0]-temp$Temp[temp$Date==peakend & temp$Depth==0], digits=1)`&deg;C over the duration of peak spawn. The maximum number of spawning clusters (`r max(spawn.num$FREQ)`) occured on `r format(peakdays[1],"%B %e")`. The majority of all spawning clusters (`r round(dim(spawn[spawn$DATE>=peakstart & spawn$DATE<=peakend,])[1]/dim(spawn)[1]*100)`%) and high activity spawning clusters (`r round(dim(spawn[spawn$DATE>=peakstart & spawn$DATE<=peakend & spawn$SPAWN==1,])[1]/dim(spawn[spawn$SPAWN==1,])[1]*100)`%) occurred during peak spawn. The frequency distributions of total, high activity, and low activity spawning cluster were similar (Fig. 8). Peak spawning for male and female lake trout occurred at similar dates, with peaks occurring on `r format(malepeak, "%B %e")` and `r format(femalepeak, "%B %e")`, respectively  (Fig. 9). Males collectively remained on the spawning shoals for maximum of `r as.numeric(mdays)` consecutive days while females collectively spent a maximum of `r as.numeric(fdays)` consecutive days.
 
```{r figure_8, echo=FALSE, message=FALSE, fig.height=7, fig.width=7, fig.cap="**Fig. 8:** Frequency distribution of spawning clusters over the duration of the 2013 spawning season in Alexie Lake, NWT. Light grey bars represent the cumulative frequency of spawning clusters (both low and high activity), dark grey bars represent high activity spawning clusters, the horizontal black line represents daily surface water temperature (&deg;C) and the vertical black lines represent the onset and conclusion of peak spawning."}


plot(FREQ~DATE,data=spawn.num,col=rgb(0,0,0,0.3),ylim=c(0,25),type="h",xlim=c(min(dateseq),max(dateseq)),main="",ylab="Frequency & Temperature ",xlab="",axes=FALSE,lwd=7,lend=1)

lines(SPAWN~DATE,data=spawn.num,col=rgb(0,0,0,0.7),type="h",lwd=7,lend=1)

lines(c(peakstart-12*60*60,peakstart-12*60*60),c(0,25))
lines(c(peakend+12*60*60,peakend+12*60*60),c(0,25))

axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d",pos=0)
axis(2,las=1,pos=min(dateseq))


lines(Temp~Date,data=temp[temp$Depth==0 & temp$Date >="2013-08-15" & temp$Date <="2013-11-01",],lwd=2)

```

```{r figure_9, echo=FALSE, message=FALSE, fig.height=7, fig.width=7, fig.cap="**Fig. 9:** Frequency distributions of male (blue) and female (red) spawning clusters over the duration of the 2013 spawning season in Alexie Lake, NWT. Light bars represent daily cumulative frequency of spawning clusters (both low and high activity), dark bars represent high activity spawning clusters and vertical black lines represent the onset and conclusion of peak spawning."}

plot(FREQ~DATE,data=spawn.num,col=rgb(0,0,0,0.3),ylim=c(-15,15),type="n",xlim=c(min(dateseq),max(dateseq)),main="",ylab="Frequency",xlab="",axes=FALSE)

lines(FREQ_MALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(0,0,1,0.3))
lines(SPAWN_MALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(0,0,1,0.7))
lines(-FREQ_FEMALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(1,0,0,0.3))
lines(-SPAWN_FEMALE~DATE,data=spawn.num,type="h",lwd=7,lend=1,col=rgb(1,0,0,0.7))

lines(c(peakstart-12*60*60,peakstart-12*60*60),c(-10,15))
lines(c(peakend+12*60*60,peakend+12*60*60),c(-10,15))

axis.POSIXct(1,at=c(dateseq[1],dateseq[79]),labels=FALSE,lwd.ticks=0,pos=0)

ylabs <- pretty(c(-10,15))
axis(2, at = ylabs, labels = abs(ylabs), pos=min(dateseq),las=1)


axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d",pos=-10)

```

```{r figure_10_stats,echo=FALSE, message=FALSE, }
require(lubridate)
ltspawnclust3$HOURPLOT<-hour(ltspawnclust3$STARTTIME)
require(psych)
timemean<-circadian.mean(ltspawnclust3$HOUR)
minutemean<-floor((timemean-floor(timemean))*60)
timemean<-paste(floor(timemean),":",minutemean,sep="")
timesd<-circadian.sd(ltspawnclust3$HOUR)$sd

timemean_s<-circadian.mean(ltspawnclust3$HOUR[ltspawnclust3$SPAWN==1])
minutemean<-floor((timemean_s-floor(timemean_s))*60)
timemean_s<-paste(floor(timemean_s),":",minutemean,sep="")
timesd_s<-circadian.sd(ltspawnclust3$HOUR[ltspawnclust3$SPAWN==1])$sd


#y <- sqrt(table(factor(ltspawnclust3[ltspawnclust3$SPAWN==1,"HOURPLOT"], levels=0:23)))
#y <- table(factor(ltspawnclust3$HOURPLOT, levels=0:23))
y <- table(factor(ltspawnclust3$HOURPLOT, levels=1:24))
y<-as.data.frame(y)
y$Var1<-as.numeric(y$Var1)

maxtime_all<-y$Var1[y$Freq==max(y$Freq)]
maxtime_all<-paste(maxtime_all,":00--",maxtime_all+1,":00",sep="")
maxfreq_all<-max(y$Freq)

mintime_all<-y$Var1[y$Freq==min(y$Freq)]
mintime_all<-paste(mintime_all,":00--",mintime_all+1,":00",sep="")
minfreq_all<-min(y$Freq)

#y$Var1<-as.numeric(y$Var1)+1
y.smooth<-loess(Freq~Var1,data=y)

y.smooth<-predict(y.smooth, se=T)

y <- table(factor(ltspawnclust3$HOURPLOT[ltspawnclust3$SPAWN==1],levels=1:24))
y<-as.data.frame(y)
y$Var1<-as.numeric(y$Var1)

maxtime_s<-y$Var1[y$Freq==max(y$Freq)]
maxtime_s<-paste(maxtime_s,":00--",maxtime_s+1,":00",sep="")
maxfreq_s<-max(y$Freq)

mintime_s<-y$Var1[y$Freq==min(y$Freq)]
mintime_s<-paste(mintime_s,":00--",mintime_s+1,":00",sep="")
minfreq_s<-min(y$Freq)

#y$Var1<-as.numeric(y$Var1)+1
y.smooth2<-loess(Freq~Var1,data=y)

y.smooth2<-predict(y.smooth2, se=T)

spawntime_all<-table(ltspawnclust3$TIME)/dim(ltspawnclust3)[1]*100
spawntime_s<-table(ltspawnclust3$TIME[ltspawnclust3$SPAWN==1])/dim(ltspawnclust3[ltspawnclust3$SPAWN==1,])[1]*100
a<-chisq.test(table(ltspawnclust3$TIME[ltspawnclust3$SPAWN==1]))

```

  Lake trout showed a significant preference for spawning during the night; `r round(spawntime_s[2])`% of high activity spawning clusters occured at night (Pearson's Chi-squared Test: *X^2^* = `r round(a$statistic,digits=2)`; *p* = `r round(a$p.value,digits=3)`). However, when considering all spawning clusters (both high and low activity clusters), lake trout did not show a significant preference for night spawning (Fig. 3.10). Mean cluster formation occurred at `r timemean` &plusmn; `r round(timesd,digits=1)` h (mean &plusmn; SD), and high activity spawning clusters occurred at `r timemean_s` &plusmn; `r round(timesd_s,digits=1)` h. The maximum number of spawning clusters (*n* = `r maxfreq_all`) occurred between `r maxtime_all`, similar to the greatest number of high activity spawning clusters (*n* = `r maxfreq_s`), which also occured between `r maxtime_s`. The minimum number of spawning clusters (*n* = `r minfreq_all`) occurred between morning hours of `r mintime_all`, and no high activity spawning occurred between the hours of `r mintime_s`. 

```{r figure_10, echo=FALSE, message=FALSE,fig.width=8,fig.height=7,fig.cap="**Fig. 10:** Lake trout spawning cluster formation start times over the duration of the 2013 spawning season (left panel), and summarized as total counts per hour (right panel), in Alexie Lake, NWT. Both low activity (open circles) and high activity (closed circles) are displayed in both panels. Light grey polygons represent dawn and dusk, dark grey polygons represent night time and vertical black lines represent the onset and conclusion of peak spawning. The lines in the right panel represent a LOESS trendline with 95% confidence interval error bars."}
#quartz()
layout(matrix(c(1,1,1,2,1,1,1,2),2,4,byrow=TRUE))
par(mar=c(5.1,4.1,4.1,0),
    oma=c(2,2,0,0))

plot(HOURPLOT~STARTTIME, data=ltspawnclust3, col=rgb(0,0,0,0),pch=19, xlab="Date",ylab="Hour", xlim=c(min(dateseq),max(dateseq)),ylim=c(0.5,24.5), axes=FALSE,xaxs="i",yaxs="i")



axis.POSIXct(1,at=c(dateseq[1],dateseq[18],dateseq[32],dateseq[48],dateseq[62],dateseq[79]),format="%b-%d",pos=0.5)


axis(2,at=seq(0,24,by=4),
     las=1)

#sunrise
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$SUNRISE,rep(0.5,dim(suntimes)[1])), density = -1, col = "grey",border="grey")
#dawn
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$DAWN,rep(0.5,dim(suntimes)[1])), density = -1, col = "darkgrey", border = "darkgrey")


#sunset
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$SUNSET,rep(24.5,dim(suntimes)[1])), density = -1, col = "grey",border="grey")
#dusk
polygon(c(suntimes$DATE,rev(suntimes$DATE)),c(suntimes$DUSK,rep(24.5,dim(suntimes)[1])), density = -1, col = "darkgrey", border = "darkgrey")

points(HOURPLOT~STARTTIME, data=ltspawnclust3,col=rgb(1,0,0,1))
#points(HOUR~STARTTIME, data=ltspawnclust3[ltspawnclust3$SPAWN==0,], pch=19, col=rgb(0,0,1,1))
points(HOURPLOT~STARTTIME, data=ltspawnclust3[ltspawnclust3$SPAWN==1,],pch=19,col=rgb(1,0,0,1))

lines(c(peakstart-12*60*60,peakstart-12*60*60),c(0.5,24.5))
lines(c(peakend+12*60*60,peakend+12*60*60),c(0.5,24.5))

box(lwd=2)


par(mar=c(5.1,0,4.1,1))


plot(y.smooth$fit,y$Var1,pch=19,ylim=c(0.5,24.5), xlim=c(0,16), xlab="", ylab="",type="n",axes=FALSE,xaxs="i",yaxs="i")
axis(1,at=c(0,16),labels=FALSE,lwd.ticks=0,pos=0.5,lwd=2)
axis(2,at=c(0.5,24.5),labels=FALSE,lwd.ticks=0,pos=0,lwd=2)
axis(3,at=seq(0,16,4), pos=24.5,lwd=2)
axis(4,at=c(0.5,24.5),labels=FALSE,lwd.ticks=0,pos=16,lwd=2)
mtext("Count",side=3,line=3,cex=0.6)


arrows(x0=y.smooth$fit+ qt(0.975,y.smooth$df)*y.smooth$se,
       y0=y$Var1,
       x1=y.smooth$fit- qt(0.975,y.smooth$df)*y.smooth$se, 
       angle=90,code=3,length=0.05)


arrows(x0=y.smooth2$fit+ qt(0.975,y.smooth2$df)*y.smooth2$se,
       y0=y$Var1,
       x1=y.smooth2$fit- qt(0.975,y.smooth2$df)*y.smooth2$se, 
       angle=90,code=3,length=0.05)

points( y.smooth$fit,y$Var1,pch=19, col="white")
lines( y.smooth$fit,y$Var1,type="b",lwd=2, col="red")
lines( y.smooth2$fit,y$Var1,type="b",lwd=2, pch=19,col="red")


```


```{r table1_stats, echo=FALSE,warning=FALSE,comment=FALSE}
#format scientific notation
knitr::knit_hooks$set(inline = function(x) {
  knitr:::format_sci(x, 'md')
})

#hist(sqrt(daily$STEPLENGTH))
#shapiro.test(sqrt(daily$STEPLENGTH))
#bartlett.test(sqrt(STEPLENGTH)~SEX,data=daily)
#shapiro.test(daily$PERSINDEX[ !daily$TRANSMITTER=="LT-32"]^())
#hist(daily$STEPLENGTH[!daily$TRANSMITTER=="LT-32"]^(1/2))

#shapiro.test(log(daily$ACCEL[!daily$TRANSMITTER=="LT-32"]))

#hist(log(daily$ACCEL[!daily$TRANSMITTER=="LT-32"]))
#displacement test
dist.test<-wilcox.test(daily$STEPLENGTH[daily$SEX=="M"],daily$STEPLENGTH[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])



#PI test
pi.test<-wilcox.test(daily$PERSINDEX[daily$SEX=="M"],daily$PERSINDEX[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])



#ACCEL test
accel.test<-wilcox.test(daily$ACCEL[daily$SEX=="M"],daily$ACCEL[daily$SEX=="F" & !daily$TRANSMITTER=="LT-32"])



```

###Lake trout spawning behaviour

  Over the course of the spawning season lake trout daily distance travelled was not significantly different between spawning males and females with both sexes moving ~ 7 km per day (table 1). The non-spawning female did show a higher daily displace of `r round(mean(daily$STEPLENGTH[daily$TRANSMITTER=="LT-32"])/1000,digits=2)` (&plusmn;`r round(sd(daily$STEPLENGTH[daily$TRANSMITTER=="LT-32"])/1000,digits=2)`; table 1). Spawning females made more linear movements than males, as indicated by females having a significantly higher persistence index (Wilcoxon Rank Sum Test: *W* = `r pi.test$statistic`; *p* `r ifelse(pi.test$p.value<0.001,"< 0.001",paste("= ",round(pi.test$p.value,digits=3),sep=""))`; table 1). The non-spawning female also made more directed movement consisitent with spawning females but with a greater PI (`r round(mean(daily$PERSINDEX[daily$TRANSMITTER=="LT-32"]),digits=2)` &plusmn; `r round(sd(daily$PERSINDEX[daily$TRANSMITTER=="LT-32"]),digits=2)`; table 1). Spawning males exerted significantly more energy than females over the course of the spawning season in spite of similar daily displacements, exhibiting daily acceleration rates `r round(spawntable$ACCEL[1]-spawntable$ACCEL[2],digits=2)` m s^-2^ greater than females (*W*= `r accel.test$statistic`; *p* `r ifelse(accel.test$p.value<0.001,"< 0.001",paste("= ",round(accel.test$p.value,digits=3),sep=""))`; table 3.1). The non-spawning female also showed similer acceleration as spawning females (`r round(mean(daily$ACCEL[daily$TRANSMITTER=="LT-32"]),digits=2)` &plusmn; `r round(sd(daily$ACCEL[daily$TRANSMITTER=="LT-32"]),digits=2)`; table 1).  
    

**Table 1:** Summary statistics for spawning male, spawning female and non-spawning female lake trout movement metrics during the 2013 spawning season in Alexie Lake, NWT. I report the number of fish (n), means and standard deviations of daily distance travelled (km), persistence index (PI), and acceleration (m s^-2^).

|Sex| n | Daily Distance (km) | PI | Acceleration (m s^-2^) |
|:-:|:-:|:-:|:-:|:-:|
|Spawning Male | 5 |`r round(spawntable$DAILY_DISPLACE[1]/1000,digits=2)` (`r round(spawntable$SD_DAILY_DISPLACE[1]/1000,digits=2)`)  |`r round(spawntable$PI[1],digits=2)` (`r round(spawntable$SD_PI[1],digits=2)`)|`r round(spawntable$ACCEL[1],digits=2)` (`r round(spawntable$SD_ACCEL[1],digits=2)`)|
|Spawning Female | 5  |`r round(spawntable$DAILY_DISPLACE[2]/1000,digits=2)` (`r round(spawntable$SD_DAILY_DISPLACE[2]/1000,digits=2)`)  |`r round(spawntable$PI[2],digits=2)` (`r round(spawntable$SD_PI[2],digits=2)`)|`r round(spawntable$ACCEL[2],digits=2)` (`r round(spawntable$SD_ACCEL[2],digits=2)`)|
|Non-spawning Female | 1  |`r round(mean(daily$STEPLENGTH[daily$TRANSMITTER=="LT-32"])/1000,digits=2)` (`r round(sd(daily$STEPLENGTH[daily$TRANSMITTER=="LT-32"])/1000,digits=2)`) |`r round(mean(daily$PERSINDEX[daily$TRANSMITTER=="LT-32"]),digits=2)` (`r round(sd(daily$PERSINDEX[daily$TRANSMITTER=="LT-32"]),digits=2)`)|`r round(mean(daily$ACCEL[daily$TRANSMITTER=="LT-32"]),digits=2)` (`r round(sd(daily$ACCEL[daily$TRANSMITTER=="LT-32"]),digits=2)`)|


```{r table2_stats_freq, include=FALSE}

require(car)
spawn$SPAWN<-as.factor(spawn$SPAWN)

#total clust
clust<-as.data.frame(table(spawn$SEX,spawn$TRANSMITTER))
clust<-clust[clust$Freq!=0,-2]
m_min<-min(clust$Freq[clust$Var1=="M"])
m_max<-max(clust$Freq[clust$Var1=="M"])

f_min<-min(clust$Freq[clust$Var1=="F"])
f_max<-max(clust$Freq[clust$Var1=="F"])


#clust sex and activity
clust<-as.data.frame(table(spawn$SEX,spawn$SPAWN,spawn$TRANSMITTER))
clust<-clust[clust$Freq!=0,]

#hist(log10(clust$Freq))

#shapiro.test(log1(clust$Freq))
#leveneTest(log(clust$Freq)~clust$Var1*clust$Var2)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

clustmod<-lm(log(Freq)~Var1*Var2,data=clust)



aclust.test<-Anova(clustmod, type=3)

```
Telemetry-tagged male lake trout created a total of `r dim(spawn[spawn$SEX=="M",])[1]` clusters, `r round(dim(spawn[spawn$SEX=="M" & spawn$SPAWN==1,])[1]/dim(spawn[spawn$SEX=="M",])[1]*100)`% of which were high activity clusters. Females on the otherhand created only `r dim(spawn[spawn$SEX=="F",])[1]` clusters, `r round(dim(spawn[spawn$SEX=="F" & spawn$SPAWN==1,])[1]/dim(spawn[spawn$SEX=="F",])[1]*100)`% of which were high activity clusters. Male individuals formed significantly more clusters on average than females (ANOVA: *F*~(1,16)~ = `r round(aclust.test[2,3],digits=2)`; *p* = `r round(aclust.test[2,4],digits=3)`; table 3.2). No significant difference was found between activity states (*p* = `r round(aclust.test[3,4],digits=3)`) or the interaction of sex and activity (*p* = `r round(aclust.test[4,4],digits=3)`).

```{r table2_stats_dur, include=FALSE}

#total duration

totaldur<-aggregate(TOTALTIME~SEX+SPAWN+TRANSMITTER,data=spawn,sum)


#hist(log(totaldur$TOTALTIME))
#shapiro.test(log(totaldur$TOTALTIME))
require(car)
#leveneTest(log(totaldur$TOTALTIME)~totaldur$SEX*totaldur$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

tdurmod<-lm(log(TOTALTIME)~SEX*SPAWN,data=totaldur)

tdurp.test<-Anova(tdurmod, type=3)



#duration per cluster

dur<-data.frame(SEX=spawn$SEX,
                SPAWN=as.factor(spawn$SPAWN),
                TRANSMITTER=spawn$TRANSMITTER,
                DURATION=spawn$TOTALTIME)


#hist(log(dur$DURATION))
#shapiro.test(log(dur$DURATION))
#bartlett.test(log(dur$DURATION)~dur$SEX)

require(car)
#leveneTest(log(dur$DURATION)~dur$SEX*dur$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

durmod<-lm(log(DURATION)~SEX*SPAWN,data=dur)

adurp.test<-Anova(durmod, type=3)



#unique site visits
require(car)
visit<-aggregate(SITE~SEX+SPAWN+TRANSMITTER,data=spawn,function(x) length(unique(x)))

#hist(log(visit$SITE))
#shapiro.test(log(visit$SITE))
#leveneTest(log(visit$SITE)~visit$SEX*visit$SPAWN)

#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

vismod<-lm(log(SITE)~SEX*SPAWN,data=visit)

avisit.test<-Anova(vismod, type=3)

```
  
Cluster duration differed significantly between activity states with lake trout spending `r round(mean(spawn$TOTALTIME[spawn$SPAWN==1])-mean(spawn$TOTALTIME[spawn$SPAWN==0]),digits=2)` more minutes in high activty cluster than low activity clusters (*F*~(1,222)~ = `r round(adurp.test[3,3],digits=2)` ;*p* = `r round(adurp.test[3,4],digits=3)`). No significant differences in cluster duration were found between sexes (*p* = `r round(adurp.test[2,4],digits=2)`) or the interaction of sex and activity (*p* = `r round(adurp.test[4,4],digits=2)`), even though males spent on average `r round(spawntable$TIMEMEAN[1]-spawntable$TIMEMEAN[2],digits=1)` minutes longer in clusters than females (table 3.2). Males in total spent significantly more time in spawning clusters (`r round(spawntable$TOTALTIME[1]/60/24, digits=1)` &plusmn; `r round(spawntable$TOTALTIMESD[1]/60/24, digits=1)` days) compared to females (`r round(spawntable$TOTALTIME[2]/60/24,digits=1)` &plusmn; `r round(spawntable$TOTALTIMESD[2]/60/24,digits=1)` days; *F*~(1,16)~ = `r round(tdurp.test[2,3],digits=2)` ; *p* = `r round(tdurp.test[2,4],digits=2)`). No significant differences in total cluster duration were found between activity state (*p* = `r round(tdurp.test[3,4],digits=2)`) or the interaction of sex and activity (*p* = `r round(tdurp.test[4,4],digits=2)`).
Males visited significantly more spatialy unique spawning sites, visiting `r spawntable$SITEVISIT[1]-spawntable$SITEVISIT[2]` more unique sites on average than females (*F*~(1,16)~ = `r round(avisit.test[2,3],digits=2)`; *p* = `r round(avisit.test[2,4],digits=3)`). No significant differences in number of site visits were found between activity state (*p* = `r round(avisit.test[3,4],digits=2)`) or the interaction of sex and activity (*p* = `r round(avisit.test[4,4],digits=2)`).
```{r table2_stats_betw, include=FALSE}

#dist between

dist<-data.frame(SEX=spawn$SEX,
                SPAWN=as.factor(spawn$SPAWN),
                TRANSMITTER=spawn$TRANSMITTER,
                DISTCLUST =spawn$DISTCLUST,
                TIMECLUST=spawn$TIMECLUST)



#hist(log(dist$DISTCLUST))

require(car)
#leveneTest((log(dist$DISTCLUST))~dist$SEX*dist$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

distmod<-lm(log(DISTCLUST)~SEX*SPAWN,data=dist)

adist.test<-Anova(distmod,type=3)


#time between

#hist(1/(dist$TIMECLUST)^(1/3))

#shapiro.test(1/(dist$TIMECLUST)^(1/3))
#qqnorm(1/((dist$TIMECLUST)^(1/3)))
require(car)
#leveneTest(log(dur$DURATION)~dur$SEX*dur$SPAWN)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))

timemod<-lm(1/(TIMECLUST)^(1/3)~SEX*SPAWN,data=dist)

atime.test<-Anova(timemod, type=3)


#sites prior to spawn
prior<-data.frame(SEX=spawn$SEX,
                SPAWN=as.factor(spawn$SPAWN),
                TRANSMITTER=spawn$TRANSMITTER,
                CLUSTPRIOR =spawn$CLUSTPRIOR)
prior<-aggregate(CLUSTPRIOR~SEX+TRANSMITTER,data=spawn,mean)


prior.test<-wilcox.test(CLUSTPRIOR~SEX,data=prior)





#MCP

smcp<-mcpclust[,-c(3,5,6)]
smcp$SPAWN<-1
names(smcp)<-c("TRANSMITTER","SEX","MCP","SPAWN")
nmcp<-mcpclust[,-c(3,4,6)]
nmcp$SPAWN<-0
names(nmcp)<-c("TRANSMITTER","SEX","MCP","SPAWN")
mcp<-rbind(smcp,nmcp)
mcp$SPAWN<-as.factor(mcp$SPAWN)

#shapiro.test(log(mcp$MCP))
require(car)
#leveneTest(log(MCP)~SEX*SPAWN,data=mcp)
#majic anova code...ANOVA type III
options(contrasts=c("contr.sum","contr.poly"))
areamod<-lm(log(MCP)~SEX*SPAWN,data=mcp)
aarea.test<-Anova(areamod,type=3)


```
Although the distance between consecutive spawning clusters for females was `r round(spawntable$DISTBETW[2],digits=2)- round(spawntable$DISTBETW[1],digits=2)` m greater than males (table 3.2), it did not significantly differ between sexes, activity state or the interaction of sex and spawning state. Similarly, despite the duration between female spawning clusters taking on average `r round(spawntable$TIMEBETW[2],digits=1)-round(spawntable$TIMEBETW[1],digits=1)`  more hours than males, no significant difference was found for duration between consecutive clusters between sexes, activity state or the interaction of sex and activity state. Low activity site visits prior to spawning (high activity site visits) were not found to be significantly different between males and females, each visiting roughly 2 sites prior to spawning (table 3.2). Male spawning site MCP area was on average `r round(spawntable$AREA[1]-spawntable$AREA[2],digits=1)` km^2^ greater than females (table 3.2), but was not found to be significantly different between sexes, activity state or the interaction between sex and activity state. 

**Table 2:** Summary statistics for male and female lake trout spawning cluster data during the 2013 spawning season in Alexie Lake, NWT. For all male and female clusters (both high and low activity clusters) and spawning male and female clusters (high activity clusters), I report the number of fish (n), means and standard deviations of clusters formed, cluster duration, unique cluster sites visited, distance between consecutive clusters (m), duration between subsequent clusters (hr), site visits prior to first spawn and minimum convex polygon area of all cluster points (km^2^).

|Sex| n | Clusters | Duration (min)| Sites | Distance Between Clusters (m)| Duration Between Clusters (hr)| Site Visits Prior to Spawn | MCP Area (km^2^) |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|All Male | 5 |`r round(spawntable$MEANCLUST[1],digits=1)` (`r round(spawntable$MCSD[1],digits=1)`)|`r round(spawntable$TIMEMEAN[1],digits=1)` (`r round(spawntable$TIMESD[1],digits=1)`)|`r round(spawntable$SITEVISIT[1],digits=1)` (`r round(spawntable$SVSD[1],digits=1)`)|`r round(spawntable$DISTBETW[1],digits=2)` (`r round(spawntable$DSD[1],digits=2)`)|`r round(spawntable$TIMEBETW[1],digits=1)` (`r round(spawntable$TSD[1],digits=1)`)|-|`r round(spawntable$AREA[1],digits=1)` (`r round(spawntable$ASD[1],digits=1)`)|
|Low Activity Male| 5 |`r round(spawntable$NOMEANCLUST[1],digits=1)` (`r round(spawntable$NCSD[1],digits=1)`)|`r round(spawntable$NOTIMEMEAN[1],digits=1)` (`r round(spawntable$NOTIMESD[1],digits=1)`)|`r round(spawntable$NOSITEVISIT[1],digits=1)` (`r round(spawntable$NSSD[1],digits=1)`)|`r round(spawntable$NODISTBETW[1],digits=2)` (`r round(spawntable$NDSD[1],digits=2)`)|`r round(spawntable$NOTIMEBETW[1],digits=1)` (`r round(spawntable$NTSD[1],digits=1)`)|-|`r round(spawntable$NOAREA[1],digits=1)` (`r round(spawntable$NASD[1],digits=1)`)|
|High Activity Male| 5 |`r round(spawntable$MEANSPAWN[1],digits=1)` (`r round(spawntable$MSSD[1],digits=1)`) |`r round(spawntable$SPAWNTIMEMEAN[1],digits=1)` (`r round(spawntable$SPAWNTIMESD[1],digits=1)`) | `r round(spawntable$SITESPAWN[1],digits=1)` (`r round(spawntable$SSSD[1],digits=1)`)|`r round(spawntable$SPAWNDISTBETW[1],digits=2)` (`r round(spawntable$SDSD[1],digits=2)`)|`r round(spawntable$SPAWNTIMEBETW[1],digits=1)` (`r round(spawntable$STSD[1],digits=1)`)|`r round(spawntable$CLUSTPRIOR[1],digits=1)` (`r round(spawntable$CPSD[1],digits=1)`)|`r round(spawntable$SPAWNAREA[1],digits=1)` (`r round(spawntable$SASD[1],digits=1)`)|
|All Female  | 5  |`r round(spawntable$MEANCLUST[2],digits=1)` (`r round(spawntable$MCSD[2],digits=1)`)|`r round(spawntable$TIMEMEAN[2],digits=1)` (`r round(spawntable$TIMESD[2],digits=1)`)|`r round(spawntable$SITEVISIT[2],digits=1)` (`r round(spawntable$SVSD[2],digits=1)`)|`r round(spawntable$DISTBETW[2],digits=2)` (`r round(spawntable$DSD[2],digits=2)`)|`r round(spawntable$TIMEBETW[2],digits=1)` (`r round(spawntable$TSD[2],digits=1)`)|-|`r round(spawntable$AREA[2],digits=1)` (`r round(spawntable$ASD[2],digits=1)`)|
|Low Activity Female| 5 |`r round(spawntable$NOMEANCLUST[2],digits=1)` (`r round(spawntable$NCSD[2],digits=1)`)|`r round(spawntable$NOTIMEMEAN[2],digits=1)` (`r round(spawntable$NOTIMESD[2],digits=1)`)|`r round(spawntable$NOSITEVISIT[2],digits=1)` (`r round(spawntable$NSSD[2],digits=1)`)|`r round(spawntable$NODISTBETW[2],digits=2)` (`r round(spawntable$NDSD[2],digits=2)`)|`r round(spawntable$NOTIMEBETW[2],digits=1)` (`r round(spawntable$NTSD[2],digits=1)`)|-|`r round(spawntable$NOAREA[2],digits=1)` (`r round(spawntable$NASD[2],digits=1)`)|
|High Activity Female| 5 |`r round(spawntable$MEANSPAWN[2],digits=1)` (`r round(spawntable$MSSD[2],digits=1)`)|`r round(spawntable$SPAWNTIMEMEAN[2],digits=1)` (`r round(spawntable$SPAWNTIMESD[2],digits=1)`) | `r round(spawntable$SITESPAWN[2],digits=1)` (`r round(spawntable$SSSD[2],digits=1)`)|`r round(spawntable$SPAWNDISTBETW[2],digits=2)` (`r round(spawntable$SDSD[2],digits=2)`)|`r round(spawntable$SPAWNTIMEBETW[2],digits=1)` (`r round(spawntable$STSD[2],digits=1)`)|`r round(spawntable$CLUSTPRIOR[2],digits=1)` (`r round(spawntable$CPSD[2],digits=1)`)|`r round(spawntable$SPAWNAREA[2],digits=1)` (`r round(spawntable$SASD[2],digits=1)`)|

##Discussion  

I have provided a detailed description of the lake trout mating system at the whole lake scale in a typical northern boreal lake. Using an acoustic telemetry monitoring system and a novel spatial temporal clustering analysis, I was able to quantify lake trout spawning movements and behaviours over the course of an entire spawning season. Lake trout were found to cluster on spawning shoals virtually around the entire nearshore region of Alexie Lake, as well as around several islands, which appears to further confirm previous findings that subtable spawning habitat is abundant in Alexie Lake (Chapter 2). Consistent with other studies, males arrived earlier than females and spent longer durations on spawning shoals over the course of the spawning season [@miller_observations_1948; @royce_breeding_1951; @martin_lake_1980; @muir_lake_2012]. Males formed >4 times as many spawing clusters and visited more sites than females. Spawning clusters predominantly were formed at night but were also observed during daylight hours, especially during the peak spawning season. Although daily travel distances were similar between sexes, higher activity rates and longer periods spent on spawning shoals by males suggests that males may exert more energy than females during the spawning season. Overall, females performed more linear movements over the course of the spawning season suggesting a searching behaviour, while males were less persitent and more random in there movements. 

Suitable spawning habitat in Alexie Lake was found to be abundant and widespread (Chapter 2), so it was not surprising that spawning clusters were also found dispersed throughout the entire nearshore region of Alexie Lake. Individual lake trout utlized a wide range of spawning areas (0.01--1.9 km^2^), with males visiting more than twice as many sites. The area of spawning site use was not found to differ between males and females in spite of increased male spawning site visits. Further, no individual lake trout visited all spawning sites or formed spawning clusters over the full extent of the collective lake trout spawning area (4.2 km^2^), suggesting indidvidual lake trout only select subset of all possible spawning shoals over the duration of the spawning season. Whether this is a function of territories that individuals stay within [although no evidence has been found to date in other lake trout populations; @muir_lake_2012] or, more likely, the result of lake trout selecting spawning shoals with preferred environmental conditions [i.e. wind events; @martin_lake_1980; @muir_lake_2012] or mates [@binder_new_2015] is unknown and requires further investigation. In chapter 2, I found that high variability in wind direction within and across spawning seasons suggests wind, and the resulting wave induced currents, are not predictable during the lake trout spawning season. This unpredictability in wind and the abundance of suitable habitat may favour spawning across multiple sites by lake trout as part of a bet-hedging strategy [Chapter 2; @fitzsimons_relationship_2014]. Diversifying the spatial, physical, and chemical characteristics of spawning habitat, to increase the overall portfolio performance of reproductive success [@moore_synchronization_2010] may be a productive tactic to buffer the unpredictable nature of wind and weather conditions.

Consistent with other studies, males formed spawning clusters earlier than females and spent longer durations on spawning shoals over the course of the spawning season [@miller_observations_1948; @royce_breeding_1951; @martin_lake_1980; @muir_lake_2012]. Males formed four times as many spawing clusters and visited twice as many sites than females. These data support the hypothesis that males maximize reproductive fitness by spending as much time as possible on spawning shoals, thus maximizing possible mate encounters and fertilizations [@muir_lake_2012]. Another explanation for the early arrival of males on spawning grounds is for signalling and attracting females [@gunn_spawning_1995; @esteve_lake_2008; @muir_lake_2012]. Recently, It has been suggested that males use display courtship behaviour such as "finning" in Great Bear Lake [@muir_lake_2012] or tactile courtship behaviour "hovering" in Lake Huron [@binder_new_2015]. Both behaviours involve small groups of males performing a display behaviour to attract females with the commonality of both occuring at the outer edge of the spawning shoal. This suggests part of the male reproductive strategy is to not only remain on active spawning shoals but also attract females to the spawning shoals it occupies, but how exactly this is accomplished remains largley unknown and requires further study. 
 
Lake trout migration onto spawning shoals coincided with surface temperatures of ~15&deg;C. Thermal preference of lake trout is between 5&deg;--15&deg;C [@plumb_performance_2009], therefore migration onto spawning shoals and the onset of spawning appears to be strongly regulated by thermal access to spawning habitat. In 2013, spawning onset began on August 29 and lasted over a 55 d period, concluding on October 23 when surface water temperature had declined to 5.5&deg;C. The length of the 2013 spawning period in Alexie Lake is quite long in comparison to other systems, where the spawning season typically lasts <14 d [@martin_lake_1980], although peak spawn duration (10 d) fell within this range. Lake trout predominately spawned at night, typically around 23:00, but spawning also occured frequently throughout the day especially during the peak season. This finding is consisitent with an increasing number of studies showing that lake trout do not only spawn at night and that day spawning does occur and may not be as rare as once believed [@esteve_lake_2008; @muir_lake_2012; @binder_new_2015]. 

Female lake trout invest large amounts of energy into their eggs, thus selecting high quality habitat and/or mates should be a priority. In contrast, males invest very little energy in their gametes, thus energy should be put towards maximizing mating events by finding and remaing on active spawning shoals [@muir_lake_2012]. I found supporting evidence for this with males expending greater energy on average than females over the course of the spawning season. This higher energy expenditure for males is in spite of travelling similar daily distances as females. The discrepacy in energy output is likely due to males increased frequency and duration on spawning shoals. Thus it appears that the disparity in gamete investment, where males invest less energy in gamete production, results in males expending more energy in attaining mates and egg fertilizations. 

Females displayed greater linearaity in movements than males over the course of the spawning season which appears to be more energeticly efficient than males non-linear movemements. I interpret the increased linearity as the result of a female  search strategy for mates and/or spawning habitat [@laidre_females_2012]. Females also visited less sites than males which could indicate that females have stricter requirements for spawning site selection than males; or, simply a byproduct of shorter duration on spawning shoals.

Although agonistic behaviour has not been observed in lake trout, evolutionary theory predicts male-male competition for mates shoud occur, especially when sex-ratios are skewed towards males [@clutton-brock_potential_1992] as is the case in lake trout [@martin_lake_1980; @muir_lake_2012; @binder_new_2015]. During spawning events, several males will travel with a female [@esteve_lake_2008; personal observation]. Typically, males will travel on one side of a female while jockeying for position next to her, even when the other side of the female remains unoccupied by other males [@binder_new_2015]. Binder et al. (2015) suggests that jockeying is likely a form of competition where the most fit individual should remain in closest proximity to the female and that closer proximity should result in an increase in fertilizations [@blanchfield_breeding_2003;@mjolnerod_mate_1998 ]. This form of endurance competition is also likely the reason for increased energy expenditure over females, although both sexes participate in the same spawning behaviours [i.e. travelling, sinking, and gamete release; @esteve_lake_200], males remain on spawning shoals for a longer duration than females and likely perform greater frequencies of these behaviours.

The deviation from the general salmonine nest building strategy to an itinerant strategy may be a reflection of the predictability of habitat quality. Mating system theory would suggest that when there is no parental care, female reproductive success is driven by habitat and mate quality, whereas male reproductive success is driven by spawning frequency [@degaudemar_sexual_1998]. Habitat quality in salminine systems is driven by water currents [@chapman_critical_1988], with females fighting for territories with high quality redd sites and males competing for proximity to females [@esteve_observations_2005]. But in lakes in the absence of predictable current flow (i.e. river or  groundwater inflow), wave induced currents are belived to drive habitat quality, although little quantitiative evidence supports this [@martin_lake_1980; @marsden_recognition_1995] and may be due to winds unpredictability within and among spawning seasons [Chapter 2]. Thus, investing energy in building a nest on a site with unpredictable habitat quality may have been abondoned in favour of moveing among several spawning sites with suitable substrate but unknown and varying habitat quality. 



#References
