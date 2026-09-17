GENERAL INFORMATION
1. Title of Dataset: Multi-model Ensemble for Robust Verification of hydrological modeling in Japan (MERV-Jp) 

2. Author Information:
   Principal investigator 
      Name: Yohei Sawada (yoheisawada@g.ecc.u-tokyo.ac.jp)
               Associate Professor, Department of Civil Engineering, The University of Tokyo

   Co-investigator
      Name: Shinichi Okugawa (okugawa@g.ecc.u-tokyo.ac.jp)
               Technical staff, Department of Civil Engineering, The University of Tokyo

3. Overview:
	MERV-Jp is the dataset of meteorological forcing and multi-model runoff simulation in 135 Japanese basins (ver1.1) / 87 Japanese basins (ver2.1), and contributes to carrying out a large sample rainfall-runoff simulation in Japan. In addition, MERV-Jp can be used as a benchmark to evaluate user's hydrological modeling. 
	The detailed description of MERV-Jp can be found at "Y. Sawada, S. Okugawa and T. Kimizuka (2022): Multi-model ensemble benchmark data for hydrological modeling in Japanese river basins, Hydrological Research Letters, 16, 73-79"  (https://doi.org/10.3178/hrl.16.73).

4. Precipitation data correction
  Precipitation data for the period from 2012 to 2015 of ver2_0 was incorrect and is corrected in ver2_1

DATA RESOURCE & METHODOLOGY
1. Summary
  To simulate runoff in 135 Japanese basins (ver1.1) / 87 Japanese basins (ver2.1), the precipitation data, temparature data, potential evapotranspiration data and observed runoff data during 11 years (1993 - 2003 : ver1.1) / 25~30 years (1986 - 2015 : ver2.1) are collected and processed. 
  Then, rainfall-runoff simulation was performed by 44 hydrological models in Modular Assessment of Rainfall-Runoff Models toolbox: Knoben et al., 2019

2. DATA resource
 (1) Basin information
   The geographical data of basins were retrieved from Global Runoff Data Centre (https://www.bafg.de/GRDC/EN/Home/homepage_node.html) [Date of retrieval: 2021-07-12]. 

 (2) Observed Runoff data
   [ver1.1]
   The observed Runoff data were retrieved from Global Runoff Data Centre (https://www.bafg.de/GRDC/EN/Home/homepage_node.html) [Date of retrieval: 2021-07-12]. 
   [ver2.1]
   The observed Runoff data were retrieved from Water Information System of Ministry of Land, Infrastructure, Transport and Tourism (http://www1.river.go.jp/) [Date of retrieval: 2023-02-15 ~ 2023-04-12] and 
   Global Runoff Data Centre (https://www.bafg.de/GRDC/EN/Home/homepage_node.html) [Date of retrieval: 2021-07-12]. 

 (3) Precipitation data
   The average precipitation data of each basin were calculated from APHRODITE data (http://aphrodite.st.hirosaki-u.ac.jp/) 
   The source of precipitation data from 1986 to 2011 were got from APHRO_JP V1207 (0.05mesh)
   The source of precipitation data from 2012 to 2015 were got from APHRO_JP V1207_R3 (0.05mesh)

 (4) Temparature data
   The average temparature data of each basin were calculated from APHRODITE data (http://aphrodite.st.hirosaki-u.ac.jp/). 
   The source of temparature data were got from APHRO_MA_TAVE_025deg_V1808 (0.25mesh)

 (5) Potential evapotranspiration (PET) data
   The average PET data of each basin were calculated from hPET data of the University of Bristol (https://data.bris.ac.uk/data/dataset/qb8ujazzda0s2aykkv0oq0ctp).


3. Methodological Information
   The meteorological forcing data shown above were processed and used as input data for hydrological models.
   We drove 44 hydrological models. The model parameters were calibrated by the first 8 years data (ver1.1) / the first 5 years data (ver2.1) of observed runoff data. The simulated runoff for 11 years (ver1.1) / 25~30 years (ver2.1) by the calibrated models was provided as MERV-Jp.
   The processed meteorological forcing and observed runoff data were also provided for each basin in MERV-Jp.
   

DATASET FILE OVERVIEW
1. Summary
  The main data of MERV-Jp contain 135 csv files (ver1.1) / 87 csv files (ver2.1) in which meteorological forcing data and observed/simulated runoff data are included.
  In addition, we have 4 additional files (including this readme) as the supporting information. 

2. File Specific Information
 (1) Meteorological forcing and observed & simulated runoff data
   Number of file: 135 (ver1.1) / 87 (ver2.1)
   Zip file name: varssim.zip   => includes below directory
      varssim
        |- ver1_1   => includes 135 csv files(*) of ver1.1
        |- ver2_1   => includes 87 csv files(*) of ver2.1
      *csv file name: 'varssim'+(basin ID)+'.csv'

   Data contents of each file:
     Daily data during 1993 to 2003 (ver1.1) / 25~30 years of 1986 to 2015 (ver2.1) 
     (i) Flag of period: '0' for calibration period and '1' for evaluation period
     (ii) Date (Year/Month/Day)
     (iii) Basin-averaged Precipitation (mm/day)
     (iv) Basin-averaged Temperature (degree)
     (v) Basin-averaged PET (mm/day)
     (vi) Observed runoff rate (mm/day)
     (vii~) Simulated runoff rate (mm/day) of 44 models

 (2) Basin list
   File name: 'basinlist_all.xlsx'
   Data contents: 
     (i) Basin index number of ver1.1
     (ii) Basin index number of ver2.1
     (iii) GRDC station number
     (iv) Name of river
     (v) Name of river (Japanese Kanji)
     (vi) Downstream GRDC station number
           *If there is a downstream basin which includes its basin, the station number of downstream basin is provided. If not, this number is '0'. 
     (vii) Basin area (km^2)
     (viii) Name of station
     (ix) Name of station (Japanese Kanji)
     (x) Total years of observed runoff data during 1986 to 2015 (Total number of xi ~ xl)
     (xi~) Flag indicating the presence of observed runoff data in that year
           *"0": None of daily observed runoff data in that year, "1": At least one daily observed runoff data in that year. 

 (3) Basin map of ver1.1
   File name: 'basinmap_135.png'
   Data contents: 
     (i) Geographical map of each basin with Basin index number

 (4) Basin map of ver2.1
   File name: 'basinmap_87.png'
   Data contents: 
     (i) Geographical map of each basin with Basin index number

 (5) Readme
   File name: 'readme_21.txt'
   This file
