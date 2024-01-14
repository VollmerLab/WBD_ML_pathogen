Electronic Sea Surface Temperature has been recording using several different sensors. Currently the Campbell Scientific Model 107 temperature sensor is being used.

Data represent the summary of the period ending at the time stamp. 

Currently, measurements are made every 10 seconds and recorded every 15 minutes.
Measurements were originally made hourly and changed to every 15 minutes on June 2, 2005

File Structure:
Datetime (dd/mm/yyyy HH:MM:SS), Date (dd/mm/yyyy), Curated Values (C), Raw Values (C), CHK_NOTE, CHK_FAIL

CHK_NOTE and CHK_FAIL provide Quality Assurance information. 

CHK_NOTE indicates quality of record. Possible values are:
	"adjusted"		-> Value of datum adjusted or gap filled
	"bad"			-> Datum failed one or more test and should not be used
	"doubtful"		-> Datum suspect and should be used with caution
	"good"			-> Datum acceptible
	"missing"		-> Datum missing
	"nc"	  	    -> Datum not yet checked

CHK_FAIL indicates the reason that a datum failed or was adjusted. Possible values are:
	"Calibration"	-> Sensor calibration corrected
	"Persistence"	-> Variation of data below expected range
	"Range"			-> Datum outside of 99.9 percentile
	"Spike"			-> Short-term increase or decrease of data significantly outside normal variance range
	"Gap Fill"		-> Data estimated or supplied from another source
	

Quality Assurance Procedures:
	- Data are loaded into a custom-built program which premits easy visual examination of the data (including simultaneous viewing of rainfall data)
	- A gap-filling procedure is carried out to automatically fill small (<= 3 consecutive records) gaps using simple linear regression
	- Larger gaps (<24 hours) are filled manually by drawing a curve that best approximates the missing data based on available, nearby data, precipitation during the gap period, and
		experience with the data by the analyst. The size of the largest gaps that are filled depends on the variability of the variable and its predictability.
		This is obviously subjective, but is believed to provide very reasonable results.
	- Data are visually examined for anomalies and are flagged, and adjusted when possible (ie. when confidence is high that the adjustment provides good results).
	
Missing Data
	Missing data are indicated by the value of -999.
	If a missing datum has been gap filled then the value of the variable RAW will be -999, but the curated variable (eg. at, rh, sr, etc) will have the estimated value. 
		CHK_NOTE will have the value "adjusted" and CHK_FAIL will have the value "Gap Fill"
	If a missing datum has not been gap filled, both RAW and the Curated variable will have the value -999
		CHK_NOTE will have the value "missing" and CHK_FAIL will be empty. 
