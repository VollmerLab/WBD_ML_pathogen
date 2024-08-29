Supplementary Files
SupplementaryFile1.csv.gz – Metadata for field collected samples with columns:
	•	“sample_id” – individual sample names.
	•	“health” – “H” healthy and “D” diseased fragments.
	•	“year” – year fragment collected.
	•	“season” – season fragment collected (“S” July and “W” January)
	•	“site” – location fragment collected from
	•	“lib.size” – total number of sequenced reads
	•	“norm.factors” – factor used to normalize read counts of ASVs

SupplementaryFile2.csv.gz – Metadata for tank collected samples with columns: 
	•	“sample_id” – individual sample names.
	•	“geno” – fragment genotype
	•	“fragment_id” – fragment identification tracked through repeated sampling
	•	“tank_id” – tank identification
	•	“time_treat” – concatenated metric for sampling time, exposure, and disease outcome separated by “_”
		o	Time – 0, 2, 8
		o	Exposure – “D” Diseased, “N” Healthy
		o	Disease Outcome - “D” Diseased, “H” Healthy
	•	“lib.size” – total number of sequenced reads
	•	“norm.factors” – factor used to normalize read counts of ASVs

SupplementaryFile3.fasta – FASTA file including complete 16s sequences named with ASV identifier and taxonomy.

SupplementaryFile4.csv.gz – Matrix of the number of reads of each ASV sequenced in each sample. Combined both field and tank samples.

SupplementaryFile5.csv.gz – Matrix of the log2 CPM of each ASV sequenced in each sample. Combined both field and tank samples.

SupplementaryFile6.csv.gz – Complete results for each ASV association. 
	•	“top_classification” – lowest taxonomic classification with more than 80% confidence.
	•	“taxonomy” – Full taxonomy including confidence in each taxonomic level.
	•	“passedFilter” – indicates taxa filtered from analysis due to rarity and/or lack of observations across sample times. 
	•	“rank_*” – machine learning model rankings, median ranking, and model estimated ranking along with standard error, confidence interval, and FDR adjusted p-value used to identify important ASVs.
	•	“ml_retained” – Indicates if the ASV was of above average importance to ML models. NA values indicate ASVs which were filtered prior to ML modelling.
	•	“fieldModel_*” – ANOVA table results for each ASV testing the effects of health, year, season and all possible interactions indicating:
		o	Sums of squares, mean squares, numerator and denominator degrees of freedom, F statistic, p-value, and FDR corrected p-value. 
		o	NA values are filled for ASVs filtered prior to differential abundance analysis.
	•	“diffAbundance_healthAssociation” – Marks the health association of ASVs from differential abundance analysis of field samples: “H” health, “D” diseased, “N” none, NA – filtered prior to differential abundance analysis.
	•	“fieldLogFC_*” – Post-hoc contrasts for ML retained ASVs testing the significance of the log2 fold-change between disease and healthy fragments within each sampling time (year: 2016, 2017 & season: “S” July, “W” January) 		showing:	
		o	Mean estimate, standard error, degrees of freedom, lower and upper 95% confidence interval, t-statistic, p-value, FDR adjusted p-value.
		o	NA values are filled for ASVs which were not marked as important by ML models.
	•	“field_consistent” – Indicates if the ASV was consistently healthy or disease associated across sampling times. NA values are filled for ASVs which were not marked as important by ML models.
	•	“tankModel_*” – ANOVA table results for each ASV testing the effect of the combination of time, disease exposure, and disease outcome, indicating: 
		o	Sums of squares, mean squares, numerator and denominator degrees of freedom, F statistic, p-value, and FDR corrected p-value. 
		o	NA values are filled for ASVs filtered prior to tank experimental analysis.
	•	“tankLogFC_*” – Post-hoc contrasts for ASVs tested in tank exposure experiments. 
		o	Contrasts include: 
		Post-exposure diseased vs healthy outcome regardless of exposure (DvH)
		Post-exposure diseased vs healthy exposure regardless of outcome (DvN)
		Post-exposure disease exposed corals with disease symptoms compared to disease exposed but still healthy corals (DDvDH)
		Post-exposure disease exposed corals with disease symptoms compared to healthy exposed and still healthy corals (DDvNH)
		Post-exposure disease exposed corals which stay healthy compared to healthy exposed and still healthy corals (DHvNH)
		Pre-exposure compared to Post-exposure in corals with the disease regardless of exposure (PostvPreD)
		Pre-exposure compared to Post-exposure in corals without the disease regardless of exposure (PostvPreH)
		o	Mean estimate, standard error, degrees of freedom, lower and upper 95% confidence interval, t-statistic, p-value, FDR adjusted p-value.
		o	NA values are filled for ASVs which were not consistently associated with healthy or diseased corals in the field experiment.
	•	“pathogen_classification” – Indicates the predicted microbial classification based on the tank results. Pathogen, Opportunist, Commensal
		o	NA values are filled for ASVs which were not consistently associated with healthy or diseased corals in the field experiment.
