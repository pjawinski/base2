* Encoding: UTF-8.

* =====================
* === Get SPSS data ===
* =====================.

* Set working directory.
cd '/Users/philippe/Desktop/base2/'.

* Get MASTER data.
GET FILE='data/BASE_II_for_Markett_jd.sav'.
DATASET NAME MASTER.
DATASET ACTIVATE MASTER.

* Get CESDetc.
GET FILE='data/BASE_II_CFGG_red+MMSE-GDS-CESD-Labor_jd.sav'.
DATASET NAME CESDetc.
DATASET ACTIVATE CESDetc.

* Get AGE.
GET FILE='data/BASE_II_MRT_age_jd.sav'.
DATASET NAME AGE.
DATASET ACTIVATE AGE.

* Get Income.
GET FILE='data/HHinc_BASE_II.sav'.
DATASET NAME Income.
DATASET ACTIVATE Income.
RENAME VARIABLE (id= ID).
EXECUTE.
SELECT IF (NOT(ID EQ "")).
EXECUTE.


* Merge Master and CESDetc.
DATASET ACTIVATE MASTER.
SORT CASES BY ID.
DATASET ACTIVATE CESDetc.
SORT CASES BY ID.
ALTER TYPE ID (A30).
DATASET ACTIVATE MASTER.
MATCH FILES /FILE=*
  /FILE='CESDetc'
  /BY ID.
EXECUTE.

* Merge Master and AGE.
DATASET ACTIVATE MASTER.
SORT CASES BY ID.
DATASET ACTIVATE AGE.
SORT CASES BY ID.
ALTER TYPE ID (A30).
DATASET ACTIVATE MASTER.
MATCH FILES /FILE=*
  /FILE='AGE'
  /BY ID.
EXECUTE.

* Merge Master and Income.
DATASET ACTIVATE MASTER.
SORT CASES BY ID.
DATASET ACTIVATE Income.
SORT CASES BY ID.
ALTER TYPE ID (A30).
DATASET ACTIVATE MASTER.
MATCH FILES /FILE=*
  /FILE='Income'
  /BY ID.
EXECUTE.

* Close datasets.
DATASET CLOSE CESDetc.
DATASET CLOSE AGE.
DATASET CLOSE Income.

* ===========================
* === Calculate variables ===
* ===========================.

* set missings.
missing values MMSE_Jahr to MMSE_Summe (998 thru hi).
missing values CESD_01 to CESD_20 (998 thru hi).
missing values GDS_01 to GDS_15 (998 thru hi).
missing values Alkohol_Menge(998 thru hi).
missing values Alkohol_6Glaser(998 thru hi).

* RECODE Rauchen_aktuell INTO Rauchen_aktuell_inverted.
RECODE Rauchen_aktuell  (4 = 0) (3 = 1) (2 = 2) (1 = 3) (MISSING=SYSMIS) (SYSMIS=SYSMIS) INTO Rauchen_aktuell_inverted.
EXECUTE.
* Define Variable Properties.
*Rauchen_aktuell_inverted.
VARIABLE LEVEL  Rauchen_aktuell_inverted(ORDINAL).
VARIABLE LABELS  Rauchen_aktuell_inverted 'Rauchen Sie zur Zeit?'.
VALUE LABELS Rauchen_aktuell_inverted
  1.00 'nein, noch nie '
  2.00 'nein, seit >1 Jahr nicht mehr'
  3.00 'nein, seit < 1 Jahr nicht mehr'
  4.00 'ja'.
EXECUTE.
FORMATS Rauchen_aktuell_inverted(f8.0).

* RECODE Educ_final.
RECODE Educ_final  (-9 thru -1 = SYSMIS) (99 = SYSMIS) (MISSING=SYSMIS) (SYSMIS=SYSMIS) (ELSE = COPY).
EXECUTE.
MISSING VALUES Educ_final().

* RECODE MMSE_SUMME.
RECODE MMSE_SUMME  (MISSING=SYSMIS) (SYSMIS=SYSMIS) (ELSE = COPY).
EXECUTE.
MISSING VALUES MMSE_SUMME().

* RECODE pnett.
RECODE pnett (-3 thru -1 = SYSMIS) (MISSING=SYSMIS) (SYSMIS=SYSMIS) (ELSE = COPY).
EXECUTE.

* RECODE hhnetto.
RECODE hnetto (-3 thru -1 = SYSMIS) (MISSING=SYSMIS) (SYSMIS=SYSMIS) (ELSE = COPY).
EXECUTE.
MISSING VALUES Alkohol_6Glaser().

* RECODE Alkohol_haeufigkeit.
VALUE LABELS Alkohol_haufigkeit
  0.00 'nie'
  1.00 'einmal im Monat oder seltender'
  2.00 'zwei- bis viermal im Monat'
  3.00 'zwei- bis dreimal die Woche'
  4.00 'Viermal die Woche oder öfter'.
EXECUTE.

* RECODE Alkohol_Menge.
RECODE Alkohol_Menge (998 thru high = SYSMIS) (MISSING=SYSMIS) (SYSMIS=SYSMIS) (ELSE = COPY).
EXECUTE.
MISSING VALUES Alkohol_Menge().
VALUE LABELS Alkohol_Menge
  0.00 '1 bis 2 Gläser pro Tag'
  1.00 '3 bis 4 Gläser pro Tag'
  2.00 '5 bis 6 Gläser pro Tag'
  3.00 '7 bis 9 Gläser pro Tag'
  4.00 '10 oder mehr Gläser pro Tag'.
EXECUTE.

* RECODE Alkohol_6Glaser.
RECODE Alkohol_6Glaser (998 thru high = SYSMIS) (MISSING=SYSMIS) (SYSMIS=SYSMIS) (ELSE = COPY).
EXECUTE.
MISSING VALUES Alkohol_6Glaser().
VALUE LABELS Alkohol_6Glaser
  0.00 'nie'
  1.00 'seltener als einmal im Monat'
  2.00 'einmal im Monat'
  3.00 'einmal in der Woche'
  4.00 'täglich oder fast täglich'.
EXECUTE.

* Calculate GDS_score.
RECODE GDS_01 GDS_05 GDS_07 GDS_11 GDS_13 (1=0) (2=1) (MISSING=SYSMIS) (999=SYSMIS) (998=SYSMIS) 
    INTO GDS_01u GDS_05u GDS_07u GDS_11u GDS_13u.
EXECUTE.

RECODE GDS_02 GDS_03 GDS_04 GDS_06 GDS_08 GDS_09 GDS_10 GDS_12 GDS_14 GDS_15 (1=1) (2=0) (MISSING=SYSMIS) (999=SYSMIS) (998=SYSMIS) 
    INTO GDS_02u GDS_03u GDS_04u GDS_06u GDS_08u GDS_09u GDS_10u GDS_12u GDS_14u GDS_15u.
EXECUTE.

COUNT GDS_missing=GDS_01u GDS_02u GDS_03u GDS_04u GDS_05u GDS_06u GDS_07u GDS_08u GDS_09u GDS_10u 
     GDS_11u GDS_12u GDS_13u GDS_14u GDS_15u (SYSMIS).
EXECUTE.

IF  (GDS_missing < 2) GDS_Summe = MEAN(GDS_01u,GDS_02u,GDS_03u,GDS_04u,GDS_05u,GDS_06u,GDS_07u,GDS_08u,GDS_09u,GDS_10u,
      GDS_11u,GDS_12u,GDS_13u,GDS_14u,GDS_15u)*15.
EXECUTE.

* Calculate CESD_score.
COUNT CESD_missing=CESD_01 CESD_02 CESD_03 CESD_04 CESD_05 CESD_06 CESD_07 CESD_08 CESD_09 CESD_10
 CESD_11 CESD_12 CESD_13 CESD_14 CESD_15 CESD_16 CESD_17 CESD_18 CESD_19 CESD_20 (SYSMIS,MISSING,999,998).
EXECUTE.

RECODE CESD_01 CESD_02 CESD_03 CESD_05 CESD_06 CESD_07 CESD_09 CESD_10
 CESD_11 CESD_13 CESD_14 CESD_15 CESD_17 CESD_18 CESD_19 CESD_20 (1=0) (2=1) (3=2) (4=3) (MISSING=SYSMIS) (999=SYSMIS) (998=SYSMIS) 
    INTO CESD_01u CESD_02u CESD_03u CESD_05u CESD_06u CESD_07u CESD_09u CESD_10u
 CESD_11u CESD_13u CESD_14u CESD_15u CESD_17u CESD_18u CESD_19u CESD_20u.
EXECUTE.

RECODE CESD_04 CESD_08 CESD_12 CESD_16 (1=3) (2=2) (3=1) (4=0) (MISSING=SYSMIS) (999=SYSMIS) (998=SYSMIS) 
    INTO CESD_04u CESD_08u CESD_12u CESD_16u.
EXECUTE.

IF  (CESD_missing < 3) CESD_Summe = MEAN(CESD_01u,CESD_02u,CESD_03u,CESD_04u,CESD_05u,CESD_06u,CESD_07u,CESD_08u,CESD_09u,CESD_10u,
 CESD_11u,CESD_12u,CESD_13u,CESD_14u,CESD_15u,CESD_16u,CESD_17u,CESD_18u,CESD_19u,CESD_20u)*20.
EXECUTE.

* rename uric acid.
RENAME VARIABLE (HarnsäuremgdL = HarnsaeuremgdL ).
EXECUTE.

* ================================================================
* === Data extraction for statistical analysis in Matlab and R ===
* ================================================================.

* select if age at MRI date is available.
SELECT IF (NOT(SYSMIS(age_T1_LIP_MRT))).
EXECUTE.

** Save dataset with variables of interest in .sav format.
SAVE OUTFILE='data/01_phenotypes_master.sav'
   /KEEP=ID sex age_T1_LIP_MRT Educ_final hnetto MMSE_Summe GDS_Summe CESD_Summe Rauchen_aktuell_inverted
Alkohol_haufigkeit Alkohol_Menge Alkohol_6Glaser BMI RRdi RRsy finalMetLscore CH_Diabetes HOMAIR HbA1c BZP1 BZP2 GammaGTGGTUL HarnsaeuremgdL TNF1 DS2_corr EM_final WM_final Gf_final futi_mean cfc_mean.

** Save dataset for analyses in Matlab and R in tab-delimited .txt format.
GET FILE='data/01_phenotypes_master.sav'.
DATASET NAME phenotypes.
DATASET ACTIVATE phenotypes.
DATASET CLOSE MASTER.
EXECUTE.

** check for missings and unusual values.
FREQUENCIES VARIABLES=sex age_T1_LIP_MRT Educ_final hnetto MMSE_Summe GDS_Summe CESD_Summe 
    Rauchen_aktuell_inverted Alkohol_haufigkeit Alkohol_Menge Alkohol_6Glaser BMI RRdi RRsy finalMetLscore CH_Diabetes 
    HOMAIR HbA1c BZP1 BZP2 GammaGTGGTUL HarnsaeuremgdL TNF1 DS2_corr EM_final WM_final 
    Gf_final futi_mean cfc_mean
  /FORMAT=NOTABLE
  /NTILES=4
  /STATISTICS=STDDEV MINIMUM MAXIMUM MEAN MEDIAN SKEWNESS SESKEW KURTOSIS SEKURT
  /ORDER=ANALYSIS.

SAVE TRANSLATE OUT = 'data/01_phenotypes_master.txt'
   / TYPE=TAB
   / FIELDNAMES
   / REPLACE
   / MAP. 
