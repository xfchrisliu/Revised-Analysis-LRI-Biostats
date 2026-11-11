options ls=200 ps=60 nodate nonumber msglevel=I validvarname=any 
	validmemname=extend;
libname xf "/home/liux3/Manuscript/26/Revised Manuscript";
****** Read excel sheets into SAS;
libname myexcel xlsx 
	"/home/liux3/Manuscript/26/Revised Manuscript/Data/Final new data.xlsx";

/*proc contents data=myexcel._all_;
run;*/

********** Macro pipeline to conditionally execute the correct GLM ***************;

%macro GLM_GermDiet;
	%put Note: Normality test p value is &pvalue_norm;
	*ods trace on;
	ods select none;
	ods output ModelANOVA=effect OverallANOVA=overall LSMeans=germdietdiff 
		SlicedANOVA=germdiff;

	%if %sysevalf(&pvalue_norm >=0.05, BOOLEAN) %then
		%do;
			%put NOTE: Data is normally distributed. Running Normal model.;

			proc glm data=fig_1;
				class 'germ status'n(ref="MPF") diet(ref="Palm");
				* Your categorical independent variable;
				model &var='germ status'n diet 'germ status'n * diet;
				lsmeans 'germ status'n*diet /slice=diet pdiff=control('MPF' 'Palm') 
					adjust=dunnett;
				run;
			quit;

			ods output LSMeans=dietdiff;

			proc sort data=fig_1;
				by 'germ status'n;
			run;

			proc glm data=fig_1;
				by 'germ status'n;
				class diet(ref="Palm");
				model &var=diet;
				lsmeans diet /pdiff=control('Palm') adjust=dunnett;
				run;
			quit;

		%end;
	%else
		%do;
			%put NOTE: Data is non-normal. Running log normal model.;

			data fig_1;
				set fig_1;
				*&var=log(&var+0.00000001);
				&var=log10(&var+0.00000001);
			run;

			proc glm data=fig_1;
				class 'germ status'n(ref="MPF") diet(ref="Palm");
				* Your categorical independent variable;
				model &var='germ status'n diet 'germ status'n * diet;
				lsmeans 'germ status'n*diet /slice=diet pdiff=control('MPF' 'Palm') 
					adjust=dunnett;
				run;
			quit;

			ods output LSMeans=dietdiff;

			proc sort data=fig_1;
				by 'germ status'n;
			run;

			proc glm data=fig_1;
				by 'germ status'n;
				class diet(ref="Palm");
				model &var=diet;
				lsmeans diet /pdiff=control('Palm') adjust=dunnett;
				run;
			quit;

		%end;
	*ods trace off;
	ods select all;
%mend GLM_GermDiet;

%macro GLM_Diet;
	%put Note: Normality test p value is &pvalue_norm;
	*ods trace on;
	ods select none;
	ods output ModelANOVA=effect OverallANOVA=overall LSMeans=dietdiff;

	%if %sysevalf(&pvalue_norm >=0.05, BOOLEAN) %then
		%do;
			%put NOTE: Data is normally distributed. Running Normal model.;

			%if &clean_fig=%str(Fig 6 ) %then
				%do;

					proc glm data=fig_1;
						class diet(ref="Control");
						model &var=diet;
						lsmeans diet /pdiff=control('Control') adjust=dunnett;
						run;
					quit;

				%end;
			%else %if &clean_fig=Fig S1 %then
				%do;

					proc glm data=fig_1;
						class diet(ref="Palm");
						* Your categorical independent variable;
						model &var=diet;
						lsmeans diet /pdiff=control('Palm') adjust=dunnett;
						run;
					quit;

				%end;
		%end;
	%else
		%do;
			%put NOTE: Data is non-normal. Running log normal model.;

			data fig_1;
				set fig_1;
				*&var=log(&var+0.00000001);
				&var=log10(&var+0.00000001);
			run;

			%if &clean_fig=%str(Fig 6 ) %then
				%do;

					proc glm data=fig_1;
						class diet(ref="Control");
						model &var=diet;
						lsmeans diet /pdiff=control('Control') adjust=dunnett;
						run;
					quit;

				%end;
			%else %if &clean_fig=Fig S1 %then
				%do;

					proc glm data=fig_1;
						class diet(ref="Palm");
						model &var=diet;
						lsmeans diet /pdiff=control('Palm') adjust=dunnett;
						run;
					quit;

				%end;
		%end;
	*ods trace off;
	ods select all;
%mend GLM_Diet;



****** Micro funtion for effects p values with one variable in a figure;

%macro one_variable(fig, var);
	data fig_1;
		set myexcel.&fig;

		if diet=' ' then
			delete;
	run;

	options nonotes nosource nomprint;
	ods select none;
	ods graphics off;

	****** Test for normality and store the p-value in a macro variable;
	*** Run GLM with only fixed effects and output standardized Pearson residuals;
	ods select none;

	%if %sysfunc(substr(%sysfunc(dequote(&fig)), 1, 6))=%str(Fig 6 ) %then
		%do;

			proc glm data=fig_1;
				class diet(ref="Control");
				model &var=diet;
				*** Capture the standardized residuals here ***;
				output out=glm_norm_check r=raw_res stdr=std_res student=student_res;
			run;
			quit;

		%end;
	%else %if %sysfunc(substr(%sysfunc(dequote(&fig)), 1, 6))=Fig S1 %then
		%do;

			proc glm data=fig_1;
				class diet(ref="Palm");
				model &var=diet;
				*** Capture the standardized residuals here ***;
				output out=glm_norm_check r=raw_res stdr=std_res student=student_res;
				run;
			quit;

		%end;
	%else
		%do;

			proc glm data=fig_1;
				class 'germ status'n(ref="MPF") diet(ref="Palm");
				model &var='germ status'n diet 'germ status'n * diet;
				*** Capture the standardized residuals here ***;
				output out=glm_norm_check r=raw_res stdr=std_res student=student_res;
				run;
			quit;

		%end;
		
	**************** Check the fixed-effects residuals for normality ***************;
	*ods select TestsForNormality;
	*ods output TestsForNormality=NormalityResults;

	proc univariate data=glm_norm_check normal plots outtable=RawNormalityResults 
			noprint;
		var student_res;
		histogram student_res / normal;
		qqplot student_res / normal(mu=est sigma=est);
	run;

	data pvalue_normality;
		length fig_id $25;
		retain fig_id variable;
		set RawNormalityResults(rename=(_PROBN_=pvalue_SW));
		length variable $32;
		fig_id="%sysfunc(dequote(&fig))";
		variable="%sysfunc(dequote(&var))";
		keep fig_id variable pvalue_SW;
		format pvalue_SW best10.;
	run;

	data _null_;
		set RawNormalityResults;
		call symputx('pvalue_norm', _PROBN_);
		* Store a pvalue into &pvalue_norm;
	run;
	
	proc delete data=glm_norm_check RawNormalityResults;
	run;

	ods select all;
	
	*** Run macro pipeline for executing the correct GLM based on normality;
	%let clean_fig = %sysfunc(substr(%sysfunc(dequote(&fig)), 1, 6));

	%if &clean_fig=%str(Fig 6 ) or &clean_fig=Fig S1 %then
		%do;
			%GLM_Diet;
			
			****** Select p values of diet main effects and multiple comparisons **********;
			data effect;
				set effect;

				if HypothesisType=1 then
					delete;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc transpose data=effect out=wide_effect;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id source;
				* Uses values in this column to name the new columns;
				var probF;
				* The actual data values you want to move into columns;
			run;

			data pvalue_effect;
				retain fig_id variable;
				set wide_effect;
				length variable $32;
				drop _name_ _label_;
				variable="%sysfunc(dequote(&var))";
				*proc print;
			run;

			data overall;
				set overall;

				if source^="Model" then
					delete;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc transpose data=overall out=wide_overall;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id source;
				* Uses values in this column to name the new columns;
				var probF;
				* The actual data values you want to move into columns;
			run;

			data pvalue_overall(rename=(Model='OVERALL MODEL'n));
				retain fig_id variable;
				set wide_overall;
				length variable $32;
				drop _name_ _label_;
				variable="%sysfunc(dequote(&var))";
				*proc print;
			run;

			data pvalueFig6_FigS1_effect;
				length fig_id $25;
				merge pvalue_effect pvalue_overall;
				by fig_id;
				format diet 'OVERALL MODEL'n best10.;
				*proc print;
			run;

			****** Select p values of diet multiple comparisons ******;

			data dietdiff;
				set dietdiff;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc transpose data=dietdiff out=wide_dietdiff;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id diet;
				* Uses values in this column to name the new columns;
				var ProbtDiff;
				* The actual data values you want to move into columns;
			run;

			%if &clean_fig=%str(Fig 6 ) %then
				%do;

					data pvalueFig6_diet;
						length fig_id $25;
						retain fig_id variable;
						set wide_dietdiff;
						length variable $32;
						drop _name_ _label_ Control;
						variable="%sysfunc(dequote(&var))";
						format '12-HSA'n '12-Oxo'n RA best10.;
						*proc print;
					run;

				%end;
			%else %if &clean_fig=Fig S1 %then
				%do;

					data pvalueFigS1_diet;
						length fig_id $25;
						retain fig_id variable;
						set wide_dietdiff;
						length variable $32;
						drop _name_ _label_ Palm;
						variable="%sysfunc(dequote(&var))";
						format borage Echium Fish HLS HOS best10.;
						*proc print;
					run;

				%end;
				
				proc delete data=wide_effect wide_overall pvalue_overall wide_dietdiff;
				run;
		%end;
	%else
		%do;
			%GLM_GermDiet;
			
			****** Select p values of effects and multiple comparisons **********;

			data effect;
				set effect;

				if HypothesisType=1 then
					delete;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc transpose data=effect out=wide_effect;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id source;
				* Uses values in this column to name the new columns;
				var probF;
				* The actual data values you want to move into columns;
			run;

			data pvalue_effect;
				retain fig_id variable;
				set wide_effect;
				length variable $32;
				drop _name_ _label_;
				variable="%sysfunc(dequote(&var))";
				*proc print;
			run;

			data overall;
				set overall;

				if source^="Model" then
					delete;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc transpose data=overall out=wide_overall;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id source;
				* Uses values in this column to name the new columns;
				var probF;
				* The actual data values you want to move into columns;
			run;

			data pvalue_overall(rename=(Model='OVERALL MODEL'n));
				retain fig_id variable;
				set wide_overall;
				length variable $32;
				drop _name_ _label_;
				variable="%sysfunc(dequote(&var))";
				*proc print;
			run;

			data pvalue_effect;
				length fig_id $25;
				merge pvalue_effect pvalue_overall;
				by fig_id;
				format 'germ status'n diet 'germ status*diet'n 'OVERALL MODEL'n best10.;
				*proc print;
			run;

			****** Select p values of germ simple effects given each level of diet ******;

			data germdiff;
				set germdiff;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc transpose data=germdiff out=wide_germdiff;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id diet;
				* Uses values in this column to name the new columns;
				var probF;
				* The actual data values you want to move into columns;
			run;

			data pvalue_germ;
				length fig_id $25;
				retain fig_id variable;
				set wide_germdiff;
				length variable $32;
				drop _name_ _label_;
				variable="%sysfunc(dequote(&var))";
				format Borage Echium Fish HLS HOS Palm best10.;
				*proc print;
			run;

			****** Select p values of diet simple effects given each level of germ ******;

			data dietdiff;
				set dietdiff;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc sort data=dietdiff;
				by fig_id 'germ status'n;
			run;

			proc transpose data=dietdiff out=wide_dietdiff;
				by fig_id 'germ status'n;
				* Keeps this variable as a single row per subject;
				id diet;
				* Uses values in this column to name the new columns;
				var ProbtDiff;
				* The actual data values you want to move into columns;
			run;

			data pvalue_diet;
				length fig_id $25;
				retain fig_id variable;
				set wide_dietdiff;
				length variable $32;
				drop _name_ _label_ Palm;
				variable="%sysfunc(dequote(&var))";
				format borage Echium Fish HLS HOS best10.;
				*proc print;
			run;

			****** Select p values of germ and diet multiple comparisons ******;

			data germdietdiff;
				set germdietdiff;

				if ProbtDiff<0 then
					delete;
				fig_id="%sysfunc(dequote(&fig))";
			run;

			proc sort data=germdietdiff;
				by fig_id;
			run;

			proc transpose data=germdietdiff out=wide_germdietdiff;
				by fig_id;
				* Keeps this variable as a single row per subject;
				id 'germ status'n diet;
				* Uses values in this column to name the new columns;
				var ProbtDiff;
				* The actual data values you want to move into columns;
			run;

			data pvalue_germdiet;
				length fig_id $25;
				retain fig_id variable;
				set wide_germdietdiff;
				length variable $32;
				drop _name_ _label_;
				variable="%sysfunc(dequote(&var))";
				format GFBorage GFEchium GFFish GFHLS GFHOS GFPalm MPFBorage MPFEchium 
					MPFFish MPFHLS MPFHOS best10.;
				*proc print;
			run;
            
            proc delete data=wide_effect wide_overall pvalue_overall wide_germdiff wide_dietdiff wide_germdietdiff;
            run;
            
		%end;
		
*		options notes source mprint;
		ods select all;
		ods graphics on;
		
%mend one_variable;

/*%one_variable('Fig 1 D'n, 'Liver Weight (g)'n);

%one_variable('FIGURE 1 E-J'n, 'Plasma Active GLP-1'n);

%one_variable('Fig 6 B-L revised'n, 'Liver 12-Hydroxystearic Acid'n);

%one_variable('Fig S1 revised'n, 'Diet TG 16:0_20:4_20:5 minor'n);*/


************************************************************************************************************
*********************** Macro pipeline for pvalue with all variables in each figure ************************
************************************************************************************************************;

/*%macro all_variable(fid);
%local num_fid i this_id fig;
%let num_fid=%sysfunc(countw(&fid, %str(|)));

%do i=1 %to &num_fid;
%let this_id=%scan(&fid, &i, %str(|));

%if &i=1 %then
%let fig=&this_id;
%else
%do;
%one_variable(&fig, &this_id);
%let clean_fig= %sysfunc(substr(%sysfunc(dequote(&fig)), 1, 6));

%if &i=2 %then
%do;

data pvalue_normality1;
set pvalue_normality;
run;

%if &clean_fig=%str(Fig 6 ) or &clean_fig=Fig S1 %then
%do;

data pvalueFig6_FigS1_effect1;
set pvalueFig6_FigS1_effect;
run;

%if &clean_fig=%str(Fig 6 ) %then
%do;

data pvalueFig6_diet1;
set pvalueFig6_diet;
run;

%end;
%else %if &clean_fig=Fig S1 %then
%do;

data pvalueFigS1_diet1;
set pvalueFigS1_diet;
run;

%end;
%end;
%else
%do;

data pvalue_effect1;
set pvalue_effect;
run;

data pvalue_germ1;
set pvalue_germ;
run;

data pvalue_diet1;
set pvalue_diet;
run;

data pvalue_germdiet1;
set pvalue_germdiet;
run;

%end;
%end;
%else
%do;

data pvalue_normality1;
set pvalue_normality1 pvalue_normality;
run;

%if &clean_fig=%str(Fig 6 ) or &clean_fig=Fig S1 %then
%do;

data pvalueFig6_FigS1_effect1;
set pvalueFig6_FigS1_effect1 pvalueFig6_FigS1_effect;
run;

%if &clean_fig=%str(Fig 6 ) %then
%do;

data pvalueFig6_diet1;
set pvalueFig6_diet1 pvalueFig6_diet;
run;

%end;
%else %if &clean_fig=Fig S1 %then
%do;

data pvalueFigS1_diet1;
set pvalueFigS1_diet1 pvalueFigS1_diet;
run;

%end;
%end;
%else
%do;

data pvalue_effect1;
set pvalue_effect1 pvalue_effect;
run;

data pvalue_germ1;
set pvalue_germ1 pvalue_germ;
run;

data pvalue_diet1;
set pvalue_diet1 pvalue_diet;
run;

data pvalue_germdiet1;
set pvalue_germdiet1 pvalue_germdiet;
run;

%end;
%end;
%end;
%end;

proc print data=pvalue_normality1;
%if &clean_fig=%str(Fig 6 ) or &clean_fig=Fig S1 %then
%do;

proc print data=pvalueFig6_FigS1_effect1;
%if &clean_fig=%str(Fig 6 ) %then
%do;

proc print data=pvalueFig6_diet1;
%end;
%else %if &clean_fig=Fig S1 %then
%do;

proc print data=pvalueFigS1_diet1;
%end;
%end;
%else
%do;

proc print data=pvalue_effect1;
proc print data=pvalue_germ1;
proc print data=pvalue_diet1;
proc print data=pvalue_germdiet1;
%end;
%mend all_variable;

%all_variable('Fig 1 D'n|'Liver Weight (g)'n);
%all_variable('Fig 1 E-J'n|'Plasma Active GLP-1'n|'Plasma Glucagon'n|'Plasma PYY'n|'Plasma Insulin'n|'Plasma MCP-1'n|'Plasma Leptin'n);
%all_variable('Fig 1 K-L revised'n|'gWAT Weight'n|'Fast Glucose'n);
%all_variable('Fig 1 M revised'n|'Steatosis %'n);
%all_variable('Fig 3 A-H revised'n|'Plasma FA 22:5'n|'Plasma TG 18:1_20:4_22:6'n|'Plasma TG 16:0_18:1_18:1 isomer'n|'Plasma KDdiA-PC'n|
'Plasma Cer d 18:1_24:1'n|'Plasma Cer.d.18:1-22:0'n|'Plasma Acylcarnitine 16:1'n|'Plasma Cholesterol'n);
%all_variable('Fig 3 I-P'n|'Liver PE 38:4'n|'Liver PE 38:5'n|'Liver PC38:5 IsomerA'n|'Liver PC40:6'n|'Liver PI 38:3'n|
'Liver Cer d42-2'n|'Liver LPE-N  36:4'n|'Liver TG O57:3'n);
%all_variable('Fig 4 A-H'n|'LIVER aHYA'n|'LIVER aKETOD'n|'LIVER CLA1'n|'LIVER HYE'n|'LIVER gKETOB'n|'LIVER gHYA'n|'LIVER HYARA'n|'LIVER 12-OH'n);
%all_variable('Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n|'LIVER Protectin dx'n|'LIVER LIPOXIN'n|'LIVER RESOLVIN D5'n|'LIVER 20-HETE'n|
'LIVER N-DGLA EA'n|'LIVER 14,15 DIHETE'n|'LIVER 12,13 DIHOME'n);

%all_variable('Fig 6 B-L revised'n|'Liver Ricinoleic Acid'n|'Liver 12-Hydroxystearic Acid'n|'Liver 12-Oxostearic Acid'n|'Body Weight (g)'n|
'Liver Weight (g)'n|'gWAT Weight (g)'n|'Necropsy Blood Glucose (mg/dL)'n|'Liver TG 14:0_16:0_16:0'n|'Liver PC 20:1_14:1'n|
'Liver 20-Hydroxydocosa-hexaenoic'n|'Liver Cer d38:1'n);
%all_variable('Fig 6 P-S revised'n|'Rsph9'n|'Ciart'n|'Slc01a1'n|'Cyp4a14'n);
%all_variable('Fig S1 revised'n|'Diet TG 12:0_12:0_14:0'n|'Diet  TG 18:1_18:1_18:1'n|'Diet TG 18:1_18:2_18:2'n|'Diet TG 18:2_20:3_22:0'n|
'Diet  TG 18:3_18:3_18:3'n|'Diet TG 12:0_18:2_22:6'n|'Diet TG 16:0_20:4_20:5 minor'n);

%all_variable('Fig S2 revised'n|'Bacteroidetes'n|'Deferribacteres'n|'Verrucomicrobia'n|'Desulfobacterota'n|'Delta and gammaproteobacteria'n|'Firmicutes'n);
%all_variable('Fig S3 revised done'n|'Hepatic Srebf1'n|'Hepatic ACC1'n|'Hepatic FASN'n|'Hepatic CD68'n|'gWAT SREBF1'n|'gWAT ACC1'n|'gWAT FASN'n|
'gWAT CD68'n|'COLON SREBF1'n|'COLON ACC1'n|'COLON FASN'n|'COLON CD68'n);
%all_variable('Fig S5 A-H revised done'n|'PLASMA HYA'n|'PLASMA HYB'n|'PLASMA aHYA+aHYC'n|'PLASMA CALA2+CGLA1'n|'PLASMA CGLA2'n|
'PLASMA gHYD'n|'PLASMA 12-HAS'n|'PLASMA 12-Oxo'n);
%all_variable('Fig S5 I-P revised done'n|'ILEUM HYA'n|'ILEUM KETOA'n|'ILEUM aHYA+aHYC'n|'ILEUM aKETOD'n|'ILEUM gHYD'n|
'ILEUM gLA+CALA3'n|'ILEUM 12-HAS'n|'ILEUM 12-Oxo'n);
%all_variable('Fig S6 revised done'n|'Leukotriene B4'n|'11,12-DiHETre'n|'9.10-DiHOME'n|'14,15-DiHETre'n|'16,17-DeHDoPE'n|'12-HODE'n|'8,9-DiHETre'n|
'11,12-DiHETE'n|'5-HEPE'n|'10-HODE'n|'15-deoxy-PGJ2'n|'9-HOTE'n|'8,9-DiHETE'n|'17,18-DiHETE'n|'8,15-DiHETE'n|'Prostaglandin E2'n|
'N-C22:6n3_EA'n|'15-HETE'n|'13-KODE'n|'9-HODE'n);
%all_variable('Fig S9 revised done 'n|'Glycine C18:1'n|'Glycine C12:2'n|'Methionine C2:0'n|'Lysine C18:1'n|
'Tyrosine C 18:2'n|'Alanine/Sarcosine - C18:2'n);
*/


**********************************************************************************************************************
********************** Macro pipeline for pvalues with all variables together with all figures ***********************
**********************************************************************************************************************;

%macro all_figure(fad);
	%local num_fad a this_ad num_fid i this_id fig;
	ods exclude all;
	
	%let num_fad=%sysfunc(countw(&fad, %str(~)));

	%do a=1 %to &num_fad;
		%let this_ad=%scan(&fad, &a, %str(~));
		%let num_fid=%sysfunc(countw(&this_ad, %str(|)));

		%do i=1 %to &num_fid;
			%let this_id=%scan(&this_ad, &i, %str(|));

			%if &i=1 %then
				%let fig=&this_id;
			%else
				%do;
					%one_variable(&fig, &this_id);
					%let clean_fig= %sysfunc(substr(%sysfunc(dequote(&fig)), 1, 6));

					%if &i=2 %then
						%do;

							data pvalue_normality&a;
								set pvalue_normality;
							run;

							%if &clean_fig=%str(Fig 6 ) or &clean_fig=Fig S1 %then
								%do;

									data pvalueFig6_FigS1_effect&a;
										set pvalueFig6_FigS1_effect;
									run;

									%if &clean_fig=%str(Fig 6 ) %then
										%do;

											data pvalueFig6_diet&a;
												set pvalueFig6_diet;
											run;

										%end;
									%else %if &clean_fig=Fig S1 %then
										%do;

											data pvalueFigS1_diet&a;
												set pvalueFigS1_diet;
											run;

										%end;
								%end;
							%else
								%do;

									data pvalue_effect&a;
										set pvalue_effect;
									run;

									data pvalue_germ&a;
										set pvalue_germ;
									run;

									data pvalue_diet&a;
										set pvalue_diet;
									run;

									data pvalue_germdiet&a;
										set pvalue_germdiet;
									run;

								%end;
						%end;
					%else
						%do;

							data pvalue_normality&a;
								set pvalue_normality&a pvalue_normality;
							run;

							%if &clean_fig=%str(Fig 6 ) or &clean_fig=Fig S1 %then
								%do;

									data pvalueFig6_FigS1_effect&a;
										set pvalueFig6_FigS1_effect&a pvalueFig6_FigS1_effect;
									run;

									%if &clean_fig=%str(Fig 6 ) %then
										%do;

											data pvalueFig6_diet&a;
												set pvalueFig6_diet&a pvalueFig6_diet;
											run;

										%end;
									%else %if &clean_fig=Fig S1 %then
										%do;

											data pvalueFigS1_diet&a;
												set pvalueFigS1_diet&a pvalueFigS1_diet;
											run;

										%end;
								%end;
							%else
								%do;

									data pvalue_effect&a;
										set pvalue_effect&a pvalue_effect;
									run;

									data pvalue_germ&a;
										set pvalue_germ&a pvalue_germ;
									run;

									data pvalue_diet&a;
										set pvalue_diet&a pvalue_diet;
									run;

									data pvalue_germdiet&a;
										set pvalue_germdiet&a pvalue_germdiet;
									run;

								%end;
						%end;
				%end;
		%end;
	    
	    proc delete data=fig_1;
	    run;
	    
    %end;

	ods exclude none;
	
	data pvalue_normality_all;
		set pvalue_normality1-pvalue_normality&num_fad;

	proc print;
	run;
	
	proc delete data=pvalue_normality1-pvalue_normality&num_fad;
	run;

%mend all_figure;

*ods rtf file="/home/liux3/Manuscript/26/Mark Brown/Results/Interaction and simple effects.rtf";
%all_figure('Fig 1 D'n|'Liver Weight (g)'n~

	'Fig 1 E-J'n|'Plasma Active GLP-1'n|'Plasma Glucagon'n|'Plasma PYY'n|'Plasma Insulin'n|'Plasma MCP-1'n|'Plasma Leptin'n~
	
	'Fig 1 K-L revised'n|'gWAT Weight'n|'Fast Glucose'n~
	
	'Fig 1 M revised'n|'Steatosis %'n~
	
	'Fig 3 A-H revised'n|'Plasma FA 22:5'n|'Plasma TG 18:1_20:4_22:6'n|'Plasma TG 16:0_18:1_18:1 isomer'n|'Plasma KDdiA-PC'n|
              'Plasma Cer d 18:1_24:1'n|'Plasma Cer.d.18:1-22:0'n|'Plasma Acylcarnitine 16:1'n|'Plasma Cholesterol'n~
              
	'Fig 3 I-P'n|'Liver PE 38:4'n|'Liver PE 38:5'n|'Liver PC38:5 IsomerA'n|'Liver PC40:6'n|'Liver PI 38:3'n|
                  'Liver Cer d42-2'n|'Liver LPE-N  36:4'n|'Liver TG O57:3'n~
                  
	'Fig 4 A-H'n|'LIVER aHYA'n|'LIVER aKETOD'n|'LIVER CLA1'n|'LIVER HYE'n|'LIVER gKETOB'n|'LIVER gHYA'n|'LIVER HYARA'n|'LIVER 12-OH'n~
	
	'Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n|'LIVER Protectin dx'n|'LIVER LIPOXIN'n|'LIVER RESOLVIN D5'n|'LIVER 20-HETE'n|
                  'LIVER N-DGLA EA'n|'LIVER 14,15 DIHETE'n|'LIVER 12,13 DIHOME'n~
                  
                  
	'Fig 6 B-L revised'n|'Liver Ricinoleic Acid'n|'Liver 12-Hydroxystearic Acid'n|'Liver 12-Oxostearic Acid'n|'Body Weight (g)'n|
	               'Liver Weight (g)'n|'gWAT Weight (g)'n|'Necropsy Blood Glucose (mg/dL)'n|'Liver TG 14:0_16:0_16:0'n|'Liver PC 20:1_14:1'n|
	               'Liver 20-Hydroxydocosa-hexaenoic'n|'Liver Cer d38:1'n~
	               
	'Fig 6 P-S revised'n|'Rsph9'n|'Ciart'n|'Slc01a1'n|'Cyp4a14'n~
	
	'Fig S1 revised'n|'Diet TG 12:0_12:0_14:0'n|'Diet  TG 18:1_18:1_18:1'n|'Diet TG 18:1_18:2_18:2'n|'Diet TG 18:2_20:3_22:0'n|
	              'Diet  TG 18:3_18:3_18:3'n|'Diet TG 12:0_18:2_22:6'n|'Diet TG 16:0_20:4_20:5 minor'n~
	              
	              
	'Fig S2 revised'n|'Bacteroidetes'n|'Deferribacteres'n|'Verrucomicrobia'n|'Desulfobacterota'n|'Delta and gammaproteobacteria'n|'Firmicutes'n~
	
	'Fig S3 revised done'n|'Hepatic Srebf1'n|'Hepatic ACC1'n|'Hepatic FASN'n|'Hepatic CD68'n|'gWAT SREBF1'n|'gWAT ACC1'n|'gWAT FASN'n|
              'gWAT CD68'n|'COLON SREBF1'n|'COLON ACC1'n|'COLON FASN'n|'COLON CD68'n~
              
	'Fig S5 A-H revised done'n|'PLASMA HYA'n|'PLASMA HYB'n|'PLASMA aHYA+aHYC'n|'PLASMA CALA2+CGLA1'n|'PLASMA CGLA2'n|
		              'PLASMA gHYD'n|'PLASMA 12-HAS'n|'PLASMA 12-Oxo'n~
              
	'Fig S5 I-P revised done'n|'ILEUM HYA'n|'ILEUM KETOA'n|'ILEUM aHYA+aHYC'n|'ILEUM aKETOD'n|'ILEUM gHYD'n|
              'ILEUM gLA+CALA3'n|'ILEUM 12-HAS'n|'ILEUM 12-Oxo'n~
              
	'Fig S6 revised done'n|'Leukotriene B4'n|'11,12-DiHETre'n|'9.10-DiHOME'n|'14,15-DiHETre'n|'16,17-DeHDoPE'n|'12-HODE'n|'8,9-DiHETre'n|
              '11,12-DiHETE'n|'5-HEPE'n|'10-HODE'n|'15-deoxy-PGJ2'n|'9-HOTE'n|'8,9-DiHETE'n|'17,18-DiHETE'n|'8,15-DiHETE'n|'Prostaglandin E2'n|
		              'N-C22:6n3_EA'n|'15-HETE'n|'13-KODE'n|'9-HODE'n~
              
	'Fig S9 revised done 'n|'Glycine C18:1'n|'Glycine C12:2'n|'Methionine C2:0'n|'Lysine C18:1'n|
              'Tyrosine C 18:2'n|'Alanine/Sarcosine - C18:2'n);
              

********************************************************************************************************              
********************************************** Normality check *****************************************
********************************************************************************************************;

proc export data=pvalue_normality_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_normality_all.csv" 
		dbms=csv replace;
run;              


********************************************************************************************************
********************************** Main and interaction effects **************************************** 
********************************************************************************************************;

************* Main and interaction effects of germ and diet **********************;

data pvalue_effect_all;
	set pvalue_effect1-pvalue_effect8 pvalue_effect12-pvalue_effect17;

proc print;
run;

proc delete data=pvalue_effect1-pvalue_effect8 pvalue_effect12-pvalue_effect17;
run;

%macro adjp_effect (pvalue, fig, meff);
	%local num_fig i this_ig this_ff this_ffn;

	data &pvalue;
		set pvalue_effect_all;
		array target_cols 'GERM STATUS'n DIET 'GERM STATUS*DIET'n 'OVERALL MODEL'n;
		%let num_fig= %sysfunc(countw(&fig, %str(|)));

		%do i=1 %to &num_fig;
			%let this_ig=%scan(&fig, &i, %str(|));
			%let this_ff=%scan(&meff, &i, %str( ));
			%let this_ffn=%sysfunc(inputn(&this_ff, best12.));

			%if &i<=9 or &i>=13 %then
				%do;

					if fig_id="%sysfunc(dequote(&this_ig))" then
						do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				%end;
			%else %if &i=10 %then
				%do;

					if fig_id="%sysfunc(dequote(&this_ig))" and substr(variable, 1, 
						7)="Hepatic" then
							do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				%end;
			%else %if &i=11 %then
				%do;

					if fig_id="%sysfunc(dequote(&this_ig))" and substr(variable, 1, 4)="gWAT" 
						then
							do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				%end;
			%else %if &i=12 %then
				%do;

					if fig_id="%sysfunc(dequote(&this_ig))" and substr(variable, 1, 5)='COLON' 
						then
							do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				%end;
		%end;
		drop j;
	run;

proc print;

%mend adjp_effect;

%adjp_effect(pvalue_effect_all_adj, 'Fig 1 D'n|'Fig 1 E-J'n|'Fig 1 K-L revised'n|'Fig 1 M revised'n|'Fig 3 A-H revised'n|'Fig 3 I-P'n|'Fig 4 A-H'n|'Fig 4 I-P'n|
     'Fig S2 revised'n|'Fig S3 revised done'n|'Fig S3 revised done'n|'Fig S3 revised done'n|'Fig S5 A-H revised done'n|
     'Fig S5 I-P revised done'n|'Fig S6 revised done'n|'Fig S9 revised done 'n, 1 2 1 1 2 3 3 3 2 1 2 1 3 4 4 2);
     
proc export data=pvalue_effect_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_effect_all.csv" 
		dbms=csv replace;
run;

proc export data=pvalue_effect_all_adj 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_effect_all_adj.csv" 
		dbms=csv replace;
run; 

             
******************* Main effects of diet for Fig 6 and Fig S1 *********************;

data pvalueFig6_FigS1_effect_all;
	set pvalueFig6_FigS1_effect9-pvalueFig6_FigS1_effect11;

proc print;
run;

proc delete data=pvalueFig6_FigS1_effect9-pvalueFig6_FigS1_effect11;
run;

%macro adjpFig6_FigS1 (pvalue, fig, meff);
	%local num_fig i this_ig this_ff this_ffn;

	data &pvalue;
		set pvalueFig6_FigS1_effect_all;
		array target_cols DIET 'OVERALL MODEL'n;
		%let num_fig= %sysfunc(countw(&fig, %str(|)));

		%do i=1 %to &num_fig;
			%let this_ig=%scan(&fig, &i, %str(|));
			%let this_ff=%scan(&meff, &i, %str( ));
			%let this_ffn=%sysfunc(inputn(&this_ff, best12.));

					if fig_id="%sysfunc(dequote(&this_ig))" then
						do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				
		%end;
		drop j;
	run;
 
proc print;

%mend adjpFig6_FigS1;

%adjpFig6_FigS1(pvalueFig6_FigS1_effect_all_adj, 'Fig 6 B-L revised'n|'Fig 6 P-S revised'n|'Fig S1 revised'n, 2 2 3);

proc export data=pvalueFig6_FigS1_effect_all_adj 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalueFig6_FigS1_effect_all_adj.csv" 
		dbms=csv replace;
run; 
		

**************** Multiple comparisons of diet for Fig 6 **********************;

data pvalueFig6_diet_all;
	set pvalueFig6_diet9 pvalueFig6_diet10;

proc print;
run;

proc delete data=pvalueFig6_diet9 pvalueFig6_diet10;
run;

%macro adjpFig6_diet (pvalue, fig, meff);
	%local num_fig i this_ig this_ff this_ffn;

	data &pvalue;
		set pvalueFig6_diet_all;
		array target_cols '12-HSA'n '12-Oxo'n RA;
		%let num_fig= %sysfunc(countw(&fig, %str(|)));

		%do i=1 %to &num_fig;
			%let this_ig=%scan(&fig, &i, %str(|));
			%let this_ff=%scan(&meff, &i, %str( ));
			%let this_ffn=%sysfunc(inputn(&this_ff, best12.));

					if fig_id="%sysfunc(dequote(&this_ig))" then
						do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				
		%end;
		drop j;
	run;
 
proc print;

%mend adjpFig6_diet;

%adjpFig6_diet(pvalueFig6_diet_all_adj, 'Fig 6 B-L revised'n|'Fig 6 P-S revised'n, 2 2);

proc export data=pvalueFig6_diet_all_adj 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalueFig6_diet_all_adj.csv" 
		dbms=csv replace;
run;
 

**************** Multiple comparisons of diet for Fig S1 **********************;

data pvalueFigS1_diet_all;
	set pvalueFigS1_diet11;

proc print;
run;

%macro adjpFigS1_diet (pvalue, fig, meff);
	%local num_fig i this_ig this_ff this_ffn;

	data &pvalue;
		set pvalueFigS1_diet_all;
		array target_cols Borage Echium Fish HLS HOS;
		%let num_fig= %sysfunc(countw(&fig, %str(|)));

		%do i=1 %to &num_fig;
			%let this_ig=%scan(&fig, &i, %str(|));
			%let this_ff=%scan(&meff, &i, %str( ));
			%let this_ffn=%sysfunc(inputn(&this_ff, best12.));

					if fig_id="%sysfunc(dequote(&this_ig))" then
						do;
							Meff=&this_ffn;

							do j=1 to dim(target_cols);
								target_cols[j]=target_cols[j] * &this_ffn;

								if target_cols[j] >1 then
									target_cols[j]=0.9999;
							end;
						end;
				
		%end;
		drop j;
	run;
 
proc print;

%mend adjpFigS1_diet;

%adjpFigS1_diet(pvalueFigS1_diet_all_adj, 'Fig S1 revised'n, 3);

proc export data=pvalueFigS1_diet_all_adj 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalueFigS1_diet_all_adj.csv" 
		dbms=csv replace;
run;


***************** Simple effects of germ given each diet ***********************;

data pvalue_germ_all;
	set pvalue_germ1-pvalue_germ8 pvalue_germ12-pvalue_germ17;

proc print;
run;

proc export data=pvalue_germ_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_germ_all.csv" 
		dbms=csv replace;
run;

***************** Simple effects of diet given each germ ***********************;

data pvalue_diet_all;
	set pvalue_diet1-pvalue_diet8 pvalue_diet12-pvalue_diet17;

proc print;
run;

proc export data=pvalue_diet_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_diet_all.csv" 
		dbms=csv replace;
run;

********************* Multiple comparisons of germ and diet combinations **************************;

data pvalue_germdiet_all;
	set pvalue_germdiet1-pvalue_germdiet8 pvalue_germdiet12-pvalue_germdiet17;

proc print;
run;

proc export data=pvalue_germdiet_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_germdiet_all.csv" 
		dbms=csv replace;
run;

		
		
****** Calculated adjusted p vales directly without a macro pipeline;

	/*data pvalue_effect_adj;
	set pvalue_effect_all;
	array target_cols 'GERM STATUS'n DIET 'GERM STATUS*DIET'n 'OVERALL MODEL'n;
	
	if fig_id='FIGURE 1 E-J' then do;
	Meff=3.52;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig S1 A-L' and substr(variable, 1, 7)='Hepatic' then do;
	Meff=2.67;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	else if fig_id='Fig S1 A-L' and substr(variable, 1, 4)='gWAT' then do;
	Meff=2.57;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	else if fig_id='Fig S1 A-L' and substr(variable, 1, 5)='COLON' then do;
	Meff=1.82;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig 4 A-H' then do;
	Meff=4.41;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig 4 I-P' then do;
	Meff=4.85;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig S4 A-H' then do;
	Meff=5.68;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig S4 I-P' then do;
	Meff=5.04;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='FigS5 liver oxylipins' then do;
	Meff=7.97;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig3 plasma FA' then do;
	Meff=6.41;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig3 liver lipidomics' then do;
	Meff=8.77;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	drop i;
	proc print;
	run;
	
	
	data pvalue_effect_adj2;
	set pvalue_effect_all;
	array target_cols 'GERM STATUS'n DIET 'GERM STATUS*DIET'n 'OVERALL MODEL'n;
	
	if fig_id='FIGURE 1 E-J' then do;
	Meff=2;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig S1 A-L' and substr(variable, 1, 7)='Hepatic' then do;
	Meff=1;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	else if fig_id='Fig S1 A-L' and substr(variable, 1, 4)='gWAT' then do;
	Meff=2;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	else if fig_id='Fig S1 A-L' and substr(variable, 1, 5)='COLON' then do;
	Meff=1;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig 4 A-H' then do;
	Meff=3;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig 4 I-P' then do;
	Meff=3;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig S4 A-H' then do;
	Meff=3;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig S4 I-P' then do;
	Meff=4;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='FigS5 liver oxylipins' then do;
	Meff=4;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig3 plasma FA' then do;
	Meff=4;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	if fig_id='Fig3 liver lipidomics' then do;
	Meff=5;
	do i = 1 to dim(target_cols);
	target_cols[i] = target_cols[i] * Meff;
	if target_cols[i] >1 then target_cols[i] =0.9999;
	end;
	end;
	drop i;
	proc print;
	run;
	*/
	*ods rtf close;



****** Export SAS datasets into split excel sheets;
*** Split SAS data;

/*data pvalue_normality_main;
set pvalue_normality_all;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_normality_supp;
set pvalue_normality_all;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;

data pvalue_effect_main;
set pvalue_effect_all;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_effect_supp;
set pvalue_effect_all;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;

data pvalue_germ_main;
set pvalue_germ_all;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_germ_supp;
set pvalue_germ_all;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;

data pvalue_diet_main;
set pvalue_diet_all;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_diet_supp;
set pvalue_diet_all;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;

data pvalue_germdiet_main;
set pvalue_germdiet_all;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_germdiet_supp;
set pvalue_germdiet_all;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;

data pvalue_effect_adj_main;
set pvalue_effect_adj;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_effect_adj_supp;
set pvalue_effect_adj;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;

data pvalue_effect_adj2_main;
set pvalue_effect_adj2;
if fid_id in ("FIGURE 1 E-J", "Fig 4 A-H", "Fig 4 I-P", "Fig3 plasma FA", "Fig3 liver lipidomics");
run;

data pvalue_effect_adj2_supp;
set pvalue_effect_adj2;
if fid_id in ("Fig S1 A-L", "Fig S4 A-H", "Fig S4 I-P", "FigS5 liver oxylipins");
run;


 *** Export data into excell sheets;
proc export data=pvalue_normality_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_normality_all.csv"
dbms=csv
replace;
sheet="pvalue_normality_main";
run;
proc export data=pvalue_normality_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_normality_all.csv"
dbms=csv;
sheet="pvalue_normality_supp";
run;


proc export data=pvalue_effect_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_effect_all.csv"
dbms=csv
replace;
sheet="pvalue_effect_main";
run;
proc export data=pvalue_effect_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_effect_all.csv"
dbms=csv;
sheet="pvalue_effect_supp";
run;


proc export data=pvalue_germ_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_germ_all.csv"
dbms=csv
replace;
sheet="pvalue_germ_main";
run;
proc export data=pvalue_germ_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_germ_all.csv"
dbms=csv;
sheet="pvalue_germ_supp";
run;


proc export data=pvalue_diet_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_diet_all.csv"
dbms=csv
replace;
sheet="value_diet_main";
run;
proc export data=pvalue_diet_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_diet_all.csv"
dbms=csv;
sheet="value_diet_supp";
run;


proc export data=pvalue_germdiet_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_germdiet_all.csv"
dbms=csv
replace;
sheet="pvalue_germdiet_main";
run;
proc export data=pvalue_germdiet_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_germdiet_all.csv"
dbms=csv;
sheet="pvalue_germdiet_supp";
run;


proc export data=pvalue_effect_adj_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_effect_adj.csv"
dbms=csv
replace;
sheet="pvalue_effect_adj_main";
run;
proc export data=pvalue_effect_adj_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_effect_adj.csv"
dbms=csv;
sheet="pvalue_effect_adj_supp";
run;


proc export data=pvalue_effect_adj2_main
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_effect_adj2.csv"
dbms=csv
replace;
sheet="pvalue_effect_adj2_main";
run;
proc export data=pvalue_effect_adj2_supp
outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_effect_adj2.csv"
dbms=csv;
sheet="pvalue_effect_adj2_supp";
run;*/